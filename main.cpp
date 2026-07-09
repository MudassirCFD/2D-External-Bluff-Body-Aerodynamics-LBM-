/*
 * 2D External Bluff Body Aerodynamics
 * D2Q9 TRT Lattice Boltzmann solver
 *
 * Cases
 *   Circle
 *   Square
 *   Droplet
 *
 * This file is the main C++ solver used for the wake shape study.
 * The basic LBM method follows standard D2Q9 examples and published
 * Lattice Boltzmann references. I adapted the setup for this shape
 * comparison, added the three body masks, force output, smoke field,
 * sponge damping, time snapshots and the case control used in the repo.
 *
 * Main physics
 *   Re = U_inf D / nu
 *   D2Q9 lattice
 *   TRT collision
 *   uniform inlet
 *   pressure outlet with rho = 1
 *   slip style top and bottom treatment
 *   far field sponge to reduce reflection and density drift
 *   momentum exchange force calculation
 *   passive scalar smoke for wake visualisation
 *
 * Output folders
 *   circle_snapshots
 *   square_snapshots
 *   droplet_snapshots
 *
 * Snapshot columns
 *   i,j,x,y,rho,ux,uy,curl,p,Cp,solid,smoke
 *
 * Force columns
 *   step,t_star,Cd,Cl,Fx_raw,Fy_raw
 *
 * Author
 *   Saiyed Mudassir
 *
 * Build
 *   Visual Studio x64 Release
 *   g++ -O3 -march=native -std=c++17 -o lbm_aero main.cpp
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <direct.h>
#define MAKE_DIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define MAKE_DIR(path) mkdir(path, 0755)
#endif


// ============================================================================
// Solver constants
// ============================================================================

static constexpr int NX = 1000;
static constexpr int NY = 320;

static constexpr int D = 40;
static constexpr int R = D / 2;
static constexpr int CX = 180;
static constexpr int CY = NY / 2;

static constexpr double U_INF = 0.06;
static constexpr double U_REF = U_INF;
static constexpr double RHO_REF = 1.0;
static constexpr double RE = 200.0;

static constexpr int NWARM = 30000;
static constexpr int NRECORD = 60000;
static constexpr int SNAP_EVERY = 200;
static constexpr int FORCE_EVERY = 10;

static constexpr int RAMP_STEPS = 8000;
static constexpr int SEED_STEPS = 12000;

static constexpr bool RUN_CIRCLE  = true;
static constexpr bool RUN_SQUARE  = true;
static constexpr bool RUN_DROPLET = true;

// Passive scalar smoke.
static constexpr double SMOKE_TAU = 0.58;
static constexpr int SMOKE_SOURCE_X = CX - R - 6;
static constexpr double SMOKE_SOURCE_SIGMA = 0.10 * static_cast<double>(D);

// Numerical farfield sponge.
static constexpr int SPONGE_X_WIDTH = 160;
static constexpr int SPONGE_Y_WIDTH = 36;
static constexpr double SPONGE_STRENGTH = 0.018;


// ============================================================================
// LBM constants
// ============================================================================

static constexpr int Q = 9;

static constexpr double PI = 3.14159265358979323846;
static constexpr double CS2 = 1.0 / 3.0;
static constexpr double CS = 0.57735026918962576451;

static constexpr double NU = U_REF * static_cast<double>(D) / RE;

static constexpr double TAU_PLUS = 0.5 + 3.0 * NU;
static constexpr double OMEGA_PLUS = 1.0 / TAU_PLUS;

// TRT magic parameter.
static constexpr double TRT_LAMBDA = 3.0 / 16.0;
static constexpr double TAU_MINUS = 0.5 + TRT_LAMBDA / (TAU_PLUS - 0.5);
static constexpr double OMEGA_MINUS = 1.0 / TAU_MINUS;

static constexpr double SMOKE_OMEGA = 1.0 / SMOKE_TAU;
static constexpr double MA = U_INF / CS;

static constexpr std::array<int, Q> EX = { {
    0,  1,  0, -1,  0,  1, -1, -1,  1
} };

static constexpr std::array<int, Q> EY = { {
    0,  0,  1,  0, -1,  1,  1, -1, -1
} };

static constexpr std::array<int, Q> OPP = { {
    0, 3, 4, 1, 2, 7, 8, 5, 6
} };

// Mirror across horizontal boundaries for slip.
static constexpr std::array<int, Q> MIRROR_Y = { {
    0, 1, 4, 3, 2, 8, 7, 6, 5
} };

static constexpr std::array<double, Q> W = { {
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0
} };


// ============================================================================
// Types and indexing
// ============================================================================

using Field = std::vector<double>;
using Mask = std::vector<std::uint8_t>;

inline std::size_t cell_id(const int i, const int j)
{
    return static_cast<std::size_t>(i) * static_cast<std::size_t>(NY)
        + static_cast<std::size_t>(j);
}

inline std::size_t dist_id(const int i, const int j, const int q)
{
    return (
        static_cast<std::size_t>(i) * static_cast<std::size_t>(NY)
        + static_cast<std::size_t>(j)
        ) * static_cast<std::size_t>(Q)
        + static_cast<std::size_t>(q);
}

inline bool finite_number(const double x)
{
    return std::isfinite(x) != 0;
}


// ============================================================================
// Smooth startup
// ============================================================================

double ramp_factor(const int step)
{
    if (step >= RAMP_STEPS) {
        return 1.0;
    }

    const double s = static_cast<double>(step) / static_cast<double>(RAMP_STEPS);

    return 0.5 * (1.0 - std::cos(PI * s));
}

double current_u_inf(const int step)
{
    return U_INF * ramp_factor(step);
}


// ============================================================================
// Equilibrium functions
// ============================================================================

inline double feq(const int q, const double rho, const double ux, const double uy)
{
    const double eu = static_cast<double>(EX[q]) * ux
        + static_cast<double>(EY[q]) * uy;

    const double u2 = ux * ux + uy * uy;

    return W[q] * rho * (
        1.0
        + eu / CS2
        + eu * eu / (2.0 * CS2 * CS2)
        - u2 / (2.0 * CS2)
        );
}

inline double geq(const int q, const double phi, const double ux, const double uy)
{
    const double eu = static_cast<double>(EX[q]) * ux
        + static_cast<double>(EY[q]) * uy;

    return W[q] * phi * (1.0 + eu / CS2);
}


// ============================================================================
// Geometry
// ============================================================================

enum class Shape
{
    Circle,
    Square,
    Droplet
};

std::string shape_name(const Shape shape)
{
    if (shape == Shape::Circle) {
        return "circle";
    }

    if (shape == Shape::Square) {
        return "square";
    }

    return "droplet";
}

std::string shape_label(const Shape shape)
{
    if (shape == Shape::Circle) {
        return "CIRCLE";
    }

    if (shape == Shape::Square) {
        return "SQUARE";
    }

    return "DROPLET";
}

Mask build_body_mask(const Shape shape)
{
    Mask body(static_cast<std::size_t>(NX) * static_cast<std::size_t>(NY), 0);

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const int dx = i - CX;
            const int dy = j - CY;

            bool inside = false;

            if (shape == Shape::Circle) {
                inside = (dx * dx + dy * dy) <= R * R;
            }
            else if (shape == Shape::Square) {
                inside = (std::abs(dx) <= R) && (std::abs(dy) <= R);
            }
            else {
                // Rigid teardrop: rounded upstream nose, tapered downstream tail.
                // Flow direction is left to right.
                if (dx <= 0) {
                    inside = (dx * dx + dy * dy) <= R * R;
                }
                else {
                    const double tail_length = 3.0 * static_cast<double>(R);

                    if (static_cast<double>(dx) <= tail_length) {
                        const double s = static_cast<double>(dx) / tail_length;
                        const double half_height =
                            static_cast<double>(R) * std::pow(1.0 - s, 0.70);

                        inside = std::abs(static_cast<double>(dy)) <= half_height;
                    }
                }
            }

            body[cell_id(i, j)] = inside ? 1u : 0u;
        }
    }

    return body;
}


// ============================================================================
// Macroscopic fields
// ============================================================================

void compute_macros(
    const Field& f,
    const Mask& body,
    Field& rho,
    Field& ux,
    Field& uy
)
{
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);

            if (body[id] != 0u) {
                rho[id] = RHO_REF;
                ux[id] = 0.0;
                uy[id] = 0.0;
                continue;
            }

            const std::size_t base = id * static_cast<std::size_t>(Q);

            double r = 0.0;
            double px = 0.0;
            double py = 0.0;

            for (int q = 0; q < Q; ++q) {
                const double fq = f[base + static_cast<std::size_t>(q)];

                r += fq;
                px += static_cast<double>(EX[q]) * fq;
                py += static_cast<double>(EY[q]) * fq;
            }

            rho[id] = r;

            if (r > 1.0e-14) {
                ux[id] = px / r;
                uy[id] = py / r;
            }
            else {
                ux[id] = 0.0;
                uy[id] = 0.0;
            }
        }
    }
}

void compute_smoke(
    const Field& g,
    const Mask& body,
    Field& smoke
)
{
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);

            if (body[id] != 0u) {
                smoke[id] = 0.0;
                continue;
            }

            const std::size_t base = id * static_cast<std::size_t>(Q);

            double phi = 0.0;

            for (int q = 0; q < Q; ++q) {
                phi += g[base + static_cast<std::size_t>(q)];
            }

            smoke[id] = phi > 0.0 ? phi : 0.0;
        }
    }
}


// ============================================================================
// Controlled wake seeding
// ============================================================================

void apply_wake_seed(
    Field& ux,
    Field& uy,
    const Mask& body,
    const int step
)
{
    if (step >= SEED_STEPS) {
        return;
    }

    const double amp =
        0.010 * U_INF * (1.0 - static_cast<double>(step) / static_cast<double>(SEED_STEPS));

    const double omega = 2.0 * PI / 900.0;

    for (int i = CX; i < CX + 8 * D && i < NX - 1; ++i) {
        for (int j = CY - 3 * D; j <= CY + 3 * D; ++j) {
            if (j <= 1 || j >= NY - 2) {
                continue;
            }

            const std::size_t id = cell_id(i, j);

            if (body[id] != 0u) {
                continue;
            }

            const double dx = static_cast<double>(i - (CX + R)) / static_cast<double>(D);
            const double dy = static_cast<double>(j - CY) / static_cast<double>(D);

            const double envelope =
                std::exp(-0.5 * (dx * dx / 5.0 + dy * dy / 1.2));

            const double lateral =
                std::sin(omega * static_cast<double>(step))
                + 0.35 * std::sin(0.63 * omega * static_cast<double>(step) + 1.3);

            uy[id] += amp * envelope * lateral;

            // Tiny streamwise relaxation near the seeded region.
            ux[id] = std::max(0.0, ux[id]);
        }
    }
}


// ============================================================================
// Collision
// ============================================================================

void collide_fluid_trt(
    const Field& f,
    Field& f_post,
    const Field& rho,
    const Field& ux,
    const Field& uy,
    const Mask& body
)
{
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);
            const std::size_t base = id * static_cast<std::size_t>(Q);

            if (body[id] != 0u) {
                for (int q = 0; q < Q; ++q) {
                    f_post[base + static_cast<std::size_t>(q)] =
                        feq(q, RHO_REF, 0.0, 0.0);
                }
                continue;
            }

            const double r = rho[id];
            const double u = ux[id];
            const double v = uy[id];

            for (int q = 0; q < Q; ++q) {
                const int oq = OPP[q];

                const double fq = f[base + static_cast<std::size_t>(q)];
                const double foq = f[base + static_cast<std::size_t>(oq)];

                const double feq_q = feq(q, r, u, v);
                const double feq_oq = feq(oq, r, u, v);

                const double f_plus = 0.5 * (fq + foq);
                const double f_minus = 0.5 * (fq - foq);
                const double eq_plus = 0.5 * (feq_q + feq_oq);
                const double eq_minus = 0.5 * (feq_q - feq_oq);

                f_post[base + static_cast<std::size_t>(q)] =
                    fq
                    - OMEGA_PLUS * (f_plus - eq_plus)
                    - OMEGA_MINUS * (f_minus - eq_minus);
            }
        }
    }
}

void collide_smoke_bgk(
    const Field& g,
    Field& g_post,
    const Field& smoke,
    const Field& ux,
    const Field& uy,
    const Mask& body
)
{
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);
            const std::size_t base = id * static_cast<std::size_t>(Q);

            if (body[id] != 0u) {
                for (int q = 0; q < Q; ++q) {
                    g_post[base + static_cast<std::size_t>(q)] = 0.0;
                }
                continue;
            }

            const double phi = smoke[id];

            for (int q = 0; q < Q; ++q) {
                const std::size_t k = base + static_cast<std::size_t>(q);
                const double eq = geq(q, phi, ux[id], uy[id]);

                g_post[k] = g[k] - SMOKE_OMEGA * (g[k] - eq);
            }
        }
    }
}


// ============================================================================
// Streaming and bounce-back
// ============================================================================

struct Forces
{
    double Fx;
    double Fy;

    Forces() : Fx(0.0), Fy(0.0) {}
};

Forces stream_fluid(
    const Field& f_post,
    Field& f,
    const Mask& body
)
{
    Forces force;

    std::fill(f.begin(), f.end(), 0.0);

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);
            const std::size_t base = id * static_cast<std::size_t>(Q);

            if (body[id] != 0u) {
                for (int q = 0; q < Q; ++q) {
                    f[base + static_cast<std::size_t>(q)] =
                        feq(q, RHO_REF, 0.0, 0.0);
                }
                continue;
            }

            for (int q = 0; q < Q; ++q) {
                if (q == 0) {
                    f[base] = f_post[base];
                    continue;
                }

                const int si = i - EX[q];
                const int sj = j - EY[q];

                if (si < 0 || si >= NX) {
                    continue;
                }

                if (sj < 0 || sj >= NY) {
                    const int mq = MIRROR_Y[q];

                    f[base + static_cast<std::size_t>(q)] =
                        f_post[base + static_cast<std::size_t>(mq)];

                    continue;
                }

                const std::size_t sid = cell_id(si, sj);

                if (body[sid] != 0u) {
                    const int oq = OPP[q];
                    const double fb =
                        f_post[base + static_cast<std::size_t>(oq)];

                    f[base + static_cast<std::size_t>(q)] = fb;

                    force.Fx += -2.0 * fb * static_cast<double>(EX[q]);
                    force.Fy += -2.0 * fb * static_cast<double>(EY[q]);

                    continue;
                }

                f[base + static_cast<std::size_t>(q)] =
                    f_post[dist_id(si, sj, q)];
            }
        }
    }

    return force;
}

void stream_smoke(
    const Field& g_post,
    Field& g,
    const Mask& body
)
{
    std::fill(g.begin(), g.end(), 0.0);

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);
            const std::size_t base = id * static_cast<std::size_t>(Q);

            if (body[id] != 0u) {
                for (int q = 0; q < Q; ++q) {
                    g[base + static_cast<std::size_t>(q)] = 0.0;
                }
                continue;
            }

            for (int q = 0; q < Q; ++q) {
                if (q == 0) {
                    g[base] = g_post[base];
                    continue;
                }

                const int si = i - EX[q];
                const int sj = j - EY[q];

                if (si < 0 || si >= NX) {
                    continue;
                }

                if (sj < 0 || sj >= NY) {
                    const int mq = MIRROR_Y[q];

                    g[base + static_cast<std::size_t>(q)] =
                        g_post[base + static_cast<std::size_t>(mq)];

                    continue;
                }

                const std::size_t sid = cell_id(si, sj);

                if (body[sid] != 0u) {
                    const int oq = OPP[q];

                    g[base + static_cast<std::size_t>(q)] =
                        g_post[base + static_cast<std::size_t>(oq)];

                    continue;
                }

                g[base + static_cast<std::size_t>(q)] =
                    g_post[dist_id(si, sj, q)];
            }
        }
    }
}


// ============================================================================
// Boundary conditions
// ============================================================================

void apply_inlet_zou_he(Field& f, const int step)
{
    const int i = 0;
    const double u0 = current_u_inf(step);

    for (int j = 1; j < NY - 1; ++j) {
        const std::size_t base = cell_id(i, j) * static_cast<std::size_t>(Q);

        const double f0 = f[base + 0u];
        const double f2 = f[base + 2u];
        const double f4 = f[base + 4u];
        const double f3 = f[base + 3u];
        const double f6 = f[base + 6u];
        const double f7 = f[base + 7u];

        const double rho =
            (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - u0);

        f[base + 1u] = f3 + (2.0 / 3.0) * rho * u0;

        f[base + 5u] =
            f7
            - 0.5 * (f2 - f4)
            + (1.0 / 6.0) * rho * u0;

        f[base + 8u] =
            f6
            + 0.5 * (f2 - f4)
            + (1.0 / 6.0) * rho * u0;
    }
}

void apply_outlet_pressure_zou_he(Field& f)
{
    const int i = NX - 1;
    const double rho = RHO_REF;

    for (int j = 1; j < NY - 1; ++j) {
        const std::size_t base = cell_id(i, j) * static_cast<std::size_t>(Q);

        const double f0 = f[base + 0u];
        const double f1 = f[base + 1u];
        const double f2 = f[base + 2u];
        const double f4 = f[base + 4u];
        const double f5 = f[base + 5u];
        const double f8 = f[base + 8u];

        double ux_out =
            -1.0 + (f0 + f2 + f4 + 2.0 * (f1 + f5 + f8)) / rho;

        if (ux_out > 0.20) {
            ux_out = 0.20;
        }

        if (ux_out < -0.05) {
            ux_out = -0.05;
        }

        f[base + 3u] = f1 - (2.0 / 3.0) * rho * ux_out;

        f[base + 6u] =
            f8
            - 0.5 * (f2 - f4)
            - (1.0 / 6.0) * rho * ux_out;

        f[base + 7u] =
            f5
            + 0.5 * (f2 - f4)
            - (1.0 / 6.0) * rho * ux_out;
    }
}

void apply_top_bottom_slip(Field& f)
{
    for (int i = 0; i < NX; ++i) {
        f[dist_id(i, 0, 2)] = f[dist_id(i, 0, 4)];
        f[dist_id(i, 0, 5)] = f[dist_id(i, 0, 8)];
        f[dist_id(i, 0, 6)] = f[dist_id(i, 0, 7)];

        f[dist_id(i, NY - 1, 4)] = f[dist_id(i, NY - 1, 2)];
        f[dist_id(i, NY - 1, 7)] = f[dist_id(i, NY - 1, 6)];
        f[dist_id(i, NY - 1, 8)] = f[dist_id(i, NY - 1, 5)];
    }
}

void reset_body_distributions(Field& f, const Mask& body)
{
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);

            if (body[id] == 0u) {
                continue;
            }

            const std::size_t base = id * static_cast<std::size_t>(Q);

            for (int q = 0; q < Q; ++q) {
                f[base + static_cast<std::size_t>(q)] =
                    feq(q, RHO_REF, 0.0, 0.0);
            }
        }
    }
}

double sponge_strength(const int i, const int j)
{
    double sx = 0.0;
    double sy = 0.0;

    const int x_start = NX - SPONGE_X_WIDTH;

    if (i >= x_start) {
        const double s =
            static_cast<double>(i - x_start) / static_cast<double>(SPONGE_X_WIDTH);

        sx = s * s;
    }

    if (j < SPONGE_Y_WIDTH) {
        const double s =
            static_cast<double>(SPONGE_Y_WIDTH - j) / static_cast<double>(SPONGE_Y_WIDTH);

        sy = s * s;
    }
    else if (j > NY - 1 - SPONGE_Y_WIDTH) {
        const double s =
            static_cast<double>(j - (NY - 1 - SPONGE_Y_WIDTH))
            / static_cast<double>(SPONGE_Y_WIDTH);

        sy = s * s;
    }

    const double s = std::max(sx, sy);

    return SPONGE_STRENGTH * s;
}

void apply_farfield_sponge(Field& f, const Mask& body, const int step)
{
    const double ufar = current_u_inf(step);

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            const std::size_t id = cell_id(i, j);

            if (body[id] != 0u) {
                continue;
            }

            const double sigma = sponge_strength(i, j);

            if (sigma <= 0.0) {
                continue;
            }

            const std::size_t base = id * static_cast<std::size_t>(Q);

            for (int q = 0; q < Q; ++q) {
                const double target = feq(q, RHO_REF, ufar, 0.0);

                f[base + static_cast<std::size_t>(q)] =
                    (1.0 - sigma) * f[base + static_cast<std::size_t>(q)]
                    + sigma * target;
            }
        }
    }
}

void apply_smoke_boundaries_and_source(
    Field& g,
    const Field& ux,
    const Field& uy,
    const Mask& body
)
{
    for (int j = 1; j < NY - 1; ++j) {
        const std::size_t base = cell_id(0, j) * static_cast<std::size_t>(Q);

        for (int q = 0; q < Q; ++q) {
            g[base + static_cast<std::size_t>(q)] = 0.0;
        }
    }

    for (int j = 1; j < NY - 1; ++j) {
        for (int q = 0; q < Q; ++q) {
            g[dist_id(NX - 1, j, q)] = g[dist_id(NX - 2, j, q)];
        }
    }

    for (int i = 0; i < NX; ++i) {
        g[dist_id(i, 0, 2)] = g[dist_id(i, 0, 4)];
        g[dist_id(i, 0, 5)] = g[dist_id(i, 0, 8)];
        g[dist_id(i, 0, 6)] = g[dist_id(i, 0, 7)];

        g[dist_id(i, NY - 1, 4)] = g[dist_id(i, NY - 1, 2)];
        g[dist_id(i, NY - 1, 7)] = g[dist_id(i, NY - 1, 6)];
        g[dist_id(i, NY - 1, 8)] = g[dist_id(i, NY - 1, 5)];
    }

    int xs = SMOKE_SOURCE_X;

    if (xs < 2) {
        xs = 2;
    }

    if (xs > NX - 3) {
        xs = NX - 3;
    }

    for (int j = 1; j < NY - 1; ++j) {
        const std::size_t id = cell_id(xs, j);

        if (body[id] != 0u) {
            continue;
        }

        const double dy = static_cast<double>(j - CY);

        const double phi =
            std::exp(-0.5 * dy * dy / (SMOKE_SOURCE_SIGMA * SMOKE_SOURCE_SIGMA));

        const std::size_t base = id * static_cast<std::size_t>(Q);

        for (int q = 0; q < Q; ++q) {
            g[base + static_cast<std::size_t>(q)] =
                geq(q, phi, ux[id], uy[id]);
        }
    }

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);

            if (body[id] == 0u) {
                continue;
            }

            const std::size_t base = id * static_cast<std::size_t>(Q);

            for (int q = 0; q < Q; ++q) {
                g[base + static_cast<std::size_t>(q)] = 0.0;
            }
        }
    }
}


// ============================================================================
// Diagnostics
// ============================================================================

struct Diagnostics
{
    double rho_min;
    double rho_max;
    double speed_max;
    bool finite_ok;
    bool physics_ok;

    Diagnostics()
        : rho_min(std::numeric_limits<double>::max()),
        rho_max(-std::numeric_limits<double>::max()),
        speed_max(0.0),
        finite_ok(true),
        physics_ok(true)
    {
    }
};

Diagnostics compute_diagnostics(
    const Field& rho,
    const Field& ux,
    const Field& uy,
    const Mask& body
)
{
    Diagnostics d;

    const int margin_x = 8;
    const int margin_y = 8;

    for (int i = margin_x; i < NX - margin_x; ++i) {
        for (int j = margin_y; j < NY - margin_y; ++j) {
            const std::size_t id = cell_id(i, j);

            if (body[id] != 0u) {
                continue;
            }

            const double r = rho[id];
            const double u = ux[id];
            const double v = uy[id];

            if (!finite_number(r) || !finite_number(u) || !finite_number(v)) {
                d.finite_ok = false;
                d.physics_ok = false;
                continue;
            }

            d.rho_min = std::min(d.rho_min, r);
            d.rho_max = std::max(d.rho_max, r);

            const double s = std::sqrt(u * u + v * v);

            d.speed_max = std::max(d.speed_max, s);
        }
    }

    if (d.rho_min < 0.85 || d.rho_max > 1.15 || d.speed_max > 0.25) {
        d.physics_ok = false;
    }

    return d;
}

void throw_if_unstable(
    const int step,
    const std::string& label,
    const Diagnostics& d
)
{
    if (d.finite_ok && d.physics_ok) {
        return;
    }

    std::ostringstream msg;

    msg << "\nERROR: simulation left the trusted physical range in " << label << "\n"
        << "  step      = " << step << "\n"
        << "  rho_min   = " << d.rho_min << "\n"
        << "  rho_max   = " << d.rho_max << "\n"
        << "  |u|max    = " << d.speed_max << "\n\n"
        << "This is not a valid run. Reduce U_INF or increase domain/sponge width.\n";

    throw std::runtime_error(msg.str());
}


// ============================================================================
// Output
// ============================================================================

void clean_output_folder(const std::string& outdir, const std::string& name)
{
    MAKE_DIR(outdir.c_str());

    std::remove((outdir + "/" + name + "_forces.csv").c_str());
    std::remove((outdir + "/run_info.txt").c_str());

    for (int n = 0; n < 20000; ++n) {
        std::ostringstream ss;

        ss << outdir << "/snap_"
            << std::setw(4) << std::setfill('0') << n
            << ".csv";

        std::remove(ss.str().c_str());
    }
}

void write_run_info(const std::string& outdir, const Shape shape)
{
    std::ofstream file((outdir + "/run_info.txt").c_str());

    file << "LBM D2Q9 TRT external bluff-body solver\n";
    file << "Author: Saiyed Mudassir\n\n";

    file << "shape              = " << shape_name(shape) << "\n";
    file << "NX                 = " << NX << "\n";
    file << "NY                 = " << NY << "\n";
    file << "D                  = " << D << "\n";
    file << "R                  = " << R << "\n";
    file << "CX                 = " << CX << "\n";
    file << "CY                 = " << CY << "\n";
    file << "U_INF              = " << U_INF << "\n";
    file << "U_REF              = " << U_REF << "\n";
    file << "RE                 = " << RE << "\n";
    file << "NU                 = " << NU << "\n";
    file << "TAU_PLUS           = " << TAU_PLUS << "\n";
    file << "TAU_MINUS          = " << TAU_MINUS << "\n";
    file << "OMEGA_PLUS         = " << OMEGA_PLUS << "\n";
    file << "OMEGA_MINUS        = " << OMEGA_MINUS << "\n";
    file << "TRT_LAMBDA         = " << TRT_LAMBDA << "\n";
    file << "MA                 = " << MA << "\n";
    file << "blockage D/NY      = " << static_cast<double>(D) / static_cast<double>(NY) << "\n";
    file << "NWARM              = " << NWARM << "\n";
    file << "NRECORD            = " << NRECORD << "\n";
    file << "SNAP_EVERY         = " << SNAP_EVERY << "\n";
    file << "FORCE_EVERY        = " << FORCE_EVERY << "\n";
    file << "expected_snapshots = " << NRECORD / SNAP_EVERY << "\n";
    file << "expected_forces    = " << NRECORD / FORCE_EVERY << "\n";
    file << "RAMP_STEPS         = " << RAMP_STEPS << "\n";
    file << "SEED_STEPS         = " << SEED_STEPS << "\n";
    file << "SMOKE_TAU          = " << SMOKE_TAU << "\n";
    file << "SMOKE_SOURCE_X     = " << SMOKE_SOURCE_X << "\n";
    file << "SMOKE_SIGMA        = " << SMOKE_SOURCE_SIGMA << "\n";
    file << "SPONGE_X_WIDTH     = " << SPONGE_X_WIDTH << "\n";
    file << "SPONGE_Y_WIDTH     = " << SPONGE_Y_WIDTH << "\n";
    file << "SPONGE_STRENGTH    = " << SPONGE_STRENGTH << "\n";
}

Field compute_curl(
    const Field& ux,
    const Field& uy,
    const Mask& body
)
{
    Field curl(static_cast<std::size_t>(NX) * static_cast<std::size_t>(NY), 0.0);

    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            const std::size_t id = cell_id(i, j);

            if (body[id] != 0u) {
                continue;
            }

            const std::size_t il = cell_id(i - 1, j);
            const std::size_t ir = cell_id(i + 1, j);
            const std::size_t ib = cell_id(i, j - 1);
            const std::size_t it = cell_id(i, j + 1);

            if (body[il] != 0u || body[ir] != 0u ||
                body[ib] != 0u || body[it] != 0u) {
                continue;
            }

            const double duy_dx = 0.5 * (uy[ir] - uy[il]);
            const double dux_dy = 0.5 * (ux[it] - ux[ib]);

            curl[id] = duy_dx - dux_dy;
        }
    }

    return curl;
}

void write_snapshot(
    const std::string& outdir,
    const int snap_id,
    const Field& rho,
    const Field& ux,
    const Field& uy,
    const Field& smoke,
    const Mask& body
)
{
    const Field curl = compute_curl(ux, uy, body);

    double p_inf = 0.0;
    int n_inf = 0;

    const int iref = 4;

    for (int j = 2; j < NY - 2; ++j) {
        const std::size_t id = cell_id(iref, j);

        if (body[id] == 0u) {
            p_inf += rho[id] * CS2;
            ++n_inf;
        }
    }

    p_inf = n_inf > 0 ? p_inf / static_cast<double>(n_inf) : CS2;

    std::ostringstream ss;

    ss << outdir << "/snap_"
        << std::setw(4) << std::setfill('0') << snap_id
        << ".csv";

    std::ofstream file(ss.str().c_str());

    file << "i,j,x,y,rho,ux,uy,curl,p,Cp,solid,smoke\n";
    file << std::fixed << std::setprecision(9);

    const double q_ref = 0.5 * RHO_REF * U_REF * U_REF;

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);

            const double p = rho[id] * CS2;
            const double Cp = (p - p_inf) / q_ref;

            file
                << i << ","
                << j << ","
                << i << ","
                << j << ","
                << rho[id] << ","
                << ux[id] << ","
                << uy[id] << ","
                << curl[id] << ","
                << p << ","
                << Cp << ","
                << static_cast<int>(body[id]) << ","
                << smoke[id] << "\n";
        }
    }
}


// ============================================================================
// Initialisation
// ============================================================================

void initialise_fluid(Field& f, const Mask& body)
{
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const std::size_t id = cell_id(i, j);
            const std::size_t base = id * static_cast<std::size_t>(Q);

            double r = RHO_REF;
            double u = 0.0;
            double v = 0.0;

            if (body[id] != 0u) {
                u = 0.0;
                v = 0.0;
            }

            for (int q = 0; q < Q; ++q) {
                f[base + static_cast<std::size_t>(q)] = feq(q, r, u, v);
            }
        }
    }
}

void initialise_smoke(Field& g)
{
    std::fill(g.begin(), g.end(), 0.0);
}


// ============================================================================
// Solver
// ============================================================================

void run_case(const Shape shape)
{
    const std::string name = shape_name(shape);
    const std::string label = shape_label(shape);
    const std::string outdir = name + "_snapshots";

    clean_output_folder(outdir, name);
    write_run_info(outdir, shape);

    const Mask body = build_body_mask(shape);

    const std::size_t num_cells =
        static_cast<std::size_t>(NX) * static_cast<std::size_t>(NY);

    const std::size_t num_dist =
        num_cells * static_cast<std::size_t>(Q);

    Field f(num_dist, 0.0);
    Field f_post(num_dist, 0.0);

    Field g(num_dist, 0.0);
    Field g_post(num_dist, 0.0);

    Field rho(num_cells, RHO_REF);
    Field ux(num_cells, 0.0);
    Field uy(num_cells, 0.0);
    Field smoke(num_cells, 0.0);

    initialise_fluid(f, body);
    initialise_smoke(g);

    compute_macros(f, body, rho, ux, uy);
    apply_smoke_boundaries_and_source(g, ux, uy, body);

    const std::string force_path = outdir + "/" + name + "_forces.csv";

    std::ofstream force_file(force_path.c_str());

    force_file << "step,t_star,Cd,Cl,Fx_raw,Fy_raw\n";
    force_file << std::fixed << std::setprecision(9);

    std::cout << "\n=== " << label << " ===\n";
    std::cout << "  folder       = " << outdir << "\n";
    std::cout << "  expected CSV = " << NRECORD / SNAP_EVERY
        << " snapshots, " << NRECORD / FORCE_EVERY
        << " force rows\n";

    const int total_steps = NWARM + NRECORD;
    int snap_count = 0;

    for (int step = 0; step < total_steps; ++step) {
        compute_macros(f, body, rho, ux, uy);
        apply_wake_seed(ux, uy, body, step);

        collide_fluid_trt(f, f_post, rho, ux, uy, body);

        const Forces force = stream_fluid(f_post, f, body);

        apply_inlet_zou_he(f, step);
        apply_outlet_pressure_zou_he(f);
        apply_top_bottom_slip(f);
        reset_body_distributions(f, body);
        apply_farfield_sponge(f, body, step);

        compute_macros(f, body, rho, ux, uy);

        compute_smoke(g, body, smoke);
        collide_smoke_bgk(g, g_post, smoke, ux, uy, body);
        stream_smoke(g_post, g, body);
        apply_smoke_boundaries_and_source(g, ux, uy, body);

        const bool in_record = step >= NWARM;
        const int record_step = step - NWARM;

        if (step % 1000 == 0) {
            const Diagnostics d = compute_diagnostics(rho, ux, uy, body);

            std::cout
                << "  step " << std::setw(7) << step
                << "  ramp=" << std::setprecision(4) << ramp_factor(step)
                << "  rho=[" << std::setprecision(6) << d.rho_min
                << ", " << d.rho_max << "]"
                << "  |u|max=" << d.speed_max
                << "\n";

            throw_if_unstable(step, label, d);
        }

        if (in_record && record_step % FORCE_EVERY == 0) {
            const double q_ref_force =
                0.5 * RHO_REF * U_REF * U_REF * static_cast<double>(D);

            const double Cd = force.Fx / q_ref_force;
            const double Cl = force.Fy / q_ref_force;

            const double t_star =
                static_cast<double>(record_step)
                * U_REF
                / static_cast<double>(D);

            force_file
                << step << ","
                << t_star << ","
                << Cd << ","
                << Cl << ","
                << force.Fx << ","
                << force.Fy << "\n";
        }

        if (in_record && record_step % SNAP_EVERY == 0) {
            compute_macros(f, body, rho, ux, uy);
            compute_smoke(g, body, smoke);

            const Diagnostics d = compute_diagnostics(rho, ux, uy, body);
            throw_if_unstable(step, label, d);

            write_snapshot(outdir, snap_count, rho, ux, uy, smoke, body);

            if (snap_count % 20 == 0) {
                const double t_star =
                    static_cast<double>(record_step)
                    * U_REF
                    / static_cast<double>(D);

                std::cout
                    << "  saved snap " << std::setw(4) << snap_count
                    << "  step " << step
                    << "  t*=" << t_star
                    << "\n";
            }

            ++snap_count;
        }
    }

    force_file.close();

    std::cout << "  Done. " << snap_count
        << " snapshots written to " << outdir << "\n";
}


// ============================================================================
// Main
// ============================================================================

int main()
{
    try {
        std::cout << "============================================================\n";
        std::cout << " LBM D2Q9 TRT External Bluff-Body Aerodynamics\n";
        std::cout << " Circle vs Square vs Droplet\n";
        std::cout << " Author: Saiyed Mudassir\n";
        std::cout << "============================================================\n\n";

        std::cout << "RUN VERIFICATION\n";
        std::cout << "  NX, NY           = " << NX << ", " << NY << "\n";
        std::cout << "  D                = " << D << "\n";
        std::cout << "  CX, CY           = " << CX << ", " << CY << "\n";
        std::cout << "  U_INF            = " << U_INF << "\n";
        std::cout << "  Re               = " << RE << "\n";
        std::cout << "  nu               = " << NU << "\n";
        std::cout << "  tau_plus         = " << TAU_PLUS << "\n";
        std::cout << "  tau_minus        = " << TAU_MINUS << "\n";
        std::cout << "  Mach             = " << MA << "\n";
        std::cout << "  blockage D/NY    = "
            << static_cast<double>(D) / static_cast<double>(NY) << "\n";
        std::cout << "  NWARM            = " << NWARM << "\n";
        std::cout << "  NRECORD          = " << NRECORD << "\n";
        std::cout << "  SNAP_EVERY       = " << SNAP_EVERY << "\n";
        std::cout << "  FORCE_EVERY      = " << FORCE_EVERY << "\n";
        std::cout << "  expected snaps   = " << NRECORD / SNAP_EVERY << "\n";
        std::cout << "  expected forces  = " << NRECORD / FORCE_EVERY << "\n";
        std::cout << "  RAMP_STEPS       = " << RAMP_STEPS << "\n";
        std::cout << "  SEED_STEPS       = " << SEED_STEPS << "\n";
        std::cout << "  SMOKE_TAU        = " << SMOKE_TAU << "\n";
        std::cout << "  SMOKE_SOURCE_X   = " << SMOKE_SOURCE_X << "\n";
        std::cout << "  SMOKE_SIGMA      = " << SMOKE_SOURCE_SIGMA << "\n";
        std::cout << "  SPONGE_X_WIDTH   = " << SPONGE_X_WIDTH << "\n";
        std::cout << "  SPONGE_Y_WIDTH   = " << SPONGE_Y_WIDTH << "\n";
        std::cout << "  SPONGE_STRENGTH  = " << SPONGE_STRENGTH << "\n\n";

        if (MA >= 0.15) {
            throw std::runtime_error("Mach number too high. Reduce U_INF.");
        }

        if (TAU_PLUS <= 0.52) {
            throw std::runtime_error("TAU_PLUS too close to 0.5. Increase U_INF or D, or reduce Re.");
        }

        if (SMOKE_TAU <= 0.5) {
            throw std::runtime_error("SMOKE_TAU must be greater than 0.5.");
        }

        if (RUN_CIRCLE) {
            run_case(Shape::Circle);
        }

        if (RUN_SQUARE) {
            run_case(Shape::Square);
        }

        if (RUN_DROPLET) {
            run_case(Shape::Droplet);
        }

        std::cout << "\nAll done.\n";
        std::cout << "Run Post_processing.py from this same working directory.\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nSolver stopped.\n";
        std::cerr << e.what() << "\n";
        return 1;
    }
}