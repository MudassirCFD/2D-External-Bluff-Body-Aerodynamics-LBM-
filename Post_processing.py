"""
Post processing for the 2D LBM wake shape study.

The script reads the solver output folders:

    circle_snapshots
    square_snapshots
    droplet_snapshots

It writes selected portfolio figures into:

    figures
    results

"""

from pathlib import Path
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation


CASES = [
    ("circle", "Circle"),
    ("square", "Square"),
    ("droplet", "Droplet"),
]

FIG_DIR = Path("figures")
RES_DIR = Path("results")

MEAN_LAST_N_SNAPSHOTS = 30
AVERAGE_LAST_FRACTION = 0.50
WAKE_PROBE_D = 5.0
ANIMATION_MAX_FRAMES = 120
ANIMATION_FPS = 10


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------


def make_folders():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)


def find_case_folder(case_name):
    """Find the solver output folder for a case."""
    possible = [
        Path(case_name + "_snapshots"),
        Path(case_name + "_snapshots_VALID"),
        Path(case_name + "_output"),
    ]

    for folder in possible:
        if folder.exists():
            return folder

    for folder in Path.cwd().rglob(case_name + "_snapshots"):
        return folder

    for folder in Path.cwd().rglob(case_name + "_snapshots_VALID"):
        return folder

    return None


def read_run_info(folder):
    info = {}
    path = folder / "run_info.txt"

    if not path.exists():
        return info

    for line in path.read_text(errors="ignore").splitlines():
        if "=" not in line:
            continue

        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()

        try:
            if "." in value or "e" in value.lower():
                info[key] = float(value)
            else:
                info[key] = int(value)
        except ValueError:
            info[key] = value

    return info


def load_force(folder, case_name):
    path = folder / (case_name + "_forces.csv")

    if not path.exists():
        raise FileNotFoundError("Missing force file: " + str(path))

    df = pd.read_csv(path)

    if "t_star" not in df.columns:
        df["t_star"] = df["step"]

    if "Cd" not in df.columns:
        raise ValueError("Force file has no Cd column: " + str(path))

    if "Cl" not in df.columns:
        raise ValueError("Force file has no Cl column: " + str(path))

    return df


def read_snapshot(path):
    return pd.read_csv(path)


def grid_from_snapshot(df):
    """
    Convert one snapshot table into 2D arrays.

    The C++ solver writes data with i as the outer loop and j as the inner loop.
    That means a direct reshape must be done as nx by ny, then transposed to y by x.
    """
    nx = int(df["i"].max()) + 1
    ny = int(df["j"].max()) + 1

    x = np.arange(nx)
    y = np.arange(ny)

    fields = {}
    for col in df.columns:
        if col in ["i", "j", "x", "y"]:
            continue

        values = df[col].to_numpy()
        fields[col] = values.reshape((nx, ny)).T

    return x, y, fields


def body_mask(fields):
    if "solid" in fields:
        return fields["solid"] > 0.5
    if "body" in fields:
        return fields["body"] > 0.5
    raise ValueError("No solid or body column was found in the snapshot")


def masked(field, body):
    out = np.array(field, dtype=float)
    out[body] = np.nan
    return out


def safe_limits(arrays, lower=1.0, upper=99.0):
    vals = []

    for arr in arrays:
        a = np.asarray(arr, dtype=float)
        a = a[np.isfinite(a)]
        if a.size:
            vals.append(a)

    if not vals:
        return 0.0, 1.0

    vals = np.concatenate(vals)
    lo = float(np.percentile(vals, lower))
    hi = float(np.percentile(vals, upper))

    if abs(hi - lo) < 1.0e-12:
        hi = lo + 1.0

    return lo, hi


def symmetric_limits(arrays, percentile=98.5):
    vals = []

    for arr in arrays:
        a = np.asarray(arr, dtype=float)
        a = a[np.isfinite(a)]
        if a.size:
            vals.append(a)

    if not vals:
        return -1.0, 1.0

    vals = np.concatenate(vals)
    lim = float(np.percentile(np.abs(vals), percentile))

    if lim < 1.0e-12:
        lim = 1.0

    return -lim, lim


def save_current(path):
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def moving_average(values, window=21):
    values = np.asarray(values, dtype=float)

    if len(values) < window:
        return values

    kernel = np.ones(window) / float(window)
    return np.convolve(values, kernel, mode="same")


def add_body_outline(ax, x, y, body):
    try:
        ax.contour(x, y, body.astype(float), levels=[0.5], linewidths=1.2)
    except Exception:
        pass


# -----------------------------------------------------------------------------
# Load all data
# -----------------------------------------------------------------------------


def load_all_cases():
    data = {}

    for case_name, label in CASES:
        folder = find_case_folder(case_name)

        if folder is None:
            print("Missing folder for case:", case_name)
            continue

        force = load_force(folder, case_name)
        snaps = sorted(folder.glob("snap_*.csv"))
        info = read_run_info(folder)

        if not snaps:
            print("No snapshots found for case:", case_name)

        data[case_name] = {
            "label": label,
            "folder": folder,
            "force": force,
            "snaps": snaps,
            "info": info,
        }

        print(label)
        print("  folder:", folder)
        print("  force rows:", len(force))
        print("  snapshots:", len(snaps))

    if not data:
        raise RuntimeError("No solver output folders were found")

    return data


# -----------------------------------------------------------------------------
# Force statistics
# -----------------------------------------------------------------------------


def force_statistics(force):
    start = int((1.0 - AVERAGE_LAST_FRACTION) * len(force))
    d = force.iloc[start:].copy()

    cl_mean = float(d["Cl"].mean())
    cl_fluct = d["Cl"] - cl_mean

    return {
        "Mean Cd": float(d["Cd"].mean()),
        "Std Cd": float(d["Cd"].std()),
        "Mean Cl": cl_mean,
        "RMS Cl": float(np.sqrt(np.mean(cl_fluct ** 2))),
        "Std Cl": float(d["Cl"].std()),
        "Samples used": int(len(d)),
    }


def estimate_strouhal(force, info):
    """Estimate Strouhal number from the strongest lift frequency."""
    if len(force) < 20:
        return float("nan")

    t = force["t_star"].to_numpy(dtype=float)
    cl = force["Cl"].to_numpy(dtype=float)

    start = int((1.0 - AVERAGE_LAST_FRACTION) * len(force))
    t = t[start:]
    cl = cl[start:]
    cl = cl - np.mean(cl)

    if len(t) < 20 or np.allclose(cl, 0.0):
        return float("nan")

    dt = float(np.mean(np.diff(t)))

    if dt <= 0.0:
        return float("nan")

    freqs = np.fft.rfftfreq(len(cl), d=dt)
    amp = np.abs(np.fft.rfft(cl))

    if len(freqs) <= 1:
        return float("nan")

    # skip zero frequency
    idx = int(np.argmax(amp[1:]) + 1)
    return float(freqs[idx])


def write_results_summary(data):
    rows = []

    for case_name, case in data.items():
        row = {"Body": case["label"]}
        row.update(force_statistics(case["force"]))
        row["Strouhal estimate"] = estimate_strouhal(case["force"], case["info"])
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(RES_DIR / "results_summary.csv", index=False)

    with open(RES_DIR / "results_summary.txt", "w", encoding="utf-8") as f:
        f.write("2D LBM Wake Shape Study\n")
        f.write("=======================\n\n")
        f.write("Statistics are taken from the final recorded part of the force history.\n")
        f.write("All cases use the same reference size and same Reynolds number.\n\n")
        f.write(summary.to_string(index=False))
        f.write("\n")

    return summary


# -----------------------------------------------------------------------------
# Main figures
# -----------------------------------------------------------------------------


def plot_geometry(data):
    fig, axes = plt.subplots(1, len(data), figsize=(5.8 * len(data), 4.2))

    if len(data) == 1:
        axes = [axes]

    for ax, (case_name, case) in zip(axes, data.items()):
        if not case["snaps"]:
            continue

        df = read_snapshot(case["snaps"][0])
        x, y, f = grid_from_snapshot(df)
        body = body_mask(f)

        ax.imshow(body, origin="lower", extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
        ax.set_title(case["label"])
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    save_current(FIG_DIR / "00_geometry_check.png")


def final_fields(data):
    final = {}

    for case_name, case in data.items():
        if not case["snaps"]:
            continue

        df = read_snapshot(case["snaps"][-1])
        x, y, f = grid_from_snapshot(df)
        body = body_mask(f)

        if "curl" in f:
            vort = f["curl"]
        elif "vorticity" in f:
            vort = f["vorticity"]
        else:
            vort = np.zeros_like(f["ux"])

        final[case_name] = {
            "label": case["label"],
            "x": x,
            "y": y,
            "body": body,
            "ux": f["ux"],
            "uy": f["uy"],
            "speed": np.sqrt(f["ux"] ** 2 + f["uy"] ** 2),
            "curl": vort,
            "Cp": f["Cp"],
            "smoke": f.get("smoke", np.zeros_like(f["ux"])),
        }

    return final


def plot_field_group(items, key, title, filename, symmetric=False, lower=1.0, upper=99.0):
    arrays = [masked(item[key], item["body"]) for item in items.values()]

    if symmetric:
        vmin, vmax = symmetric_limits(arrays)
    else:
        vmin, vmax = safe_limits(arrays, lower, upper)

    fig, axes = plt.subplots(1, len(items), figsize=(5.8 * len(items), 4.2))

    if len(items) == 1:
        axes = [axes]

    for ax, item in zip(axes, items.values()):
        x = item["x"]
        y = item["y"]
        body = item["body"]
        arr = masked(item[key], body)

        im = ax.imshow(
            arr,
            origin="lower",
            extent=[x.min(), x.max(), y.min(), y.max()],
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )

        add_body_outline(ax, x, y, body)
        ax.set_title(item["label"] + " " + title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax)

    save_current(FIG_DIR / filename)


def plot_final_fields(data):
    items = final_fields(data)

    if not items:
        return

    plot_field_group(items, "curl", "vorticity", "01_vorticity_hero.png", symmetric=True)
    plot_field_group(items, "Cp", "pressure coefficient", "02_pressure_coefficient.png", lower=2.0, upper=98.0)
    plot_field_group(items, "speed", "speed field", "03_speed_field.png", lower=1.0, upper=99.0)


def plot_force_coefficients(data):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    for case_name, case in data.items():
        f = case["force"]
        t = f["t_star"]

        axes[0].plot(t, f["Cd"], alpha=0.35, label=case["label"] + " raw")
        axes[0].plot(t, moving_average(f["Cd"], 21), linewidth=2.0, label=case["label"] + " smooth")

        axes[1].plot(t, f["Cl"], alpha=0.35, label=case["label"] + " raw")
        axes[1].plot(t, moving_average(f["Cl"], 21), linewidth=2.0, label=case["label"] + " smooth")

    axes[0].set_ylabel("Cd")
    axes[0].set_title("Drag coefficient history")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].set_xlabel("Non dimensional time, tU/D")
    axes[1].set_ylabel("Cl")
    axes[1].set_title("Lift coefficient history")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    save_current(FIG_DIR / "04_force_coefficients.png")


def plot_engineering_summary(summary):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    bodies = summary["Body"]

    axes[0].bar(bodies, summary["Mean Cd"])
    axes[0].set_title("Mean drag coefficient")
    axes[0].set_ylabel("Mean Cd")
    axes[0].grid(True, axis="y", alpha=0.3)

    axes[1].bar(bodies, summary["RMS Cl"])
    axes[1].set_title("RMS lift coefficient")
    axes[1].set_ylabel("RMS Cl")
    axes[1].grid(True, axis="y", alpha=0.3)

    axes[2].bar(bodies, summary["Strouhal estimate"])
    axes[2].set_title("Strouhal estimate")
    axes[2].set_ylabel("St")
    axes[2].grid(True, axis="y", alpha=0.3)

    save_current(FIG_DIR / "05_engineering_summary.png")


def mean_fields(data):
    out = {}

    for case_name, case in data.items():
        snaps = case["snaps"][-MEAN_LAST_N_SNAPSHOTS:]

        if not snaps:
            continue

        speed_sum = None
        curl_sum = None
        ux_sum = None
        count = 0
        x = None
        y = None
        body = None

        for snap in snaps:
            df = read_snapshot(snap)
            x, y, f = grid_from_snapshot(df)
            body = body_mask(f)
            speed = np.sqrt(f["ux"] ** 2 + f["uy"] ** 2)
            curl = f["curl"] if "curl" in f else np.zeros_like(speed)

            if speed_sum is None:
                speed_sum = np.zeros_like(speed)
                curl_sum = np.zeros_like(curl)
                ux_sum = np.zeros_like(f["ux"])

            speed_sum += np.nan_to_num(masked(speed, body))
            curl_sum += np.nan_to_num(masked(curl, body))
            ux_sum += np.nan_to_num(masked(f["ux"], body))
            count += 1

        out[case_name] = {
            "label": case["label"],
            "x": x,
            "y": y,
            "body": body,
            "speed": speed_sum / float(count),
            "curl": curl_sum / float(count),
            "ux": ux_sum / float(count),
            "info": case["info"],
        }

    return out


def plot_wake_profiles(data):
    means = mean_fields(data)

    if not means:
        return

    plt.figure(figsize=(8, 6))

    for case_name, item in means.items():
        info = item["info"]
        x = item["x"]
        y = item["y"]
        ux = masked(item["ux"], item["body"])

        cx = int(info.get("CX", 180))
        d = float(info.get("D", 40))
        u_inf = float(info.get("U_INF", 0.06))

        probe_x = int(min(len(x) - 2, max(1, cx + WAKE_PROBE_D * d)))
        deficit = u_inf - ux[:, probe_x]

        plt.plot(deficit, y, label=item["label"] + " at x " + str(probe_x))

    plt.xlabel("Wake deficit, U_inf minus u")
    plt.ylabel("y")
    plt.title("Mean wake deficit profile")
    plt.grid(True, alpha=0.3)
    plt.legend()
    save_current(FIG_DIR / "06_wake_profiles.png")


def plot_mean_wake(data):
    means = mean_fields(data)

    if not means:
        return

    arrays = [masked(item["speed"], item["body"]) for item in means.values()]
    vmin, vmax = safe_limits(arrays, 1.0, 99.0)

    fig, axes = plt.subplots(1, len(means), figsize=(5.8 * len(means), 4.2))

    if len(means) == 1:
        axes = [axes]

    for ax, item in zip(axes, means.values()):
        x = item["x"]
        y = item["y"]
        body = item["body"]
        arr = masked(item["speed"], body)

        im = ax.imshow(
            arr,
            origin="lower",
            extent=[x.min(), x.max(), y.min(), y.max()],
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )

        add_body_outline(ax, x, y, body)
        ax.set_title(item["label"] + " mean wake speed")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax)

    save_current(FIG_DIR / "07_mean_wake.png")


# -----------------------------------------------------------------------------
# Animation
# -----------------------------------------------------------------------------


def available_writer():
    try:
        if animation.writers.is_available("ffmpeg"):
            return "ffmpeg"
    except Exception:
        pass

    try:
        if animation.writers.is_available("pillow"):
            return "pillow"
    except Exception:
        pass

    return None


def make_wake_animation(data):
    writer = available_writer()

    if writer is None:
        print("No animation writer was found. Skipping wake video.")
        return

    usable = {name: case for name, case in data.items() if case["snaps"]}

    if not usable:
        return

    # Use the same number of frames from each case.
    min_snaps = min(len(case["snaps"]) for case in usable.values())
    frame_count = min(min_snaps, ANIMATION_MAX_FRAMES)
    indices = np.linspace(0, min_snaps - 1, frame_count).astype(int)

    first_items = []
    for case in usable.values():
        df = read_snapshot(case["snaps"][indices[0]])
        x, y, f = grid_from_snapshot(df)
        body = body_mask(f)
        curl = f["curl"] if "curl" in f else np.zeros_like(f["ux"])
        first_items.append((case["label"], x, y, masked(curl, body), body))

    lim_arrays = [item[3] for item in first_items]
    vmin, vmax = symmetric_limits(lim_arrays)

    fig, axes = plt.subplots(1, len(usable), figsize=(5.8 * len(usable), 4.2))

    if len(usable) == 1:
        axes = [axes]

    images = []

    for ax, item in zip(axes, first_items):
        label, x, y, curl, body = item
        im = ax.imshow(
            curl,
            origin="lower",
            extent=[x.min(), x.max(), y.min(), y.max()],
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )
        add_body_outline(ax, x, y, body)
        ax.set_title(label + " vorticity")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        images.append(im)

    plt.tight_layout()

    def update(frame_number):
        for im, (case_name, case) in zip(images, usable.items()):
            snap = case["snaps"][indices[frame_number]]
            df = read_snapshot(snap)
            x, y, f = grid_from_snapshot(df)
            body = body_mask(f)
            curl = f["curl"] if "curl" in f else np.zeros_like(f["ux"])
            im.set_data(masked(curl, body))

        return images

    ani = animation.FuncAnimation(fig, update, frames=frame_count, interval=1000.0 / ANIMATION_FPS)

    if writer == "ffmpeg":
        ani.save(FIG_DIR / "08_wake_evolution.mp4", writer="ffmpeg", fps=ANIMATION_FPS, dpi=140)
    else:
        ani.save(FIG_DIR / "08_wake_evolution.gif", writer="pillow", fps=ANIMATION_FPS, dpi=120)

    plt.close(fig)


# -----------------------------------------------------------------------------
# Main script
# -----------------------------------------------------------------------------


def main():
    make_folders()
    data = load_all_cases()

    summary = write_results_summary(data)

    plot_geometry(data)
    plot_final_fields(data)
    plot_force_coefficients(data)
    plot_engineering_summary(summary)
    plot_wake_profiles(data)
    plot_mean_wake(data)
    make_wake_animation(data)

    print()
    print("Post processing complete")
    print("Figures folder:", FIG_DIR.resolve())
    print("Results folder:", RES_DIR.resolve())
    print()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
