# 2D LBM Wake Shape Study

Circle, square and droplet body comparison at Re 200 using a C++ D2Q9 TRT Lattice Boltzmann solver.

This repository shows a small but complete numerical simulation study. I used it to compare how simple body shape changes separation, vortex shedding, wake width, drag and unsteady lift.

This is not presented as a fully validated solver or a final design claim. It is a controlled numerical study where the same solver, domain, Reynolds number, inlet condition and post processing method are used for all three bodies. Only the body shape is changed.

## Honesty about the solver

I am not claiming that every basic idea in this solver was written from a blank page.

The starting point came from standard educational D2Q9 Lattice Boltzmann examples, cylinder wake examples and the main LBM literature. In particular, the work of Jonas Latt and other LBM researchers helped me understand the basic D2Q9 layout, streaming, collision and cylinder wake style setup [1].

From that starting point, I rebuilt and extended the project for my own study. The final repository includes my own case setup, geometry control for circle, square and droplet bodies, TRT collision structure, force extraction, passive smoke visualisation, wake seeding, sponge damping, Python post processing, force comparison, time averaged wake plots and engineering interpretation.

The aim is to show that I understand the simulation chain, not to pretend that the basic LBM method or standard cylinder example came from me.

## What this study compares

The three bodies are:

* Circle
* Square
* Droplet

The same flow condition is used in every case.

The comparison focuses on:

* wake structure
* vorticity shedding
* pressure coefficient field
* speed field
* mean drag coefficient
* RMS lift coefficient
* Strouhal estimate
* mean wake deficit
* wake profile downstream of the body

The main question is simple:

How does the shape change separation, pressure recovery, wake size and unsteady loading?

## Why this is an unsteady study

This is not a steady RANS type calculation.

A steady RANS calculation would mainly give a time averaged solution. That can be useful, but it can hide the actual vortex shedding behind a bluff body.

This project is closer in spirit to an unsteady wake calculation. It stores the force history and flow snapshots over time, so the vortex shedding and lift fluctuation can be compared directly.

It is not URANS either, because it does not solve the Reynolds averaged Navier Stokes equations with a turbulence model. It uses a Lattice Boltzmann Method. The useful comparison is only this:

A steady method would hide much of the wake motion.

This LBM run keeps the unsteady behaviour visible.

That matters because at Re 200, bluff bodies naturally form an unsteady wake with vortex shedding [4].

## Solver method

The solver uses the Lattice Boltzmann Method.

Instead of solving the Navier Stokes equations directly in pressure and velocity form, LBM works with distribution functions on a lattice. These distributions stream from one node to another and then relax towards equilibrium during collision [2].

The macroscopic density and velocity are recovered from the distribution functions:

rho = sum of f_i

rho u = sum of f_i e_i

where f_i is the particle distribution in direction i and e_i is the lattice velocity direction.

For this project I used the D2Q9 model.

D2Q9 means:

* D2 means two dimensional
* Q9 means nine velocity directions
* one rest direction
* four axial directions
* four diagonal directions

The solver uses a TRT collision model. TRT means Two Relaxation Time. It separates the symmetric and anti symmetric parts of the distribution function. This gives better numerical control than a simple BGK model, especially near wall boundaries and solid bodies [3].

## Flow setup

Main setup:

| Quantity | Value |
|---|---:|
| Method | D2Q9 TRT LBM |
| Reynolds number | 200 |
| Domain | 1000 by 320 lattice cells |
| Body reference size | 40 lattice cells |
| Inlet speed | 0.06 lattice units |
| Kinematic viscosity | 0.012 |
| Warm up steps | 30000 |
| Recorded steps | 60000 |
| Snapshot interval | 200 steps |
| Force interval | 10 steps |
| Snapshots | 300 per case |
| Force samples | 6000 per case |

The solver includes:

* smooth inlet ramp
* pressure outlet
* slip treatment at the top and bottom boundaries
* solid body bounce back
* far field sponge damping
* outlet sponge damping
* momentum exchange force calculation
* passive smoke field for wake visualisation

## Reynolds number

The Reynolds number is:

Re = U D / nu

where:

U is the inlet speed  
D is the body reference size  
nu is the kinematic viscosity  

For this case:

Re = 0.06 x 40 / 0.012

Re = 200

This Reynolds number is high enough to create a clear unsteady bluff body wake, but still small enough to keep the case simple and easy to inspect.

At Re 200, the circular cylinder wake is expected to show regular vortex shedding, so it is a useful reference case for future validation [4] [5].

## LBM viscosity

In lattice units, the viscosity is linked to the relaxation time.

For a basic LBM relation:

nu = c_s squared multiplied by tau minus 0.5

For D2Q9:

c_s squared = 1 / 3

So:

nu = one third multiplied by tau minus 0.5

This is one reason why the relaxation time must be chosen carefully. It controls the viscosity, and the viscosity controls the Reynolds number [2].

In this solver, the TRT model uses two relaxation rates rather than one. The positive relaxation rate controls the viscous part, while the second relaxation rate helps control the anti symmetric part of the distribution [3].

## Drag and lift coefficients

The solver uses momentum exchange at the solid boundary to estimate the force on the body.

The drag coefficient is:

Cd = Fx / q Aref

The lift coefficient is:

Cl = Fy / q Aref

where:

q = 0.5 rho U squared

For this two dimensional study, the reference area is treated as the body reference size multiplied by unit depth.

The important point is that the same reference value is used for all three bodies. This makes the comparison fair because the normalisation is not changed from one shape to another.

## Pressure coefficient

The pressure coefficient is used to show the pressure behaviour around each body and inside the wake.

Cp = p minus p_ref divided by q

where:

p is the local pressure  
p_ref is the reference pressure  
q is the dynamic pressure  

A strong low pressure region behind the body usually means poor pressure recovery. Poor pressure recovery is one of the main reasons bluff bodies produce high pressure drag.

## Vorticity

Vorticity shows local rotation in the flow.

For two dimensional flow:

omega = d v / d x minus d u / d y

This is one of the most important fields in this project because it shows the shear layers and the vortices released into the wake.

The vorticity plots are more useful than a single drag number because they show why the force changes.

## Strouhal number

The Strouhal number links vortex shedding frequency to the body size and inlet speed:

St = f D / U

where:

f is the shedding frequency  
D is the body reference size  
U is the inlet speed  

For bluff body flow, the Strouhal number is a useful check because it tells us whether the wake shedding rate is in a sensible range. Published circular cylinder data can be used later to check this more carefully [4] [5].

## Why the shapes behave differently

The square has sharp corners. The flow cannot turn smoothly around the corners, so separation starts early. This creates a wide wake and poor pressure recovery.

The circle is smoother than the square. The flow has a cleaner path around the body, but it is still a bluff body. It still produces a strong vortex street and unsteady lift.

The droplet has a rounded front and tapered rear. The rear taper gives the flow a better recovery path behind the body. In this setup, that reduces wake width, mean drag and lift fluctuation.

## Results

The final recorded part of the simulation gave:

| Body | Mean Cd | RMS Cl | Strouhal estimate |
|---|---:|---:|---:|
| Circle | 1.506 | 0.542 | 0.210 |
| Square | 1.698 | 0.448 | 0.174 |
| Droplet | 1.181 | 0.229 | 0.203 |

The droplet reduced mean drag by about 30.5 per cent compared with the square.

The droplet reduced mean drag by about 21.6 per cent compared with the circle.

The droplet also gave the lowest RMS lift coefficient. That means the unsteady side force from vortex shedding was weaker in this run.

## Engineering judgement

The square gives the highest drag. This makes sense because the sharp corners force early separation. The wake stays wide and pressure recovery is poor.

The circle is cleaner than the square, but it still behaves like a bluff body. The vortex street is clear and the lift fluctuation remains high.

The droplet gives the best result in this setup. The tapered rear shape helps the wake recover. The wake is narrower, mean drag is lower and RMS lift is also lower.

The main lesson is simple:

Do not only run a simulation and collect colourful plots. Read the wake, check the force history, compare the unsteady loading and then decide what actually improved.

## What I learned

This study helped me practise the full numerical workflow:

* setting up a simple solver
* controlling geometry
* keeping flow conditions consistent
* extracting forces
* checking wake behaviour
* using time averaged fields
* comparing drag with lift fluctuation
* judging whether the result makes physical sense

The most important point for me was that drag alone is not enough. A shape can only be judged properly when the wake, pressure field and unsteady loading are also checked.

## Post processing

The post processing is written in Python.

It creates the main figures and summary files:

| Output | Purpose |
|---|---|
| geometry check | confirms that the three shapes are set up correctly |
| vorticity field | shows shear layers and vortex shedding |
| pressure coefficient field | shows pressure recovery and low pressure wake regions |
| speed field | shows acceleration, wake deficit and separated flow |
| force coefficient history | shows Cd and Cl over time |
| engineering summary | gives the main comparison in one view |
| wake profiles | compares downstream velocity deficit |
| mean wake field | shows the time averaged wake |
| wake evolution video | shows the wake moving over time |
| results summary csv | stores the numerical summary |
| results summary text | stores a readable summary |

The raw snapshot folders are not included because they are large. The solver can regenerate them.

## Repository structure

2D External Bluff Body Aerodynamics LBM

README.md  
main.cpp  
Post_processing.py  
requirements.txt  
.gitignore  
LICENSE  

figures  
geometry check image  
vorticity image  
pressure coefficient image  
speed field image  
force coefficient image  
engineering summary image  
wake profile image  
mean wake image  
wake evolution video  

results  
results summary csv  
results summary text  

docs  
numerical method notes  
validation notes  
run notes  

## How to build

This project was developed in Visual Studio on Windows.

Use Release x64 for the final run because Debug is much slower.

From Developer PowerShell:

    msbuild .\2D_LBM_Unsteady.sln /t:Rebuild /p:Configuration=Release /p:Platform=x64

Then run:

    .\x64\Release\2D_LBM_Unsteady.exe

## How to run post processing

From the project root:

    python Post_processing.py

Python packages needed:

    numpy
    pandas
    matplotlib

Install them with:

    pip install numpy pandas matplotlib

## What is not uploaded

The raw snapshot folders are not uploaded because they are large.

These folders can be regenerated by running the solver:

circle snapshots  
square snapshots  
droplet snapshots  

Only the final figures, video and summary files are kept in the repository.

## Limits of this study

This is a numerical simulation study. It is not yet a final validation grade solver.

The comparison is still useful because all three cases are run with the same solver, same domain, same Reynolds number and same post processing method.

The next checks should be:

* grid independence
* domain sensitivity
* outlet sensitivity
* blockage sensitivity
* force extraction check
* comparison with published circular cylinder wake data at Re 200

## What I would improve next

The next step would be to run a proper grid study and compare the circular cylinder case with published benchmark data.

I would also check the outlet distance, top and bottom blockage, sponge strength and force extraction method. These checks are important before treating the values as validation quality results.

A stronger version of this project would include:

* coarse, medium and fine grids
* longer domain sensitivity
* wake frequency extraction using FFT
* circular cylinder benchmark comparison
* uncertainty estimate for Cd, Cl and Strouhal number
* cleaner documentation of the force extraction method

## Credits

The basic LBM ideas used in this project come from standard D2Q9 Lattice Boltzmann learning examples and the published LBM literature.

I used these resources to understand the basic method:

* Jonas Latt and the wider LBM research community for clear D2Q9 and cylinder wake style examples [1]
* Qian, d Humieres and Lallemand for the early lattice BGK formulation [2]
* Lallemand and Luo for the theory and stability behaviour of LBM schemes [3]
* Zou and He for pressure and velocity boundary condition ideas [6]
* Williamson and Roshko for bluff body wake and vortex shedding behaviour [4] [5]

My contribution in this repository was to turn the basic idea into a complete shape comparison study, add the circle, square and droplet setup, run the cases, post process the results and interpret the wake behaviour from an aerodynamic point of view.

Jonas Latt papers and profile search:  
https://arxiv.org/search/physics?searchtype=author&query=Latt,+J

## References

[1] Jonas Latt and Bastien Chopard. Lattice Boltzmann Method with regularized non equilibrium distribution functions. 2005.

[2] Qian Y H, d Humieres D and Lallemand P. Lattice BGK Models for Navier Stokes Equation. Europhysics Letters, 1992.

[3] Lallemand P and Luo L S. Theory of the lattice Boltzmann method: dispersion, dissipation, isotropy, Galilean invariance and stability. Physical Review E, 2000.

[4] Williamson C H K. Vortex dynamics in the cylinder wake. Annual Review of Fluid Mechanics, 1996.

[5] Roshko A. On the development of turbulent wakes from vortex streets. NACA Report 1191, 1954.

[6] Zou Q and He X. On pressure and velocity boundary conditions for the lattice Boltzmann BGK model. Physics of Fluids, 1997.

[7] Succi S. The Lattice Boltzmann Equation for Fluid Dynamics and Beyond. Oxford University Press, 2001.

## Author

Saiyed Mudassir

MSc Aerospace Computational Engineering  
BEng Aircraft Propulsion and Turbomachinery  

Focus areas: CFD, aerodynamics, numerical methods and engineering design.
