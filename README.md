# 2D LBM Wake Shape Study

Circle, square and droplet body comparison at Re = 200 using a C++ D2Q9 TRT Lattice Boltzmann solver.

This repository shows my numerical simulation and post processing work for a simple aerodynamic shape comparison. The aim was not to produce a final validated design, but to show the full simulation chain: solver setup, geometry control, force extraction, wake visualisation, time averaging and engineering judgement.

## What this study does

Three bodies are compared under the same flow condition:

* Circle
* Square
* Droplet

The same domain, Reynolds number, inlet condition and solver settings are used. Only the body shape is changed.

The comparison focuses on:

* Wake structure
* Vorticity shedding
* Pressure coefficient field
* Speed field
* Mean drag coefficient
* RMS lift coefficient
* Strouhal estimate
* Time averaged wake deficit
* Wake profile comparison

## Solver

The flow is solved using a two dimensional D2Q9 Lattice Boltzmann method with a two relaxation time collision model.

Main setup:

```text
Method: D2Q9 TRT LBM
Reynolds number: 200
Domain: 1000 x 320 lattice cells
Body diameter: 40 lattice cells
Inlet speed: 0.06 lattice units
Kinematic viscosity: 0.012
Warm up steps: 30000
Recording steps: 60000
Snapshots: 300 per case
Force samples: 6000 per case
```

## Why these shapes

The three shapes give a clean aerodynamic comparison.

The square has sharp corners, so the flow separates early and forms a wide wake.

The circle is smoother, but it is still a bluff body with a strong vortex street.

The droplet gives the flow a better recovery path behind the body, reducing wake width and unsteady loading.

## Main results

From the final recorded part of the simulation:

| Body    | Mean Cd | RMS Cl | Strouhal estimate |
| ------- | ------: | -----: | ----------------: |
| Circle  |   1.506 |  0.542 |             0.210 |
| Square  |   1.698 |  0.448 |             0.174 |
| Droplet |   1.181 |  0.229 |             0.203 |

The droplet reduced mean drag by about:

```text
30.5 percent compared with the square
21.6 percent compared with the circle
```

The droplet also gave the lowest lift fluctuation, which means the unsteady side force from vortex shedding was weaker.

## Engineering judgement

The square gives the poorest behaviour because its sharp corners force early separation and poor pressure recovery.

The circle is cleaner than the square, but the wake is still broad.

The droplet gives the best result in this setup because it narrows the wake, lowers drag and reduces lift fluctuation.

This is the main point of the study: do not only run a simulation. Read the wake, check the force history and judge what improved.

## Generated post processing

The Python post processing creates:

```text
geometry_check.png
vorticity_hero.png
pressure_coefficient.png
speed_field.png
force_coefficients.png
engineering_summary.png
wake_profiles.png
mean_wake.png
wake_evolution.mp4
results_summary.csv
results_summary.txt
```

## Repository structure

```text
2d-lbm-wake-shape-study/
│
├── README.md
├── main.cpp
├── Post_processing.py
├── requirements.txt
├── .gitignore
│
├── figures/
│   ├── geometry_check.png
│   ├── vorticity_hero.png
│   ├── force_coefficients.png
│   ├── engineering_summary.png
│   └── wake_evolution.mp4
│
├── results/
│   ├── results_summary.csv
│   └── results_summary.txt
│
└── docs/
    └── notes.md
```

The raw snapshot folders are not included because they are large. The solver can regenerate them.

## How to build the solver

This was developed in Visual Studio on Windows.

From Developer PowerShell:

```powershell
msbuild .\2D_LBM_Unsteady.sln /t:Rebuild /p:Configuration=Debug /p:Platform=x64
```

Then run:

```powershell
.\x64\Debug\2D_LBM_Unsteady.exe
```

## How to run post processing

From the project root:

```powershell
python Post_processing.py
```

Python packages:

```text
numpy
pandas
matplotlib
```

Install them with:

```powershell
pip install -r requirements.txt
```

## Notes on validation

This is a numerical simulation and visualisation study. It is not yet a final validation grade result.

The next technical checks would be:

* Grid independence
* Domain sensitivity
* Outlet and blockage sensitivity
* Force extraction check
* Comparison with benchmark circular cylinder data at Re = 200

## Author

Saiyed Mudassir

MSc Aerospace Computational Engineering
BEng Aircraft Propulsion and Turbomachinery

Focus areas: CFD, aerodynamics and numerical methods
