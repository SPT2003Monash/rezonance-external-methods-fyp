# Student A — Accelerated Impedance Identification

This folder contains the final files for Student A's contribution to the FYP on accelerated dq-admittance identification.

## Main investigations

- Richardson extrapolation for analytical RL/RLC cases
- Average GFMI validation at two operating points
- Switching GFMI validation
- Low-frequency pulse identification
- Accuracy-runtime Pareto analysis
- Comparison against Rezonance results

## Folder structure

- `cases/` — PSCAD model cases
- `data/student_a/` — input datasets and generated result tables
- `scripts/student_a/` — final MATLAB analysis scripts
- `figures/student_a/` — final figures
- `reports/` — final report and presentation material
- `docs/` — documentation and reproduction instructions

## Tested frequencies

The principal GFMI validation frequencies are 1, 5 and 10 Hz.

## Average GFMI timestep sequence

Richardson inputs:

- 50 µs
- 25 µs
- 12.5 µs

Held-out validation reference:

- 2 µs

## Switching GFMI timestep sequence

See `REPRODUCING_RESULTS.md` for the exact timestep sequence and scripts.

## Software

- PSCAD 5.0.1
- MATLAB