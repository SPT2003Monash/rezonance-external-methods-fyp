# Reproducing Student A Results

Run the MATLAB scripts from their locations inside:

`Student_A_Final/scripts/student_a/`

The scripts determine the repository root from their own locations, so MATLAB does not need to be opened in a particular working directory.

## Recommended run order

### 1. Analytical RL and RLC validation

Run the scripts in:

`scripts/student_a/analytical_rlc/`

The required datasets are in:

`data/student_a/analytical_rlc/`

### 2. Average GFMI Richardson validation

Run:

`average_gfmi_richardson_analysis.m`

Set the required frequency and operating-point identifier near the beginning of the script.

Available operating points:

- OP1: original operating point Pref = 1.0 pu
- OP2: Pref = 0.6 pu

Principal frequencies:

- 1 Hz
- 5 Hz
- 10 Hz

The input datasets are located in:

- `data/student_a/average_gfmi/op1/`
- `data/student_a/average_gfmi/op2/`

### 3. Average GFMI Pareto analysis

Run:

`average_gfmi_pareto_analysis.m`

This generates the accuracy-runtime Pareto figures using the final 2 µs-reference errors.

### 4. Average GFMI Rezonance comparison

Run:

`average_gfmi_rezonance_comparison.m`

### 5. Pulse identification

Run:

`average_gfmi_pulse_identification.m`

Then run:

`average_gfmi_pulse_rezonance_comparison.m`

Pulse datasets are located in:

`data/student_a/pulse/`

### 6. Switching GFMI analysis

Run:

1. `switching_gfmi_raw_vs_envelope.m`
2. `switching_gfmi_richardson_analysis.m`
3. `switching_gfmi_pareto_analysis.m`

The switching datasets are located in:

`data/student_a/switching_gfmi/`

## Notes

- Do not rename the CSV files because the scripts use fixed filename patterns.
- Generated result tables are saved into the corresponding results folders.
- Generated figures are saved under `figures/student_a/`.
