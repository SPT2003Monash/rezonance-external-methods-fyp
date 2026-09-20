# Student A Data Dictionary

## Data categories

### `analytical_rlc`

Datasets used to validate the Richardson method against analytical RL and RLC admittance models.

### `average_gfmi/op1`

Average-GFMI datasets at the original operating point.

Filename pattern:

`AVG_GFMI_<frequency>Hz_<axis>pert_<timestep>us.csv`

Example:

`AVG_GFMI_5Hz_dpert_25us.csv`

### `average_gfmi/op2`

Average-GFMI datasets at operating point 2, with Pref = 0.6 pu.

Filename pattern:

`AVG_GFMI_OP2_<frequency>Hz_<axis>pert_<timestep>us.csv`

### `pulse`

Average-GFMI pulse-response datasets, accepted-frequency results and Rezonance comparison data.

### `switching_gfmi`

Switching-GFMI datasets used for raw waveform, envelope and Richardson analyses.

Filename pattern:

`Switching_GFMI_<frequency>Hz_<axis>pert_<timestep>.csv`

## Naming conventions

- `dpert` — d-axis voltage perturbation
- `qpert` — q-axis voltage perturbation
- `50us`, `25us`, etc. — PSCAD solver timestep
- `OP2` — second operating point
- `ref` — held-out reference case

## Main measured quantities

Depending on the dataset, columns include:

- simulation time
- three-phase PCC voltage
- three-phase PCC current
- active power
- reactive power
- RMS voltage
- dq-frame voltage and current quantities

The exact column names follow the PSCAD output-channel names used by each model.