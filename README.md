HIV Transmission Dynamics Agent-Based Model
================

- [HIV Transmission Dynamics Agent-Based
  Model](#hiv-transmission-dynamics-agent-based-model)
  - [Simulation overview](#simulation-overview)
  - [Running the Simulation](#running-the-simulation)
  - [Intervention Scenarios Summary](#intervention-scenarios-summary)
  - [Note](#note)

# HIV Transmission Dynamics Agent-Based Model

This MATLAB project simulates HIV transmission dynamics within a
population of men who have sex with men (MSM) in the Netherlands using
an agent-based model (ABM). The model includes two phases:
**calibration** and **main forward simulation**, with multiple scenarios
that allow exploring various intervention simulating increase in
diagnosis rates.

## Simulation overview

1.  **Calibration Phase**: This phase calibrates the model by to fit
    observed data in 2017-2023. Each calibration run is initialized by
    specifying the `pars_num` and `traj_num` parameters.
2.  **Main Forward Simulation Phase**: After calibration, this phase
    simulates intervention scenarios by altering testing rates and
    target populations, allowing for comparisons between baseline and
    intervention scenarios. The simulation starts in 2024, the
    intervention starts in 2026 and runs for 15 years.

## Running the Simulation

### 1. Calibration Phase

To start the calibration, use the following function:

``` matlab
Main_Traj_create_State(pars_num, traj_num)
```

- **`pars_num`**: Specifies the parameter number (runs from 1 to 100)
  used to load the configuration file.
- **`traj_num`**: Specifies the trajectory number (runs from 1 to 20),
  representing different trajectories to capture model variability.

**Example**:

``` matlab
Main_Traj_create_State(1, 1)
```

### 2. Main Forward Simulation Phase

Once calibration is complete, initiate the forward simulation with
various intervention scenarios using:

``` matlab
Main_Traj_Load(pars_num, traj_num, inter_str)
```

- **`pars_num`**: Same parameter number used in the calibration phase.
- **`traj_num`**: Same trajectory number used in the calibration phase.
- **`inter_str`**: A string specifying the intervention scenario.
  Available values:
  - `'base'`: Baseline scenario with no increase in diagnosis rate.
  - `'immigr_XX_YY'`: One-time testing of incoming immigrants at the
    point of entry, with `XX.YY%` of individuals agreeing to take the
    HIV test (values of `XX.YY`: 10, 25, 50, 75). Note that here and
    elsewhere, decimal values are specified by two integer components,
    where `XX` denotes the integer part and `YY` denotes the fractional
    remainder, so that the numeric value is interpreted as `XX + YY/100`
    (e.g., `25_0` corresponds to 25.0).
  - `'parts_inc_XX_YY_WW_ZZ'`: Increase in testing frequency to at least
    `WW.ZZ` per year for all individuals whose number of non-steady
    partners in the last 6 months was equal or exceeded `XX.YY` (values
    of `XX.YY`: 5, 10; values of `WW.ZZ`: 0.2, 0.53, 1.79).
  - `'late_inc_XX_YY_1_0_WW_ZZ'`: Decrease in testing interval to at
    least `WW.ZZ` years for all individuals who tested less frequently
    that this value (values of `WW.ZZ`: 0.55, 1.88, 5.0; value of
    `XX.YY` should be strictly less than the value of `WW.ZZ`).
  - `'immigr_parts_XX_YY_WW_ZZ_AA_BB'`: Combination of one-time testing
    of incoming immigrants at the point of entry, with `XX.YY%` agreeing
    to take the HIV test (values of `XX.YY`: 25, 50) with increase in
    testing frequency to at least `AA.BB` per year for all individuals
    whose number of non-steady partners in the last 6 months was equal
    or exceeded `WW.ZZ` (values of `XX.YY`: 5, 10; values of `WW.ZZ`:
    0.2, 0.53, 1.79).
  - `'immigr_late_XX_YY_WW_ZZ'`: Combination of one-time testing of
    incoming immigrants at the point of entry, with XX.YY% agreeing to
    take the HIV test (values of `XX.YY`: 25, 50) with decrease in
    testing interval to at least `WW.ZZ` years for all individuals who
    tested less frequently that this value (values of `WW.ZZ`: 0.55,
    1.88, 5.0).

**Example**:

``` matlab
Main_Traj_Load(1, 1, 'immigr_late_25_0_0_55')
```

### Output Structure

*Simulation outputs for calibration are saved in folders structured as
follows:*

    Output_par_num_pars_num_batch_traj_num

*Simulation outputs for main scenarios are saved in folders structured
as follows:*

    Output_par_num_pars_num_batch_traj_num_inter_str

Each folder name reflects the parameter settings, trajectory number, and
scenario used in that particular run, storing output data for each
specific intervention scenario.

## Intervention Scenarios Summary

- **`base`**: Baseline scenario, no additional diagnosis interventions.
- **`immigr_XX_YY`**: One-time testing of incoming immigrants at the
  point of entry
- **`parts_inc_XX_YY_WW_ZZ`**: Increase in testing frequency to at least
  `WW.ZZ` per year for all individuals whose number of non-steady
  partners in the last 6 months was equal or exceeded `XX.YY`
- **`late_inc_XX_YY_1_0_WW_ZZ`**: Decrease in testing interval to at
  least `WW.ZZ` years for all individuals who tested less frequently
  that this value.
- **`immigr_parts_XX_YY_WW_ZZ_AA_BB`**: Combination of one-time testing
  of incoming immigrants at the point of entry, with `XX.YY%` agreeing
  to take the HIV test with increase in testing frequency to at least
  AA.BB per year for all individuals whose number of non-steady partners
  in the last 6 months was equal or exceeded `WW.ZZ`
- **`immigr_late_XX_YY_WW_ZZ`**: Combination of one-time testing of
  incoming immigrants at the point of entry, with `XX.YY%` agreeing to
  take the HIV test with decrease in testing interval to at least
  `WW.ZZ` years for all individuals who tested less frequently that this
  value

By running various scenarios, this model provides insights into the
potential impact of increased diagnosis rates on HIV transmission
dynamics within the MSM population.

## Note

Calibration rounds creates a folder **State** with subfolders
**State_pars_XX_batch_YY**:

- **`XX`**: Parameter number used in the calibration phase.
- **`YY`**: Trajectory number used in the calibration phase.

These folders contain snapshot of the population at the end of the
calibration round, which is used to initialize the model in the main
scenario simulation.
