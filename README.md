# Replication code — Targeted Local Projections (TLP)

Self-contained MATLAB code for the Monte Carlo simulations and the empirical
(BLP) application comparing inference methods for impulse response functions:
Local Projections (LP), Smooth LP (SLP), VAR, Bayesian LP (BLP), and the
**Targeted Local Projection (TLP)** — an optimal linear combination of LP and VAR.

This folder contains only the files required to run the three entry points below;
all unused/legacy code from the development repository has been removed.

## Entry points (at the folder root)

| Script | What it does |
|---|---|
| `MAIN_simulation.m` | Monte Carlo study under the Olea et al. Smets–Wouters DGP (data generation → LP/SLP/VAR/BLP/TLP estimation → bootstrap CIs → coverage/length/bias/RMSE statistics). Saves one result `.mat` per sample size into `SimulationsResults/`. |
| `ApplicationQuarterly.m` | Empirical application: 7-variable quarterly US VAR with a Cholesky-identified FFR shock, 1954Q3–2019Q4. Reads `DataNew2019.xlsx`; saves a `BLP_application_*.mat` into `SimulationsResults/`. |
| `PlotsPaper.m` | Produces every figure and table in the paper from the result `.mat` files in `SimulationsResults/`. |

## How to run

Each entry script resolves its own paths relative to its own location (via
`mfilename`), so you can run it from any working directory. Open the script and
run it:

```matlab
MAIN_simulation       % simulation study  -> SimulationsResults/<name_sim>*.mat
ApplicationQuarterly  % empirical app     -> SimulationsResults/BLP_application_*.mat
PlotsPaper            % all paper figures/tables (reads SimulationsResults/)
```

Each script runs `addpath(genpath('Subroutines'))` and adds `SimulationsResults`
to the path automatically. Figures and tables are written to a `TablesAndPlots/`
subfolder (created automatically in the current working directory).

## Directory layout

```
.
├── MAIN_simulation.m       simulation entry point (Olea SW DGP)
├── ApplicationQuarterly.m  empirical-application entry point (BLP / FFR shock)
├── PlotsPaper.m            figure/table driver
├── DataNew2019.xlsx        data for the empirical application
├── SimulationsResults/     result .mat files (written here, read by PlotsPaper)
└── Subroutines/
    ├── Shrink_to_VAR/        VAR estimation, IRFs, bootstrap CIs
    ├── Quadratic/            SLP smoothing (penalty, HAC, bootstrap)
    ├── Bayesian/subroutines/ Bayesian LP (csminwel optimizer, priors)
    ├── OleaDGP_SW/           Smets–Wouters DGP
    │   └── inputs/varma_sw_dgps_*.mat   precomputed DGP (loaded by run_Olea_DGP)
    └── other_functions/      TLP/SLP/VAR/BLP estimators, plotting, bootstrap, stats, helpers
```

## Required MATLAB toolboxes

- Statistics and Machine Learning Toolbox (`quantile`, `corr`, `nanmean`, `norminv`)
- Optimization Toolbox (`fmincon`)
- Econometrics Toolbox (`lagmatrix`)
- Financial Toolbox (`x2mdate`, used only by the empirical application for date parsing)
- Parallel Computing Toolbox (`parfor`) — optional; runs serially if unavailable.

## Reproducing the figures: the result `.mat` files

`PlotsPaper.m` does not run any estimation; it `load`s result `.mat` files from
`SimulationsResults/` and draws from them. The full set it expects is:

- `Final16June_T{200,800}_sim1000_pL10_pV8_eta{1,2,4,8,16,32,64}_DGP_OleaSW_P1.mat`
  — main simulations (`P_VAR = 8`, all `eta`).
- `Final16June_T{200,800}_sim1000_pL10_pV1_eta1_DGP_OleaSW_P1.mat`
  — `P_VAR = 1` comparison (only `eta = 1`).
- `Final16JuneGarch_T{200,800}_sim1000_pL10_pV8_eta{1,2,4}_DGP_OleaSW_P1.mat`
  — GARCH-innovation variants.
- `BLP_application_T259_pL{8,10}_pV4_FFRshock.mat`
  — empirical-application outputs.

Generate these by running `MAIN_simulation.m` and `ApplicationQuarterly.m`; the
filenames must match exactly. In `MAIN_simulation.m` the name stem and grid are set
near the top:

- `name_sim`  — filename prefix; set it to `"Final16June_T"` to match `PlotsPaper.m`.
- `eta_scale` — the `eta` grid to loop over (`[1 2 4 8 16 32 64]`).
- `P_VAR`     — set to `8` for the main runs and `1` for the `pV1` comparison.

The GARCH variants additionally require `garch_flag = 1` (and the `Garch` name
stem); the empirical-application files come from running `ApplicationQuarterly.m`
at `P_LP = 8` and `P_LP = 10`. These result files are large and are not shipped
with the code — generate or copy them into `SimulationsResults/` before running
`PlotsPaper.m`.
