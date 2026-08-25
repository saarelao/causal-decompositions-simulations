# Simulation code for graph-based causal variance decompositions

This repository reproduces the simulation results reported in the manuscript
"Graph-based causal variance decompositions: When 'variance explained' means causation" 
dated 2026-08-18. It contains the data-generating mechanism, estimators,
Monte Carlo workflows, and figure code. No individual-level or external data
are required.

## What is reproduced

The simulation varies:

- sample size: `n = 500, 1500`;
- number of baseline covariates: `p = 5, 15`;
- strength of the `Z:A` outcome interaction: `gamma_ZA = 0, 0.6, 1`.

The point-estimation study evaluates the ordinary parametric estimators, the
`Z:A`-interaction parametric estimator, its Firth-logistic variant, and XGBoost.
For modified ordering, both direct and Bayes-inversion versions are evaluated
where applicable. The uncertainty study evaluates ordinary parametric models
with and without the `Z:A` interaction under topological ordering, modified
direct ordering, and modified Bayes inversion.

The workflows generate:

- the complete point-estimation figures used in the supplement;
- the simplified point-estimation figures used in the main manuscript;
- the main-manuscript coverage figures;
- the full coverage and empirical-SD/posterior-SD diagnostic figures.

The main manuscript displays the `gamma_ZA = 0` and `gamma_ZA = 1` simplified
figures. The `gamma_ZA = 0.6` results and all full diagnostics are generated
because they appear in the supplement.

## Repository layout

```text
R/             Simulation, estimation, and plotting functions
scripts/       Executable workflows, in their intended run order
environment/   Package versions from the reference analysis environment
outputs/       Generated results (ignored by Git, except for this directory)
```

Run every command from the repository root.

## Software

The reference analysis used R 4.6.1. Package versions are recorded in
`environment/package-versions.csv`. Install the required packages with:

```r
install.packages(c(
  "ggplot2", "ggh4x", "logistf", "MASS", "Matrix",
  "nnet", "scales", "xgboost"
))
```

Then check the environment:

```sh
Rscript scripts/00_check_environment.R
```

The check reports version differences but stops only when a required package is
missing. For archival reproduction, use the versions in the CSV file or create
an `renv` lockfile on the target platform before running the simulations.

## Quick validation

The smoke test calibrates a small data-generating mechanism, checks both true
decompositions, simulates data, and runs two ordinary parametric replicates. It
does not write output files.

```sh
Rscript scripts/01_smoke_test.R
```

## Full reproduction

The following commands reproduce all simulation results and manuscript
artifacts:

```sh
Rscript scripts/02_run_point_simulations.R
Rscript scripts/03_run_uncertainty_simulations.R
Rscript scripts/04_make_figures.R
```

The first two commands are computationally expensive. In particular, the point
study repeatedly tunes XGBoost by five-fold cross-validation, and the
uncertainty study uses 500 posterior draws inside each of 500 simulation
replicates. The uncertainty workflow defaults to four parallel workers. Plan to
run these workflows on a suitable compute system rather than interactively.

Results are saved as RDS files, with CSV summaries where applicable. By
default, completed scenario files are retained and skipped on a subsequent run,
so an interrupted workflow can be resumed safely. Generated outputs are not
tracked by Git.

## Reference settings and seeds

The full-study defaults are defined at the top of the two simulation scripts:

| Setting | Point study | Uncertainty study |
|---|---:|---:|
| Monte Carlo replicates | 500 | 500 |
| Posterior draws per replicate | — | 500 |
| Truth Monte Carlo sample | 200,000 | 200,000 |
| Calibration Monte Carlo sample | 200,000 | 200,000 |
| Base seed | 202606 | 202606 |
| XGBoost folds | 5 | — |
| XGBoost maximum rounds | 2,000 | — |
| XGBoost early stopping | 50 rounds | — |

Each scenario receives a deterministic seed derived from the base seed. The
saved metadata files record the settings, scenario grid, and seeds used by a
run.

## Optional trial-run settings

Environment variables can reduce the workload or select scenarios without
editing the scripts:

| Variable | Default | Meaning |
|---|---|---|
| `SIM_B` | `500` | Monte Carlo replicates |
| `B_POST` | `500` | Posterior draws in the uncertainty study |
| `TRUTH_MC_N` | `200000` | Monte Carlo sample for the simulation truth |
| `CALIBRATION_MC_N` | `200000` | Monte Carlo sample for calibration |
| `SIM_N` | `500,1500` | Comma-separated sample sizes |
| `SIM_P` | `5,15` | Comma-separated covariate dimensions |
| `SIM_GAMMA` | `0,0.6,1` | Comma-separated interaction values |
| `N_WORKERS` | `4` | Parallel uncertainty workers |
| `SKIP_EXISTING` | `true` | Skip an existing scenario RDS file |
| `POINT_RESULTS_DIR` | default output path | Point-result directory |
| `UNCERTAINTY_RESULTS_DIR` | default output path | Uncertainty-result directory |
| `FIGURE_DIR` | default output path | Figure directory |

The point script also accepts `XGB_NFOLD`, `XGB_MAX_NROUNDS`, and
`XGB_EARLY_STOPPING`. Nondefault settings are intended for code checks; results
for the manuscript require the defaults above.

For example, this PowerShell session runs one reduced point-estimation scenario
and keeps it separate from full results:

```powershell
$env:SIM_B = "2"
$env:TRUTH_MC_N = "5000"
$env:CALIBRATION_MC_N = "5000"
$env:SIM_N = "500"
$env:SIM_P = "5"
$env:SIM_GAMMA = "0"
$env:XGB_NFOLD = "2"
$env:XGB_MAX_NROUNDS = "20"
$env:XGB_EARLY_STOPPING = "5"
$env:POINT_RESULTS_DIR = "outputs/trial_point_results"
Rscript scripts/02_run_point_simulations.R
```

Open a new shell, or remove these environment variables, before running the
full study.

## Citation

If you use this code in academic work, please cite the accompanying paper.

## License

This code is released under the MIT License. See the `LICENSE` file for details.
