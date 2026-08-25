rm(list = ls())

source("R/simulation_functions.R")
source("R/nested_ZA_interaction_workflow.R")

message("Calibrating a small gamma_ZA = 0 scenario...")
pars <- make_nested_scenario_pars(
    p_scenario = 5,
    gamma_ZA = 0,
    calibration_mc_n = 5000,
    calibration_seed = 202606
)
check_nested_calibration(pars, tol = 0.01)

truth <- true_decompositions_mc_nested(
    pars = pars,
    mc_n = 5000,
    seed_top = 1001,
    seed_mod = 2001
)
check_nested_true_decompositions(truth)

set.seed(202607)
dat <- simulate_one_nested(n = 300, pars = pars)
check_nested_simulated_prob(dat)

message("Running two ordinary parametric replicates...")
result <- run_estimator_simulation_nested(
    B = 2,
    n = 300,
    pars = pars,
    true_vals = truth,
    fit_X_given_Z = TRUE,
    outcome_interaction = "none",
    seed = 202607,
    verbose = FALSE
)
stopifnot(
    nrow(result$estimates$topological) == 2L,
    nrow(result$estimates$modified_direct) == 2L,
    nrow(result$estimates$modified_bayes) == 2L
)
message("Smoke test passed. No files were written.")

