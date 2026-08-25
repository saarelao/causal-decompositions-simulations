rm(list = ls())

suppressPackageStartupMessages({
    library(MASS)
    library(nnet)
})

source("R/simulation_functions.R")
source("R/nested_ZA_interaction_workflow.R")

env_integer <- function(name, default) {
    value <- Sys.getenv(name, unset = "")
    if (!nzchar(value)) return(as.integer(default))
    parsed <- suppressWarnings(as.integer(value))
    if (is.na(parsed) || parsed < 1L) stop(name, " must be a positive integer.")
    parsed
}

env_numeric_vector <- function(name, default) {
    value <- Sys.getenv(name, unset = "")
    if (!nzchar(value)) return(default)
    parsed <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]))
    if (any(!is.finite(parsed))) stop(name, " must be a comma-separated numeric vector.")
    parsed
}

env_flag <- function(name, default) {
    value <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
    if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
        stop(name, " must be true or false.")
    }
    value %in% c("true", "1", "yes")
}

B <- env_integer("SIM_B", 500)
truth_mc_n <- env_integer("TRUTH_MC_N", 200000)
calibration_mc_n <- env_integer("CALIBRATION_MC_N", 200000)
base_seed <- env_integer("BASE_SEED", 202606)
calibration_tol <- 0.002
gamma0_truth_tol <- 0.003
nfold <- env_integer("XGB_NFOLD", 5)
max_nrounds <- env_integer("XGB_MAX_NROUNDS", 2000)
early_stopping_rounds <- env_integer("XGB_EARLY_STOPPING", 50)
out_dir <- Sys.getenv(
    "POINT_RESULTS_DIR",
    unset = "outputs/simulation_results_nested_ZA_interaction"
)
skip_existing <- env_flag("SKIP_EXISTING", TRUE)
n_values <- env_numeric_vector("SIM_N", c(500, 1500))
p_values <- env_numeric_vector("SIM_P", c(5, 15))
gamma_values <- env_numeric_vector("SIM_GAMMA", c(0, 0.6, 1.0))

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(file.path(out_dir, "csv"))) {
    dir.create(file.path(out_dir, "csv"), recursive = TRUE)
}

scenario_grid <- expand.grid(
    n = n_values,
    p = p_values,
    gamma_ZA = gamma_values,
    KEEP.OUT.ATTRS = FALSE
)
scenario_grid <- scenario_grid[order(scenario_grid$p, scenario_grid$gamma_ZA, scenario_grid$n), ]
scenario_grid$seed <- base_seed + seq_len(nrow(scenario_grid))

xgb_grid_name <- "expanded"
xgb_grid <- default_xgb_grid_expanded()

metadata <- list(
    timestamp = as.character(Sys.time()),
    B = B,
    truth_mc_n = truth_mc_n,
    calibration_mc_n = calibration_mc_n,
    base_seed = base_seed,
    calibration_tol = calibration_tol,
    gamma0_truth_tol = gamma0_truth_tol,
    xgb_grid_name = xgb_grid_name,
    xgb_grid = xgb_grid,
    nfold = nfold,
    max_nrounds = max_nrounds,
    early_stopping_rounds = early_stopping_rounds,
    skip_existing = skip_existing,
    scenarios = scenario_grid
)

saveRDS(metadata, file = file.path(out_dir, "run_metadata.rds"))
write.csv(xgb_grid, file = file.path(out_dir, "xgb_grid_expanded.csv"), row.names = FALSE)
writeLines(
    c(
        paste0("timestamp: ", metadata$timestamp),
        paste0("B: ", B),
        paste0("truth_mc_n: ", truth_mc_n),
        paste0("calibration_mc_n: ", calibration_mc_n),
        paste0("xgb_grid_rows: ", nrow(xgb_grid)),
        paste0("max_nrounds: ", max_nrounds),
        paste0("early_stopping_rounds: ", early_stopping_rounds),
        paste0("nfold: ", nfold),
        "scenarios:",
        paste0(
            "  n=",
            scenario_grid$n,
            ", p=",
            scenario_grid$p,
            ", gamma_ZA=",
            scenario_grid$gamma_ZA,
            ", seed=",
            scenario_grid$seed
        )
    ),
    con = file.path(out_dir, "run_metadata.txt")
)

error_log <- file.path(out_dir, "error_log.txt")
if (file.exists(error_log)) {
    file.remove(error_log)
}

pars_cache <- list()
truth_cache <- list()

for (i in seq_len(nrow(scenario_grid))) {
    n_scenario <- scenario_grid$n[i]
    p_scenario <- scenario_grid$p[i]
    gamma_ZA <- scenario_grid$gamma_ZA[i]
    scenario_seed <- scenario_grid$seed[i]
    gamma_id <- format_gamma_ZA_for_id(gamma_ZA)
    cache_id <- paste(p_scenario, gamma_id, sep = "_")

    if (is.null(pars_cache[[cache_id]])) {
        cat("\n========================================\n")
        cat("Preparing nested scenario p = ", p_scenario, ", gamma_ZA = ", gamma_ZA, "\n", sep = "")
        cat("========================================\n")

        pars_nested <- make_nested_scenario_pars(
            p_scenario = p_scenario,
            gamma_ZA = gamma_ZA,
            calibration_mc_n = calibration_mc_n,
            calibration_seed = base_seed + 1000 + length(pars_cache)
        )
        check_nested_calibration(pars_nested, tol = calibration_tol)

        true_vals <- true_decompositions_mc_nested(
            pars = pars_nested,
            mc_n = truth_mc_n
        )
        check_nested_true_decompositions(true_vals)

        if (gamma_ZA == 0) {
            check_nested_gamma0_matches_original(
                pars_nested = pars_nested,
                mc_n = truth_mc_n,
                tol = gamma0_truth_tol
            )
        }

        pars_cache[[cache_id]] <- pars_nested
        truth_cache[[cache_id]] <- true_vals
    }

    pars_nested <- pars_cache[[cache_id]]
    true_vals <- truth_cache[[cache_id]]

    dat_check <- simulate_one_nested(n = 1000, pars = pars_nested)
    check_nested_simulated_prob(dat_check)

    param_specs <- list(
        no_int = "none",
        ZA_int = "ZA"
    )

    for (spec_name in names(param_specs)) {
        outcome_interaction <- param_specs[[spec_name]]
        scenario_id <- paste("nested", paste0("n", n_scenario), paste0("p", p_scenario), gamma_id, "param", spec_name, sep = "_")

        cat("\n========================================\n")
        cat("Running ", scenario_id, "\n", sep = "")
        cat("========================================\n")

        tryCatch({
            if (skip_existing && file.exists(file.path(out_dir, paste0(scenario_id, ".rds")))) {
                message("Skipping existing scenario: ", scenario_id)
            } else {
                sim_res <- run_estimator_simulation_nested(
                    B = B,
                    n = n_scenario,
                    pars = pars_nested,
                    true_vals = true_vals,
                    fit_X_given_Z = TRUE,
                    outcome_interaction = outcome_interaction,
                    seed = scenario_seed,
                    verbose = TRUE
                )
                save_simulation_result_nested(
                    sim_res = sim_res,
                    pars = pars_nested,
                    scenario_id = scenario_id,
                    out_dir = out_dir,
                    save_raw = TRUE,
                    save_csv = TRUE
                )
            }
        }, error = function(err) {
            message("Scenario failed: ", scenario_id)
            message(conditionMessage(err))
            log_nested_scenario_error(error_log, scenario_id, err)
        })
    }

    scenario_id <- paste("nested", paste0("n", n_scenario), paste0("p", p_scenario), gamma_id, "param_ZA_firth", sep = "_")

    cat("\n========================================\n")
    cat("Running ", scenario_id, "\n", sep = "")
    cat("========================================\n")

    tryCatch({
        if (skip_existing && file.exists(file.path(out_dir, paste0(scenario_id, ".rds")))) {
            message("Skipping existing scenario: ", scenario_id)
        } else {
            sim_res_firth <- run_estimator_simulation_nested_firth(
                B = B,
                n = n_scenario,
                pars = pars_nested,
                true_vals = true_vals,
                fit_X_given_Z = TRUE,
                seed = scenario_seed,
                verbose = TRUE
            )
            save_simulation_result_nested(
                sim_res = sim_res_firth,
                pars = pars_nested,
                scenario_id = scenario_id,
                out_dir = out_dir,
                save_raw = TRUE,
                save_csv = TRUE
            )
        }
    }, error = function(err) {
        message("Scenario failed: ", scenario_id)
        message(conditionMessage(err))
        log_nested_scenario_error(error_log, scenario_id, err)
    })

    scenario_id <- paste("nested", paste0("n", n_scenario), paste0("p", p_scenario), gamma_id, "xgb", xgb_grid_name, sep = "_")

    cat("\n========================================\n")
    cat("Running ", scenario_id, "\n", sep = "")
    cat("========================================\n")

    tryCatch({
        if (skip_existing && file.exists(file.path(out_dir, paste0(scenario_id, ".rds")))) {
            message("Skipping existing scenario: ", scenario_id)
        } else {
            sim_res_xgb <- run_estimator_simulation_with_xgb_nested(
                B = B,
                n = n_scenario,
                pars = pars_nested,
                true_vals = true_vals,
                xgb_grid = xgb_grid,
                nfold = nfold,
                max_nrounds = max_nrounds,
                early_stopping_rounds = early_stopping_rounds,
                seed = scenario_seed,
                verbose = TRUE,
                xgb_grid_name = xgb_grid_name
            )
            save_simulation_result_nested(
                sim_res = sim_res_xgb,
                pars = pars_nested,
                scenario_id = scenario_id,
                out_dir = out_dir,
                save_raw = TRUE,
                save_csv = TRUE
            )
        }
    }, error = function(err) {
        message("Scenario failed: ", scenario_id)
        message(conditionMessage(err))
        log_nested_scenario_error(error_log, scenario_id, err)
    })
}

cat("\nNested Z:A interaction simulation run complete.\n")
