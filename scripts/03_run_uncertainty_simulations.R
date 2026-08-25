rm(list = ls())

suppressPackageStartupMessages({
    library(MASS)
    library(nnet)
})

source("R/simulation_functions.R")
source("R/nested_ZA_interaction_workflow.R")
source("R/uncertainty_simulation_nested.R")

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
B_post <- env_integer("B_POST", 500)
truth_mc_n <- env_integer("TRUTH_MC_N", 200000)
calibration_mc_n <- env_integer("CALIBRATION_MC_N", 200000)
base_seed <- env_integer("BASE_SEED", 202606)
posterior_seed_offset <- 1000000L
n_workers <- env_integer("N_WORKERS", 4)
calibration_tol <- 0.002
gamma0_truth_tol <- 0.003
out_dir <- Sys.getenv(
    "UNCERTAINTY_RESULTS_DIR",
    unset = "outputs/uncertainty_results_nested_ZA_interaction"
)
skip_existing <- env_flag("SKIP_EXISTING", TRUE)
n_values <- env_numeric_vector("SIM_N", c(500, 1500))
p_values <- env_numeric_vector("SIM_P", c(5, 15))
gamma_values <- env_numeric_vector("SIM_GAMMA", c(0, 0.6, 1.0))

scenario_grid <- expand.grid(
    n = n_values,
    p = p_values,
    gamma_ZA = gamma_values,
    KEEP.OUT.ATTRS = FALSE
)
scenario_grid <- scenario_grid[
    order(scenario_grid$p, scenario_grid$gamma_ZA, scenario_grid$n),
    ,
    drop = FALSE
]
scenario_grid$seed <- base_seed + seq_len(nrow(scenario_grid))
scenario_indices <- seq_len(nrow(scenario_grid))

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(file.path(out_dir, "csv"))) {
    dir.create(file.path(out_dir, "csv"), recursive = TRUE)
}

metadata <- list(
    timestamp = as.character(Sys.time()),
    B = B,
    B_post = B_post,
    truth_mc_n = truth_mc_n,
    calibration_mc_n = calibration_mc_n,
    base_seed = base_seed,
    posterior_seed_offset = posterior_seed_offset,
    n_workers = n_workers,
    calibration_tol = calibration_tol,
    gamma0_truth_tol = gamma0_truth_tol,
    skip_existing = skip_existing,
    outcome_interactions = c("none", "ZA"),
    scenarios = scenario_grid,
    scenario_indices = scenario_indices
)
saveRDS(metadata, file = file.path(out_dir, "run_metadata.rds"))
writeLines(
    c(
        paste0("timestamp: ", metadata$timestamp),
        paste0("B: ", B),
        paste0("B_post: ", B_post),
        paste0("truth_mc_n: ", truth_mc_n),
        paste0("calibration_mc_n: ", calibration_mc_n),
        paste0("posterior_seed_offset: ", posterior_seed_offset),
        paste0("n_workers: ", n_workers),
        paste0("outcome_interactions: ", paste(metadata$outcome_interactions, collapse = ", ")),
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
pars_cache <- list()
truth_cache <- list()

selected_scenarios <- scenario_grid[scenario_indices, , drop = FALSE]
job_grid <- do.call(rbind, lapply(seq_len(nrow(selected_scenarios)), function(i) {
    scenario <- selected_scenarios[i, , drop = FALSE]
    data.frame(
        n = scenario$n,
        p = scenario$p,
        gamma_ZA = scenario$gamma_ZA,
        seed = scenario$seed,
        outcome_interaction = c("none", "ZA"),
        stringsAsFactors = FALSE
    )
}))
job_grid$cache_id <- paste(
    job_grid$p,
    vapply(job_grid$gamma_ZA, format_gamma_ZA_for_id, character(1)),
    sep = "_"
)
job_grid$scenario_id <- vapply(seq_len(nrow(job_grid)), function(i) {
    make_nested_uncertainty_scenario_id(
        n = job_grid$n[i],
        p = job_grid$p[i],
        gamma_ZA = job_grid$gamma_ZA[i],
        outcome_interaction = job_grid$outcome_interaction[i]
    )
}, character(1))
job_grid$result_file <- file.path(out_dir, paste0(job_grid$scenario_id, ".rds"))
job_grid$exists <- file.exists(job_grid$result_file)
pending_grid <- if (skip_existing) {
    job_grid[!job_grid$exists, , drop = FALSE]
} else {
    job_grid
}

message(
    "Selected ",
    nrow(job_grid),
    " jobs: ",
    sum(job_grid$exists),
    " existing and ",
    nrow(pending_grid),
    " pending."
)

calibration_grid <- unique(
    scenario_grid[, c("p", "gamma_ZA"), drop = FALSE]
)
calibration_grid$cache_id <- paste(
    calibration_grid$p,
    vapply(calibration_grid$gamma_ZA, format_gamma_ZA_for_id, character(1)),
    sep = "_"
)
calibration_grid$calibration_seed <-
    base_seed + 1000 + seq_len(nrow(calibration_grid)) - 1

for (cache_id in unique(pending_grid$cache_id)) {
    calibration_row <- calibration_grid[
        calibration_grid$cache_id == cache_id,
        ,
        drop = FALSE
    ]
    p_scenario <- calibration_row$p[1]
    gamma_ZA <- calibration_row$gamma_ZA[1]
    cat("\n========================================\n")
    cat(
        "Preparing pending uncertainty family p = ",
        p_scenario,
        ", gamma_ZA = ",
        gamma_ZA,
        "\n",
        sep = ""
    )
    cat("========================================\n")

    pars_nested <- make_nested_scenario_pars(
        p_scenario = p_scenario,
        gamma_ZA = gamma_ZA,
        calibration_mc_n = calibration_mc_n,
        calibration_seed = calibration_row$calibration_seed[1]
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
    set.seed(base_seed + calibration_row$calibration_seed[1])
    dat_check <- simulate_one_nested(n = 1000, pars = pars_nested)
    check_nested_simulated_prob(dat_check)
    pars_cache[[cache_id]] <- pars_nested
    truth_cache[[cache_id]] <- true_vals
}

jobs <- lapply(seq_len(nrow(pending_grid)), function(i) {
    cache_id <- pending_grid$cache_id[i]
    list(
        n = pending_grid$n[i],
        p = pending_grid$p[i],
        gamma_ZA = pending_grid$gamma_ZA[i],
        seed = pending_grid$seed[i],
        outcome_interaction = pending_grid$outcome_interaction[i],
        pars = pars_cache[[cache_id]],
        true_vals = truth_cache[[cache_id]]
    )
})

job_status <- list()
if (length(jobs) > 0) {
    workers_used <- min(n_workers, length(jobs))
    message("Running ", length(jobs), " pending jobs on ", workers_used, " workers.")
    if (workers_used == 1) {
        job_status <- lapply(jobs, run_nested_uncertainty_job,
            B = B,
            B_post = B_post,
            posterior_seed_offset = posterior_seed_offset,
            out_dir = out_dir,
            skip_existing = skip_existing,
            save_csv = TRUE,
            worker_messages = TRUE
        )
    } else {
        cluster <- parallel::makeCluster(workers_used, outfile = "")
        project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
        parallel::clusterCall(cluster, function(path) {
            setwd(path)
            source("R/simulation_functions.R")
            source("R/uncertainty_simulation_nested.R")
            NULL
        }, project_dir)
        job_status <- tryCatch(
            parallel::parLapplyLB(
                cluster,
                jobs,
                function(job, ...) run_nested_uncertainty_job(job, ...),
                B = B,
                B_post = B_post,
                posterior_seed_offset = posterior_seed_offset,
                out_dir = out_dir,
                skip_existing = skip_existing,
                save_csv = TRUE,
                worker_messages = TRUE
            ),
            finally = parallel::stopCluster(cluster)
        )
    }
}

status_rows <- lapply(job_status, function(status) {
    as.data.frame(status, stringsAsFactors = FALSE)
})
if (skip_existing && any(job_grid$exists)) {
    status_rows <- c(
        lapply(job_grid$scenario_id[job_grid$exists], function(scenario_id) {
            data.frame(
                scenario_id = scenario_id,
                status = "existing",
                error = "",
                n_success = NA_integer_,
                n_failed = NA_integer_,
                worker_pid = NA_integer_,
                stringsAsFactors = FALSE
            )
        }),
        status_rows
    )
}
if (length(status_rows) > 0) {
    status_df <- do.call(rbind, status_rows)
    write.csv(
        status_df,
        file = file.path(out_dir, "nested_uncertainty_parallel_job_status.csv"),
        row.names = FALSE
    )
    failed <- status_df[status_df$status == "failed", , drop = FALSE]
    if (nrow(failed) > 0) {
        for (i in seq_len(nrow(failed))) {
            log_nested_scenario_error(
                error_log,
                failed$scenario_id[i],
                simpleError(failed$error[i])
            )
        }
        warning(nrow(failed), " parallel uncertainty jobs failed; see ", error_log)
    }
}

result_files <- list.files(
    out_dir,
    pattern = "^nested_uncertainty_.*[.]rds$",
    full.names = TRUE
)
if (length(result_files) > 0) {
    combined_summary <- make_nested_uncertainty_summary_data(result_files)
    write.csv(
        combined_summary,
        file = file.path(out_dir, "nested_uncertainty_combined_summary.csv"),
        row.names = FALSE
    )
}

cat("\nNested Z:A interaction uncertainty simulation run complete.\n")
