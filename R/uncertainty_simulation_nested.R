bind_successful_uncertainty_rows <- function(rows, estimator) {
    keep <- !vapply(rows, is.null, logical(1))
    if (!any(keep)) {
        stop("No successful uncertainty replicates for estimator ", estimator, ".")
    }
    do.call(rbind, rows[keep])
}

summarize_uncertainty_performance_nested <- function(interval_df, B_requested) {
    out <- summarize_uncertainty_performance(interval_df)
    out$n_success <- length(unique(interval_df$replicate))
    out$B_requested <- B_requested
    out
}

run_uncertainty_simulation_nested <- function(
    B = 200,
    n = 1000,
    pars,
    true_vals = NULL,
    fit_X_given_Z = TRUE,
    outcome_interaction = c("none", "ZA"),
    B_post = 1000,
    seed = 1,
    posterior_seed_offset = 1000000L,
    verbose = TRUE
) {
    outcome_interaction <- match.arg(outcome_interaction)
    if (is.null(true_vals)) {
        true_vals <- true_decompositions_mc_nested(pars = pars, mc_n = 1e+05)
    }
    truth_top <- true_vals$topological$components[
        c("delta_Z", "delta_X", "delta_A", "delta_res", "total")
    ]
    truth_mod <- true_vals$modified$components[
        c("delta_X", "delta_Z", "delta_A", "delta_res", "total")
    ]
    top_rows <- vector("list", B)
    mod_direct_rows <- vector("list", B)
    mod_bayes_rows <- if (fit_X_given_Z) vector("list", B) else NULL
    simulation_seeds <- seed + seq_len(B)
    posterior_seeds <- seed + posterior_seed_offset + seq_len(B)
    model_diagnostics <- vector("list", B)

    for (b in seq_len(B)) {
        set.seed(simulation_seeds[b])
        dat_b <- simulate_one_nested(n = n, pars = pars)
        fits_b <- fit_decomp_models_nested(
            dat = dat_b,
            pars = pars,
            fit_X_given_Z = fit_X_given_Z,
            outcome_interaction = outcome_interaction
        )
        diag_b <- model_diagnostics_nested(
            dat = dat_b,
            fits = fits_b,
            pars = pars,
            replicate = b,
            outcome_interaction = outcome_interaction,
            simulation_seed = simulation_seeds[b]
        )
        diag_b$posterior_seed <- posterior_seeds[b]
        diag_b$uncertainty_ok <- FALSE
        diag_b$uncertainty_error <- ""

        error_message <- tryCatch({
            if (!inherits(fits_b$fit_Z_given_X, "multinom") ||
                !inherits(fits_b$fit_A_given_XZ, "multinom") ||
                !inherits(fits_b$fit_Y_given_XZA, "glm")) {
                stop("One or more required parametric models failed to fit.")
            }
            set.seed(posterior_seeds[b])
            top_hat <- estimate_topological_decomp(dat_b, fits_b, pars)
            mod_direct_hat <- estimate_modified_decomp_direct(dat_b, fits_b, pars)
            top_unc <- estimate_decomposition_uncertainty(
                dat_b,
                pars,
                fits_b,
                estimator = "topological",
                B_post = B_post
            )
            mod_direct_unc <- estimate_decomposition_uncertainty(
                dat_b,
                pars,
                fits_b,
                estimator = "modified_direct",
                B_post = B_post
            )
            top_row_b <- make_uncertainty_interval_df(
                b,
                "topological",
                top_hat,
                top_unc,
                truth_top
            )
            mod_direct_row_b <- make_uncertainty_interval_df(
                b,
                "modified_direct",
                mod_direct_hat,
                mod_direct_unc,
                truth_mod
            )
            if (fit_X_given_Z) {
                mod_bayes_hat <- estimate_modified_decomp_bayes(dat_b, fits_b, pars)
                mod_bayes_unc <- estimate_decomposition_uncertainty(
                    dat_b,
                    pars,
                    fits_b,
                    estimator = "modified_bayes",
                    B_post = B_post
                )
                mod_bayes_row_b <- make_uncertainty_interval_df(
                    b,
                    "modified_bayes",
                    mod_bayes_hat,
                    mod_bayes_unc,
                    truth_mod
                )
            }
            top_rows[[b]] <- top_row_b
            mod_direct_rows[[b]] <- mod_direct_row_b
            if (fit_X_given_Z) {
                mod_bayes_rows[[b]] <- mod_bayes_row_b
            }
            NULL
        }, error = function(err) {
            conditionMessage(err)
        })
        if (is.null(error_message)) {
            diag_b$uncertainty_ok <- TRUE
        } else {
            diag_b$uncertainty_error <- error_message
        }
        model_diagnostics[[b]] <- diag_b

        if (verbose && (b %% max(1, floor(B / 10)) == 0)) {
            message(
                "Completed nested uncertainty replicate ",
                b,
                " / ",
                B,
                " (outcome_interaction = ",
                outcome_interaction,
                ")"
            )
        }
    }

    top_interval_df <- bind_successful_uncertainty_rows(
        top_rows,
        "topological"
    )
    mod_direct_interval_df <- bind_successful_uncertainty_rows(
        mod_direct_rows,
        "modified_direct"
    )
    mod_bayes_interval_df <- if (fit_X_given_Z) {
        bind_successful_uncertainty_rows(mod_bayes_rows, "modified_bayes")
    } else {
        NULL
    }
    diagnostics <- do.call(rbind, model_diagnostics)
    successful_replicates <- diagnostics$replicate[diagnostics$uncertainty_ok]

    list(
        params = c(
            pars,
            list(
                n = n,
                outcome_interaction = outcome_interaction,
                B = B,
                B_post = B_post,
                seed = seed,
                posterior_seed_offset = posterior_seed_offset
            )
        ),
        truth = list(topological = truth_top, modified = truth_mod),
        intervals = list(
            topological = top_interval_df,
            modified_direct = mod_direct_interval_df,
            modified_bayes = mod_bayes_interval_df
        ),
        summary = list(
            topological = summarize_uncertainty_performance_nested(
                top_interval_df,
                B
            ),
            modified_direct = summarize_uncertainty_performance_nested(
                mod_direct_interval_df,
                B
            ),
            modified_bayes = if (!is.null(mod_bayes_interval_df)) {
                summarize_uncertainty_performance_nested(
                    mod_bayes_interval_df,
                    B
                )
            } else {
                NULL
            }
        ),
        convergence = diagnostics[
            ,
            c("conv_Z", "conv_A", "conv_Y"),
            drop = FALSE
        ],
        simulation_seeds = data.frame(
            replicate = seq_len(B),
            simulation_seed = simulation_seeds,
            posterior_seed = posterior_seeds
        ),
        successful_replicates = successful_replicates,
        model_diagnostics = diagnostics,
        nested_metadata = nested_calibration_metadata(
            pars,
            n = n,
            outcome_interaction = outcome_interaction
        ),
        uncertainty_metadata = list(
            B = B,
            B_post = B_post,
            seed = seed,
            posterior_seed_offset = posterior_seed_offset,
            n_success = length(successful_replicates),
            n_failed = B - length(successful_replicates)
        )
    )
}

make_nested_uncertainty_scenario_id <- function(n,
                                                p,
                                                gamma_ZA,
                                                outcome_interaction) {
    spec_id <- if (identical(outcome_interaction, "ZA")) {
        "param_ZA_int"
    } else {
        "param_no_int"
    }
    paste(
        "nested_uncertainty",
        paste0("n", n),
        paste0("p", p),
        format_gamma_ZA_for_id(gamma_ZA),
        spec_id,
        sep = "_"
    )
}

run_nested_uncertainty_job <- function(job,
                                       B,
                                       B_post,
                                       posterior_seed_offset,
                                       out_dir,
                                       skip_existing = TRUE,
                                       save_csv = TRUE,
                                       worker_messages = TRUE) {
    scenario_id <- make_nested_uncertainty_scenario_id(
        n = job$n,
        p = job$p,
        gamma_ZA = job$gamma_ZA,
        outcome_interaction = job$outcome_interaction
    )
    result_file <- file.path(out_dir, paste0(scenario_id, ".rds"))
    if (skip_existing && file.exists(result_file)) {
        return(list(
            scenario_id = scenario_id,
            status = "skipped",
            error = "",
            n_success = NA_integer_,
            n_failed = NA_integer_,
            worker_pid = Sys.getpid()
        ))
    }
    if (worker_messages) {
        message("Worker ", Sys.getpid(), " starting ", scenario_id)
    }
    status <- tryCatch({
        unc_res <- run_uncertainty_simulation_nested(
            B = B,
            n = job$n,
            pars = job$pars,
            true_vals = job$true_vals,
            fit_X_given_Z = TRUE,
            outcome_interaction = job$outcome_interaction,
            B_post = B_post,
            seed = job$seed,
            posterior_seed_offset = posterior_seed_offset,
            verbose = FALSE
        )
        save_uncertainty_simulation_result_nested(
            unc_res = unc_res,
            pars = job$pars,
            scenario_id = scenario_id,
            out_dir = out_dir,
            save_csv = save_csv
        )
        list(
            scenario_id = scenario_id,
            status = "completed",
            error = "",
            n_success = unc_res$uncertainty_metadata$n_success,
            n_failed = unc_res$uncertainty_metadata$n_failed,
            worker_pid = Sys.getpid()
        )
    }, error = function(err) {
        list(
            scenario_id = scenario_id,
            status = "failed",
            error = conditionMessage(err),
            n_success = NA_integer_,
            n_failed = NA_integer_,
            worker_pid = Sys.getpid()
        )
    })
    if (worker_messages) {
        message(
            "Worker ",
            Sys.getpid(),
            " ",
            status$status,
            " ",
            scenario_id
        )
    }
    status
}

save_uncertainty_simulation_result_nested <- function(
    unc_res,
    pars,
    scenario_id,
    out_dir = "outputs/uncertainty_results_nested_ZA_interaction",
    save_csv = TRUE
) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    metadata <- unc_res$nested_metadata
    obj_to_save <- list(
        scenario_id = scenario_id,
        scenario = list(
            n = unc_res$params$n,
            p = pars$p,
            H = pars$H,
            R = pars$R,
            sparse = pars$sparse,
            gamma_ZA = pars$gamma_ZA,
            outcome_interaction = unc_res$params$outcome_interaction,
            B = unc_res$uncertainty_metadata$B,
            B_post = unc_res$uncertainty_metadata$B_post
        ),
        pars = pars,
        truth = unc_res$truth,
        intervals = unc_res$intervals,
        summary = unc_res$summary,
        convergence = unc_res$convergence,
        simulation_seeds = unc_res$simulation_seeds,
        successful_replicates = unc_res$successful_replicates,
        model_diagnostics = unc_res$model_diagnostics,
        nested_metadata = metadata,
        uncertainty_metadata = unc_res$uncertainty_metadata
    )
    rds_file <- file.path(out_dir, paste0(scenario_id, ".rds"))
    saveRDS(obj_to_save, file = rds_file)
    saved_files <- list(rds = rds_file)

    if (save_csv) {
        csv_dir <- file.path(out_dir, "csv")
        if (!dir.exists(csv_dir)) {
            dir.create(csv_dir, recursive = TRUE)
        }
        csv_files <- list()
        for (name_stub in names(unc_res$intervals)) {
            if (!is.null(unc_res$intervals[[name_stub]])) {
                path <- file.path(
                    csv_dir,
                    paste0(
                        scenario_id,
                        "_",
                        name_stub,
                        "_interval_replicates.csv"
                    )
                )
                write.csv(
                    unc_res$intervals[[name_stub]],
                    file = path,
                    row.names = FALSE
                )
                csv_files[[paste0(name_stub, "_interval_replicates")]] <- path
            }
        }
        for (name_stub in names(unc_res$summary)) {
            if (!is.null(unc_res$summary[[name_stub]])) {
                path <- file.path(
                    csv_dir,
                    paste0(
                        scenario_id,
                        "_",
                        name_stub,
                        "_interval_summary.csv"
                    )
                )
                write.csv(
                    unc_res$summary[[name_stub]],
                    file = path,
                    row.names = FALSE
                )
                csv_files[[paste0(name_stub, "_interval_summary")]] <- path
            }
        }
        diagnostics_file <- file.path(
            csv_dir,
            paste0(scenario_id, "_model_diagnostics.csv")
        )
        write.csv(
            unc_res$model_diagnostics,
            file = diagnostics_file,
            row.names = FALSE
        )
        csv_files$model_diagnostics <- diagnostics_file
        seeds_file <- file.path(
            csv_dir,
            paste0(scenario_id, "_simulation_seeds.csv")
        )
        write.csv(
            unc_res$simulation_seeds,
            file = seeds_file,
            row.names = FALSE
        )
        csv_files$simulation_seeds <- seeds_file
        truth_df <- do.call(rbind, lapply(names(unc_res$truth), function(name) {
            data.frame(
                decomposition = name,
                component = names(unc_res$truth[[name]]),
                value = as.numeric(unc_res$truth[[name]]),
                stringsAsFactors = FALSE
            )
        }))
        truth_file <- file.path(csv_dir, paste0(scenario_id, "_truth.csv"))
        write.csv(truth_df, file = truth_file, row.names = FALSE)
        csv_files$truth <- truth_file
        calibration_df <- data.frame(
            race = pars$race_labels,
            gamma_ZA = pars$gamma_ZA,
            lambda_Z_calibrated = pars$lambda_Z_calibrated,
            calibration_target_EY_by_Z = pars$calibration_target_EY_by_Z,
            calibration_achieved_EY_by_Z = pars$calibration_achieved_EY_by_Z,
            calibration_mc_n = pars$calibration_mc_n,
            calibration_seed = pars$calibration_seed,
            stringsAsFactors = FALSE
        )
        calibration_file <- file.path(
            csv_dir,
            paste0(scenario_id, "_calibration.csv")
        )
        write.csv(calibration_df, file = calibration_file, row.names = FALSE)
        csv_files$calibration <- calibration_file
        metadata_df <- data.frame(
            scenario_id = scenario_id,
            n = unc_res$params$n,
            p = pars$p,
            gamma_ZA = pars$gamma_ZA,
            outcome_interaction = unc_res$params$outcome_interaction,
            B = unc_res$uncertainty_metadata$B,
            B_post = unc_res$uncertainty_metadata$B_post,
            seed = unc_res$uncertainty_metadata$seed,
            posterior_seed_offset = unc_res$uncertainty_metadata$posterior_seed_offset,
            n_success = unc_res$uncertainty_metadata$n_success,
            n_failed = unc_res$uncertainty_metadata$n_failed,
            stringsAsFactors = FALSE
        )
        metadata_file <- file.path(
            csv_dir,
            paste0(scenario_id, "_metadata.csv")
        )
        write.csv(metadata_df, file = metadata_file, row.names = FALSE)
        csv_files$metadata <- metadata_file
        saved_files$csv <- csv_files
    }
    invisible(saved_files)
}

load_nested_uncertainty_results <- function(
    uncertainty_results = "outputs/uncertainty_results_nested_ZA_interaction"
) {
    if (is.character(uncertainty_results)) {
        if (length(uncertainty_results) == 1 &&
            dir.exists(uncertainty_results)) {
            files <- list.files(
                uncertainty_results,
                pattern = "^nested_uncertainty_.*[.]rds$",
                full.names = TRUE
            )
        } else {
            files <- uncertainty_results
        }
        if (length(files) == 0) {
            stop("No nested uncertainty result .rds files found.")
        }
        results <- lapply(files, readRDS)
        names(results) <- sub("[.]rds$", "", basename(files))
        return(results)
    }
    if (!is.list(uncertainty_results)) {
        stop(
            "uncertainty_results must be a result list, file path vector, or directory."
        )
    }
    if (!is.null(uncertainty_results$summary) &&
        !is.null(uncertainty_results$scenario)) {
        out <- list(uncertainty_results)
        names(out) <- uncertainty_results$scenario_id
        return(out)
    }
    uncertainty_results
}

make_nested_uncertainty_summary_data <- function(uncertainty_results) {
    results <- load_nested_uncertainty_results(uncertainty_results)
    rows <- unlist(lapply(results, function(obj) {
        lapply(names(obj$summary), function(estimator) {
            summary_df <- obj$summary[[estimator]]
            if (is.null(summary_df)) {
                return(NULL)
            }
            data.frame(
                scenario_id = obj$scenario_id,
                n = obj$scenario$n,
                p = obj$scenario$p,
                H = obj$scenario$H,
                gamma_ZA = obj$scenario$gamma_ZA,
                outcome_interaction = obj$scenario$outcome_interaction,
                estimator = estimator,
                summary_df,
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        })
    }), recursive = FALSE)
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) {
        stop("No nested uncertainty summaries found.")
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out <- out[
        order(
            out$gamma_ZA,
            out$p,
            out$n,
            match(out$outcome_interaction, c("none", "ZA")),
            match(
                out$estimator,
                c("topological", "modified_direct", "modified_bayes")
            ),
            out$component
        ),
        ,
        drop = FALSE
    ]
    rownames(out) <- NULL
    out
}
