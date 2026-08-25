nested_ZA_component_labels <- function() {
    c(
        delta_Z = "Delta[Z]",
        delta_X = "Delta[X]",
        delta_A = "Delta[A]",
        delta_res = "Delta[res]"
    )
}

nested_ZA_estimator_order <- function() {
    c(
        "Parametric no Z:A outcome",
        "Parametric Z:A outcome",
        "Parametric no Z:A outcome + Bayes Z|X",
        "Parametric Z:A outcome + Bayes Z|X",
        "Firth Z:A outcome",
        "Firth Z:A outcome + Bayes Z|X",
        "XGBoost"
    )
}

nested_ZA_scenario_label <- function(n, p, gamma_ZA) {
    paste0("n==", n, "*','~p==", p, "*','~gamma[ZA]==", gamma_ZA)
}

nested_ZA_estimator_label <- function(obj, estimate_name) {
    outcome_interaction <- obj$scenario$outcome_interaction
    if (identical(outcome_interaction, "xgb")) {
        return("XGBoost")
    }
    if (identical(outcome_interaction, "ZA_firth")) {
        if (identical(estimate_name, "modified_bayes")) {
            return("Firth Z:A outcome + Bayes Z|X")
        }
        return("Firth Z:A outcome")
    }
    suffix <- if (identical(outcome_interaction, "ZA")) {
        "Z:A outcome"
    } else {
        "no Z:A outcome"
    }
    if (identical(estimate_name, "modified_bayes")) {
        paste("Parametric", suffix, "+ Bayes Z|X")
    } else {
        paste("Parametric", suffix)
    }
}

nested_ZA_estimate_plot_data <- function(mat,
                                         obj,
                                         estimate_name,
                                         decomposition,
                                         components) {
    if (is.null(mat)) {
        return(NULL)
    }
    missing <- setdiff(c(components, "total"), colnames(mat))
    if (length(missing) > 0) {
        stop("Missing estimate component(s): ", paste(missing, collapse = ", "), ".")
    }
    data.frame(
        scenario_id = obj$scenario_id,
        decomposition = decomposition,
        n = obj$scenario$n,
        p = obj$scenario$p,
        gamma_ZA = obj$scenario$gamma_ZA,
        scenario = nested_ZA_scenario_label(
            obj$scenario$n,
            obj$scenario$p,
            obj$scenario$gamma_ZA
        ),
        estimator = nested_ZA_estimator_label(obj, estimate_name),
        replicate = rep(seq_len(nrow(mat)), times = length(components)),
        component = rep(components, each = nrow(mat)),
        estimate = as.numeric(mat[, components, drop = FALSE]),
        total = rep(as.numeric(mat[, "total"]), times = length(components)),
        stringsAsFactors = FALSE
    )
}

nested_ZA_truth_plot_data <- function(obj, decomposition, components) {
    truth <- obj$truth[[decomposition]]
    missing <- setdiff(c(components, "total"), names(truth))
    if (length(missing) > 0) {
        stop("Missing truth component(s): ", paste(missing, collapse = ", "), ".")
    }
    data.frame(
        scenario_id = obj$scenario_id,
        decomposition = decomposition,
        n = obj$scenario$n,
        p = obj$scenario$p,
        gamma_ZA = obj$scenario$gamma_ZA,
        scenario = nested_ZA_scenario_label(
            obj$scenario$n,
            obj$scenario$p,
            obj$scenario$gamma_ZA
        ),
        component = components,
        truth = as.numeric(truth[components]),
        total = as.numeric(truth["total"]),
        stringsAsFactors = FALSE
    )
}

read_nested_ZA_plot_data <- function(
    result_dir = "outputs/simulation_results_nested_ZA_interaction",
    gamma_ZA_values = c(0, 0.6, 1)
) {
    files <- list.files(result_dir, pattern = "[.]rds$", full.names = TRUE)
    files <- files[basename(files) != "run_metadata.rds"]
    if (length(files) == 0) {
        stop("No nested simulation result .rds files found in ", result_dir, ".")
    }

    results <- lapply(files, readRDS)
    valid <- vapply(results, function(obj) {
        is.list(obj) &&
            !is.null(obj$scenario) &&
            !is.null(obj$estimates) &&
            !is.null(obj$truth) &&
            !is.null(obj$scenario$gamma_ZA)
    }, logical(1))
    results <- results[valid]
    selected <- vapply(results, function(obj) {
        any(abs(obj$scenario$gamma_ZA - gamma_ZA_values) < 1e-12)
    }, logical(1))
    results <- results[selected]
    if (length(results) == 0) {
        stop("No matching nested simulation result objects found in ", result_dir, ".")
    }

    top_components <- c("delta_Z", "delta_X", "delta_A", "delta_res")
    modified_components <- c("delta_X", "delta_Z", "delta_A", "delta_res")
    estimate_rows <- list()
    truth_rows <- list()
    k <- 0
    t <- 0

    for (obj in results) {
        estimate_specs <- list(
            list("topological", "topological", "topological", top_components),
            list("topological_xgb", "topological_xgb", "topological", top_components),
            list("modified_direct", "modified_direct", "modified", modified_components),
            list("modified_bayes", "modified_bayes", "modified", modified_components),
            list("modified_direct_xgb", "modified_direct_xgb", "modified", modified_components)
        )
        for (spec in estimate_specs) {
            mat <- obj$estimates[[spec[[1]]]]
            if (!is.null(mat)) {
                k <- k + 1
                estimate_rows[[k]] <- nested_ZA_estimate_plot_data(
                    mat = mat,
                    obj = obj,
                    estimate_name = spec[[2]],
                    decomposition = spec[[3]],
                    components = spec[[4]]
                )
            }
        }
        t <- t + 1
        truth_rows[[t]] <- nested_ZA_truth_plot_data(
            obj,
            "topological",
            top_components
        )
        t <- t + 1
        truth_rows[[t]] <- nested_ZA_truth_plot_data(
            obj,
            "modified",
            modified_components
        )
    }

    estimates <- do.call(rbind, estimate_rows)
    truths <- do.call(rbind, truth_rows)
    truth_key <- interaction(
        truths$decomposition,
        truths$n,
        truths$p,
        truths$gamma_ZA,
        truths$component,
        drop = TRUE
    )
    truth_priority <- !grepl("_param_no_int$", truths$scenario_id)
    truth_order <- order(truth_priority)
    truths <- truths[truth_order, , drop = FALSE]
    truth_key <- truth_key[truth_order]
    truths <- truths[!duplicated(truth_key), , drop = FALSE]

    scenario_order <- unique(
        estimates[order(estimates$p, estimates$gamma_ZA, estimates$n), "scenario"]
    )
    estimator_order <- nested_ZA_estimator_order()
    component_order <- unique(c(top_components, modified_components))
    estimates$scenario <- factor(estimates$scenario, levels = scenario_order)
    truths$scenario <- factor(truths$scenario, levels = scenario_order)
    estimates$estimator <- factor(estimates$estimator, levels = estimator_order)
    estimates$component <- factor(estimates$component, levels = component_order)
    truths$component <- factor(truths$component, levels = component_order)

    list(
        estimates = estimates,
        truths = truths,
        top_components = top_components,
        modified_components = modified_components,
        estimator_order = estimator_order
    )
}

nested_ZA_plot_scale_context <- function(plot_data) {
    by_decomposition <- lapply(c("topological", "modified"), function(name) {
        list(
            estimates = plot_data$estimates[
                plot_data$estimates$decomposition == name,
                ,
                drop = FALSE
            ],
            truths = plot_data$truths[
                plot_data$truths$decomposition == name,
                ,
                drop = FALSE
            ]
        )
    })
    names(by_decomposition) <- c("topological", "modified")
    component_absolute <- lapply(by_decomposition, function(dat) {
        simulation_y_limits_by_component(
            dat,
            show_mean_overlay = TRUE,
            expand = 0.08
        )
    })
    percent_data <- lapply(by_decomposition, function(dat) {
        scale_simulation_plot_data(dat, y_scale = "percent_total")
    })
    percent_global <- lapply(percent_data, function(dat) {
        simulation_y_limits(dat, show_mean_overlay = TRUE, expand = 0.08)
    })
    component_percent <- lapply(percent_data, function(dat) {
        simulation_y_limits_by_component(
            dat,
            show_mean_overlay = TRUE,
            expand = 0.08
        )
    })
    list(
        component_absolute = component_absolute,
        percent_global = percent_global,
        component_percent = component_percent
    )
}

make_nested_ZA_boxplot <- function(plot_data,
                                   decomposition,
                                   gamma_ZA,
                                   y_scale_mode = c("fixed", "independent_y"),
                                   component_y_limits = NULL,
                                   y_value_scale = c("absolute", "percent_total"),
                                   global_y_limits = NULL,
                                   component_labels = nested_ZA_component_labels()) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The ggplot2 package is required to make nested simulation plots.")
    }
    y_scale_mode <- match.arg(y_scale_mode)
    y_value_scale <- match.arg(y_value_scale)
    if (identical(y_scale_mode, "independent_y") &&
        !requireNamespace("ggh4x", quietly = TRUE)) {
        stop("The ggh4x package is required for independent y-axis scales.")
    }

    est <- plot_data$estimates[
        plot_data$estimates$decomposition == decomposition &
            abs(plot_data$estimates$gamma_ZA - gamma_ZA) < 1e-12,
        ,
        drop = FALSE
    ]
    tru <- plot_data$truths[
        plot_data$truths$decomposition == decomposition &
            abs(plot_data$truths$gamma_ZA - gamma_ZA) < 1e-12,
        ,
        drop = FALSE
    ]
    if (nrow(est) == 0 || nrow(tru) == 0) {
        stop("No plotting data found for gamma_ZA = ", gamma_ZA, ".")
    }

    scaled <- scale_simulation_plot_data(
        list(estimates = est, truths = tru),
        y_scale = y_value_scale
    )
    est <- scaled$estimates
    tru <- scaled$truths
    mean_data <- summarize_simulation_mean_overlay(est)
    tru$truth_label <- if (identical(y_value_scale, "percent_total")) {
        paste0("truth = ", formatC(tru$truth, format = "f", digits = 1), "%")
    } else {
        paste0("truth = ", formatC(tru$truth, format = "f", digits = 4))
    }
    tru$label_x <- factor(
        plot_data$estimator_order[1],
        levels = plot_data$estimator_order
    )
    if (identical(y_scale_mode, "fixed")) {
        if (identical(y_value_scale, "percent_total")) {
            if (is.null(global_y_limits)) {
                stop("global_y_limits must be supplied for fixed percent-total scales.")
            }
            tru$label_y <- global_y_limits[2] - 0.02 * diff(global_y_limits)
        } else {
            tru$label_y <- 0.22
        }
    } else {
        if (is.null(component_y_limits)) {
            stop("component_y_limits must be supplied for independent y-axis scales.")
        }
        tru$label_y <- vapply(as.character(tru$component), function(component) {
            limits <- component_y_limits[[component]]
            limits[2] - 0.02 * diff(limits)
        }, numeric(1))
    }

    plot <- ggplot2::ggplot() +
        ggplot2::geom_hline(
            data = tru,
            ggplot2::aes(yintercept = truth),
            linetype = "dashed",
            color = "gray20",
            linewidth = 0.3
        ) +
        ggplot2::geom_boxplot(
            data = est,
            ggplot2::aes(x = estimator, y = estimate, fill = estimator),
            width = 0.65,
            outlier.size = 0.4,
            linewidth = 0.25
        ) +
        ggplot2::geom_errorbar(
            data = mean_data,
            ggplot2::aes(
                x = estimator,
                ymin = lower,
                ymax = upper,
                color = estimator
            ),
            width = 0.18,
            linewidth = 0.35,
            position = ggplot2::position_nudge(x = 0.46),
            show.legend = FALSE
        ) +
        ggplot2::geom_point(
            data = mean_data,
            ggplot2::aes(
                x = estimator,
                y = mean_estimate,
                color = estimator,
                fill = estimator
            ),
            shape = 21,
            size = 0.8,
            stroke = 0.35,
            position = ggplot2::position_nudge(x = 0.46),
            show.legend = FALSE
        ) +
        ggplot2::geom_text(
            data = tru,
            ggplot2::aes(x = label_x, y = label_y, label = truth_label),
            hjust = 0,
            vjust = 1,
            size = 2.2,
            color = "gray20",
            inherit.aes = FALSE
        )

    facet_labeller <- ggplot2::labeller(
        scenario = ggplot2::label_parsed,
        component = ggplot2::as_labeller(
            component_labels,
            default = ggplot2::label_parsed
        )
    )
    if (identical(y_scale_mode, "fixed")) {
        plot <- plot +
            ggplot2::facet_grid(
                rows = ggplot2::vars(scenario),
                cols = ggplot2::vars(component),
                scales = "fixed",
                labeller = facet_labeller
            )
        if (identical(y_value_scale, "percent_total")) {
            plot <- plot +
                ggplot2::scale_y_continuous(
                    expand = ggplot2::expansion(mult = c(0, 0))
                ) +
                ggplot2::coord_cartesian(ylim = global_y_limits, clip = "off")
        } else {
            plot <- plot +
                ggplot2::scale_y_continuous(
                    breaks = seq(0, 0.2, by = 0.05),
                    expand = ggplot2::expansion(mult = c(0, 0))
                ) +
                ggplot2::coord_cartesian(ylim = c(0, 0.225), clip = "off")
        }
    } else {
        y_scales <- lapply(names(component_y_limits), function(component) {
            limits <- component_y_limits[[component]]
            stats::as.formula(paste0(
                "component == '", component,
                "' ~ ggplot2::scale_y_continuous(limits = c(",
                format(limits[1], digits = 17, scientific = FALSE), ", ",
                format(limits[2], digits = 17, scientific = FALSE),
                "), expand = ggplot2::expansion(mult = c(0, 0)))"
            ))
        })
        plot <- plot +
            ggh4x::facet_grid2(
                rows = ggplot2::vars(scenario),
                cols = ggplot2::vars(component),
                scales = "free_y",
                independent = "y",
                labeller = facet_labeller
            ) +
            ggh4x::facetted_pos_scales(y = y_scales)
    }

    plot +
        ggplot2::labs(
            x = NULL,
            y = if (identical(y_value_scale, "percent_total")) {
                "Variance component (% of total)"
            } else {
                "Estimated variance component"
            },
            fill = NULL
        ) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "gray92",
                color = "gray60",
                linewidth = 0.3
            ),
            strip.text = ggplot2::element_text(face = "bold", size = 8),
            axis.text.x = ggplot2::element_text(
                angle = 35,
                hjust = 1,
                vjust = 1,
                size = 7
            ),
            axis.text.y = ggplot2::element_text(size = 7),
            axis.title.y = ggplot2::element_text(size = 8),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size = 8),
            legend.key.size = grid::unit(0.35, "cm")
        )
}

format_nested_gamma_ZA_for_file <- function(gamma_ZA) {
    sub(
        "[.]0+$",
        "",
        sub("[.]", "p", format(gamma_ZA, scientific = FALSE, trim = TRUE))
    )
}

save_nested_ZA_interaction_plots <- function(
    result_dir = "outputs/simulation_results_nested_ZA_interaction",
    out_dir = "outputs/figures/nested_ZA_interaction",
    gamma_ZA_values = c(0, 0.6, 1),
    write_fixed_absolute = TRUE,
    write_independent_absolute = TRUE,
    write_fixed_percent = TRUE,
    write_independent_percent = TRUE,
    save_plot_data = TRUE,
    width = 10,
    height = 10
) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    plot_data <- read_nested_ZA_plot_data(
        result_dir = result_dir,
        gamma_ZA_values = gamma_ZA_values
    )
    scales <- nested_ZA_plot_scale_context(plot_data)
    saved <- character()

    save_plot <- function(plot, file_name) {
        path <- file.path(out_dir, file_name)
        ggplot2::ggsave(
            filename = path,
            plot = plot,
            width = width,
            height = height
        )
        saved <<- c(saved, path)
    }

    for (gamma_ZA in gamma_ZA_values) {
        gamma_file <- format_nested_gamma_ZA_for_file(gamma_ZA)
        for (decomposition in c("topological", "modified")) {
            prefix <- paste0(decomposition, "_nested_gammaZA", gamma_file)
            if (write_fixed_absolute) {
                save_plot(
                    make_nested_ZA_boxplot(
                        plot_data,
                        decomposition,
                        gamma_ZA
                    ),
                    paste0(prefix, "_boxplots.pdf")
                )
            }
            if (write_independent_absolute) {
                save_plot(
                    make_nested_ZA_boxplot(
                        plot_data,
                        decomposition,
                        gamma_ZA,
                        y_scale_mode = "independent_y",
                        component_y_limits = scales$component_absolute[[decomposition]]
                    ),
                    paste0(prefix, "_boxplots_independent_y.pdf")
                )
            }
            if (write_fixed_percent) {
                save_plot(
                    make_nested_ZA_boxplot(
                        plot_data,
                        decomposition,
                        gamma_ZA,
                        y_value_scale = "percent_total",
                        global_y_limits = scales$percent_global[[decomposition]]
                    ),
                    paste0(prefix, "_boxplots_percent_total.pdf")
                )
            }
            if (write_independent_percent) {
                save_plot(
                    make_nested_ZA_boxplot(
                        plot_data,
                        decomposition,
                        gamma_ZA,
                        y_scale_mode = "independent_y",
                        component_y_limits = scales$component_percent[[decomposition]],
                        y_value_scale = "percent_total"
                    ),
                    paste0(prefix, "_boxplots_independent_y_percent_total.pdf")
                )
            }
        }

        if (save_plot_data) {
            plot_data_path <- file.path(
                out_dir,
                paste0("nested_gammaZA", gamma_file, "_plot_data.rds")
            )
            saveRDS(
                list(
                    estimates = plot_data$estimates[
                        abs(plot_data$estimates$gamma_ZA - gamma_ZA) < 1e-12,
                        ,
                        drop = FALSE
                    ],
                    truths = plot_data$truths[
                        abs(plot_data$truths$gamma_ZA - gamma_ZA) < 1e-12,
                        ,
                        drop = FALSE
                    ]
                ),
                file = plot_data_path
            )
            saved <- c(saved, plot_data_path)
        }
    }
    invisible(saved)
}
