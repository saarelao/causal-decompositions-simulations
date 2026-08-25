manuscript_nested_estimator_levels <- function() {
    c(
        "Topological: main effects",
        "Topological: Z:A",
        "Topological: XGBoost",
        "Modified: main effects",
        "Modified: Z:A",
        "Modified: XGBoost"
    )
}

manuscript_nested_estimator_label <- function(decomposition, estimator) {
    model_labels <- c(
        "Parametric no Z:A outcome" = "main effects",
        "Parametric Z:A outcome" = "Z:A",
        "XGBoost" = "XGBoost"
    )
    model <- unname(model_labels[as.character(estimator)])
    order_label <- ifelse(decomposition == "topological", "Topological", "Modified")
    paste0(order_label, ": ", model)
}

prepare_nested_ZA_manuscript_bias_data <- function(
    result_dir = "outputs/simulation_results_nested_ZA_interaction",
    gamma_ZA_values = c(0, 0.6, 1)
) {
    plot_data <- read_nested_ZA_plot_data(
        result_dir = result_dir,
        gamma_ZA_values = gamma_ZA_values
    )
    keep_estimators <- c(
        "Parametric no Z:A outcome",
        "Parametric Z:A outcome",
        "XGBoost"
    )
    estimates <- plot_data$estimates[
        as.character(plot_data$estimates$estimator) %in% keep_estimators,
        ,
        drop = FALSE
    ]

    truth_columns <- c(
        "decomposition", "n", "p", "gamma_ZA", "component", "truth"
    )
    truths <- plot_data$truths[, truth_columns, drop = FALSE]
    estimates <- merge(
        estimates,
        truths,
        by = c("decomposition", "n", "p", "gamma_ZA", "component"),
        all.x = TRUE,
        sort = FALSE
    )
    if (anyNA(estimates$truth)) {
        stop("Truth values could not be matched to all selected estimates.")
    }

    estimates$bias <- estimates$estimate - estimates$truth
    estimates$manuscript_estimator <- manuscript_nested_estimator_label(
        estimates$decomposition,
        estimates$estimator
    )
    estimator_levels <- manuscript_nested_estimator_levels()
    estimates$manuscript_estimator <- factor(
        estimates$manuscript_estimator,
        levels = estimator_levels
    )
    estimates$model <- factor(
        sub("^[^:]+: ", "", estimates$manuscript_estimator),
        levels = c("main effects", "Z:A", "XGBoost")
    )
    estimates$order <- factor(
        ifelse(estimates$decomposition == "topological", "Topological", "Modified"),
        levels = c("Topological", "Modified")
    )
    estimates$scenario <- factor(
        paste0("n==", estimates$n, "*','~p==", estimates$p),
        levels = unique(paste0(
            "n==", c(500, 1500, 500, 1500),
            "*','~p==", c(5, 5, 15, 15)
        ))
    )

    group_columns <- c(
        "n", "p", "gamma_ZA", "scenario", "component",
        "manuscript_estimator", "model", "order"
    )
    groups <- interaction(estimates[, group_columns], drop = TRUE)
    split_estimates <- split(estimates, groups)
    summaries <- lapply(split_estimates, function(dat) {
        finite_bias <- dat$bias[is.finite(dat$bias)]
        if (length(finite_bias) == 0) {
            return(NULL)
        }
        data.frame(
            n = dat$n[1],
            p = dat$p[1],
            gamma_ZA = dat$gamma_ZA[1],
            scenario = dat$scenario[1],
            component = dat$component[1],
            estimator = dat$manuscript_estimator[1],
            model = dat$model[1],
            order = dat$order[1],
            mean_bias = mean(finite_bias),
            lower = unname(stats::quantile(finite_bias, 0.025, names = FALSE)),
            upper = unname(stats::quantile(finite_bias, 0.975, names = FALSE)),
            mcse_bias = stats::sd(finite_bias) / sqrt(length(finite_bias)),
            n_success = length(finite_bias),
            stringsAsFactors = FALSE
        )
    })
    summaries <- do.call(rbind, summaries)
    rownames(summaries) <- NULL
    summaries$scenario <- factor(
        as.character(summaries$scenario),
        levels = levels(estimates$scenario)
    )
    summaries$component <- factor(
        as.character(summaries$component),
        levels = c("delta_Z", "delta_X", "delta_A", "delta_res")
    )
    summaries$estimator <- factor(
        as.character(summaries$estimator),
        levels = estimator_levels
    )
    summaries$model <- factor(
        as.character(summaries$model),
        levels = c("main effects", "Z:A", "XGBoost")
    )
    summaries$order <- factor(
        as.character(summaries$order),
        levels = c("Topological", "Modified")
    )
    summaries$mc_lower <- summaries$mean_bias - 1.96 * summaries$mcse_bias
    summaries$mc_upper <- summaries$mean_bias + 1.96 * summaries$mcse_bias

    list(estimates = estimates, summaries = summaries)
}

make_nested_ZA_manuscript_bias_plot <- function(bias_data,
                                                gamma_ZA,
                                                y_limits = NULL,
                                                show_mc_intervals = FALSE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The ggplot2 package is required to make manuscript figures.")
    }
    dat <- bias_data$summaries[
        abs(bias_data$summaries$gamma_ZA - gamma_ZA) < 1e-12,
        ,
        drop = FALSE
    ]
    if (nrow(dat) == 0) {
        stop("No manuscript plotting data found for gamma_ZA = ", gamma_ZA, ".")
    }

    plot <- ggplot2::ggplot(
        dat,
        ggplot2::aes(
            x = estimator,
            y = mean_bias,
            ymin = lower,
            ymax = upper,
            color = order,
            shape = model
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0,
            color = "gray35",
            linetype = "dashed",
            linewidth = 0.35
        ) +
        ggplot2::geom_vline(
            xintercept = 3.5,
            color = "gray75",
            linetype = "dotted",
            linewidth = 0.3
        )

    if (show_mc_intervals) {
        plot <- plot + ggplot2::geom_errorbar(
            ggplot2::aes(ymin = mc_lower, ymax = mc_upper),
            width = 0.34,
            linewidth = 0.25
        )
    }

    plot <- plot +
        ggplot2::geom_errorbar(width = 0.18, linewidth = 0.45) +
        ggplot2::geom_point(size = 2.1, stroke = 0.65) +
        ggplot2::facet_grid(
            rows = ggplot2::vars(scenario),
            cols = ggplot2::vars(component),
            scales = "fixed",
            labeller = ggplot2::labeller(
                scenario = ggplot2::label_parsed,
                component = ggplot2::as_labeller(
                    nested_ZA_component_labels(),
                    default = ggplot2::label_parsed
                )
            )
        ) +
        ggplot2::scale_color_manual(
            values = c(Topological = "black", Modified = "gray55")
        ) +
        ggplot2::scale_shape_manual(
            values = c("main effects" = 16, "Z:A" = 17, XGBoost = 15)
        ) +
        ggplot2::labs(
            x = NULL,
            y = expression("Estimation error " * (widehat(Delta) - Delta)),
            color = "Decomposition",
            shape = "Outcome model"
        ) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
            strip.background = ggplot2::element_rect(
                fill = "gray92",
                color = "gray60",
                linewidth = 0.3
            ),
            strip.text = ggplot2::element_text(face = "bold", size = 8),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = 7),
            axis.title.y = ggplot2::element_text(size = 8),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.key.width = grid::unit(0.55, "cm")
        )

    if (!is.null(y_limits)) {
        plot <- plot + ggplot2::coord_cartesian(ylim = y_limits)
    }
    plot
}

save_nested_ZA_manuscript_bias_plots <- function(
    result_dir = "outputs/simulation_results_nested_ZA_interaction",
    out_dir = "outputs/figures/nested_ZA_interaction/manuscript",
    gamma_ZA_values = c(0, 0.6, 1),
    show_mc_intervals = FALSE,
    width = 8,
    height = 8
) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    bias_data <- prepare_nested_ZA_manuscript_bias_data(
        result_dir = result_dir,
        gamma_ZA_values = gamma_ZA_values
    )
    selected <- bias_data$summaries[
        vapply(bias_data$summaries$gamma_ZA, function(value) {
            any(abs(value - gamma_ZA_values) < 1e-12)
        }, logical(1)),
        ,
        drop = FALSE
    ]
    max_abs_error <- max(abs(c(selected$lower, selected$upper)), na.rm = TRUE)
    common_y_limits <- c(-1, 1) * max_abs_error * 1.04
    saved <- character()
    for (gamma_ZA in gamma_ZA_values) {
        path <- file.path(
            out_dir,
            paste0(
                "nested_gammaZA",
                format_nested_gamma_ZA_for_file(gamma_ZA),
                "_estimation_error.pdf"
            )
        )
        ggplot2::ggsave(
            filename = path,
            plot = make_nested_ZA_manuscript_bias_plot(
                bias_data,
                gamma_ZA,
                y_limits = common_y_limits,
                show_mc_intervals = show_mc_intervals
            ),
            width = width,
            height = height
        )
        saved <- c(saved, path)
    }
    invisible(saved)
}
