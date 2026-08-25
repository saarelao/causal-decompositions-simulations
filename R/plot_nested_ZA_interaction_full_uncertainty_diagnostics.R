full_nested_uncertainty_estimator_levels <- function(ordering) {
    ordering <- match.arg(ordering, c("topological", "modified"))
    if (ordering == "topological") {
        c(
            "Parametric no Z:A outcome",
            "Parametric Z:A outcome"
        )
    } else {
        c(
            "Parametric no Z:A outcome",
            "Parametric Z:A outcome",
            "Parametric no Z:A outcome + Bayes Z|X",
            "Parametric Z:A outcome + Bayes Z|X"
        )
    }
}

full_nested_uncertainty_estimator_palette <- function(ordering) {
    ordering <- match.arg(ordering, c("topological", "modified"))
    if (!requireNamespace("scales", quietly = TRUE)) {
        stop("The scales package is required to reproduce the point-figure colors.")
    }
    labels <- full_nested_uncertainty_estimator_levels(ordering)
    colors <- if (ordering == "topological") {
        scales::hue_pal()(4)[seq_along(labels)]
    } else {
        scales::hue_pal()(7)[seq_along(labels)]
    }
    stats::setNames(colors, labels)
}

prepare_nested_ZA_full_uncertainty_diagnostics <- function(
    summary_file = paste0(
        "outputs/uncertainty_results_nested_ZA_interaction/",
        "nested_uncertainty_combined_summary.csv"
    ),
    gamma_ZA_values = c(0, 0.6, 1)
) {
    dat <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
    required <- c(
        "n", "p", "gamma_ZA", "outcome_interaction", "estimator",
        "component", "coverage", "mean_posterior_sd", "empirical_sd_point"
    )
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0) {
        stop("Missing uncertainty summary column(s): ", paste(missing, collapse = ", "))
    }
    selected_gamma <- vapply(dat$gamma_ZA, function(value) {
        any(abs(value - gamma_ZA_values) < 1e-12)
    }, logical(1))
    dat <- dat[
        selected_gamma &
            dat$outcome_interaction %in% c("none", "ZA") &
            dat$estimator %in% c("topological", "modified_direct", "modified_bayes") &
            dat$component %in% c("delta_Z", "delta_X", "delta_A", "delta_res"),
        ,
        drop = FALSE
    ]
    if (nrow(dat) != 288L) {
        stop("Expected 288 uncertainty summary rows but found ", nrow(dat), ".")
    }
    dat$order <- ifelse(dat$estimator == "topological", "topological", "modified")
    model_label <- ifelse(
        dat$outcome_interaction == "none",
        "no Z:A outcome",
        "Z:A outcome"
    )
    bayes_suffix <- ifelse(dat$estimator == "modified_bayes", " + Bayes Z|X", "")
    dat$estimator_label <- paste0("Parametric ", model_label, bayes_suffix)
    dat$sd_ratio <- dat$empirical_sd_point / dat$mean_posterior_sd
    if (any(!is.finite(dat$sd_ratio))) {
        stop("Non-finite empirical-to-posterior SD ratio found.")
    }

    coverage_rows <- dat
    coverage_rows$diagnostic <- "Coverage"
    coverage_rows$value <- coverage_rows$coverage
    ratio_rows <- dat
    ratio_rows$diagnostic <- "Empirical SD / mean posterior SD"
    ratio_rows$value <- ratio_rows$sd_ratio
    plot_data <- rbind(coverage_rows, ratio_rows)
    plot_data$diagnostic <- factor(
        plot_data$diagnostic,
        levels = c("Coverage", "Empirical SD / mean posterior SD")
    )
    plot_data$scenario <- factor(
        paste0(
            "n==", plot_data$n,
            "*','~p==", plot_data$p,
            "*','~gamma[ZA]==", plot_data$gamma_ZA
        ),
        levels = unique(paste0(
            "n==", c(500, 1500, 500, 1500),
            "*','~p==", c(5, 5, 15, 15),
            "*','~gamma[ZA]==", rep(gamma_ZA_values[1], 4)
        ))
    )
    plot_data$component <- factor(
        plot_data$component,
        levels = c("delta_Z", "delta_X", "delta_A", "delta_res")
    )
    plot_data
}

make_nested_ZA_full_uncertainty_diagnostic_plot <- function(plot_data,
                                                            ordering,
                                                            gamma_ZA,
                                                            y_limits) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The ggplot2 package is required to make uncertainty figures.")
    }
    ordering <- match.arg(ordering, c("topological", "modified"))
    dat <- plot_data[
        plot_data$order == ordering &
            abs(plot_data$gamma_ZA - gamma_ZA) < 1e-12,
        ,
        drop = FALSE
    ]
    if (nrow(dat) == 0) {
        stop("No uncertainty diagnostics found for the requested figure.")
    }
    estimator_levels <- full_nested_uncertainty_estimator_levels(ordering)
    dat$estimator_label <- factor(dat$estimator_label, levels = estimator_levels)
    dat$scenario <- factor(
        paste0("n==", dat$n, "*','~p==", dat$p),
        levels = unique(paste0(
            "n==", c(500, 1500, 500, 1500),
            "*','~p==", c(5, 5, 15, 15)
        ))
    )

    ggplot2::ggplot(
        dat,
        ggplot2::aes(
            x = estimator_label,
            y = value,
            fill = estimator_label,
            shape = diagnostic,
            group = interaction(estimator_label, diagnostic)
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0.95,
            color = "gray55",
            linetype = "dashed",
            linewidth = 0.35
        ) +
        ggplot2::geom_hline(
            yintercept = 1,
            color = "gray20",
            linetype = "dashed",
            linewidth = 0.35
        ) +
        ggplot2::geom_point(
            size = 2.1,
            stroke = 0.65,
            color = "black",
            position = ggplot2::position_dodge(width = 0.34)
        ) +
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
        ggplot2::scale_fill_manual(
            values = full_nested_uncertainty_estimator_palette(ordering),
            drop = FALSE
        ) +
        ggplot2::scale_shape_manual(
            values = c(
                Coverage = 21,
                "Empirical SD / mean posterior SD" = 24
            )
        ) +
        ggplot2::coord_cartesian(ylim = y_limits) +
        ggplot2::labs(
            x = NULL,
            y = "Coverage or SD ratio",
            fill = "Estimator",
            shape = "Diagnostic"
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(
                override.aes = list(shape = 21, color = "black")
            ),
            shape = ggplot2::guide_legend(
                override.aes = list(fill = "white", color = "black")
            )
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
            legend.box = "vertical",
            legend.text = ggplot2::element_text(size = 7),
            legend.key.width = grid::unit(0.45, "cm")
        )
}

save_nested_ZA_full_uncertainty_diagnostic_plots <- function(
    summary_file = paste0(
        "outputs/uncertainty_results_nested_ZA_interaction/",
        "nested_uncertainty_combined_summary.csv"
    ),
    out_dir = "outputs/figures/nested_ZA_interaction/full_uncertainty_diagnostics",
    gamma_ZA_values = c(0, 0.6, 1),
    width = 10,
    height = 10
) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    plot_data <- prepare_nested_ZA_full_uncertainty_diagnostics(
        summary_file,
        gamma_ZA_values
    )
    range_values <- range(c(plot_data$value, 0.95, 1), na.rm = TRUE)
    padding <- 0.04 * diff(range_values)
    common_y_limits <- c(
        max(0, range_values[1] - padding),
        range_values[2] + padding
    )

    saved <- character()
    for (gamma_ZA in gamma_ZA_values) {
        gamma_file <- format_nested_gamma_ZA_for_file(gamma_ZA)
        for (ordering in c("topological", "modified")) {
            path <- file.path(
                out_dir,
                paste0(
                    ordering,
                    "_nested_gammaZA",
                    gamma_file,
                    "_uncertainty_diagnostics.pdf"
                )
            )
            ggplot2::ggsave(
                filename = path,
                plot = make_nested_ZA_full_uncertainty_diagnostic_plot(
                    plot_data,
                    ordering,
                    gamma_ZA,
                    y_limits = common_y_limits
                ),
                width = width,
                height = height
            )
            saved <- c(saved, path)
        }
    }
    invisible(saved)
}
