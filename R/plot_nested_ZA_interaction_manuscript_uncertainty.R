manuscript_uncertainty_estimator_levels <- function() {
    c(
        "Topological: main effects",
        "Topological: Z:A",
        "Modified: main effects",
        "Modified: Z:A"
    )
}

prepare_nested_ZA_manuscript_coverage_data <- function(
    summary_file = paste0(
        "outputs/uncertainty_results_nested_ZA_interaction/",
        "nested_uncertainty_combined_summary.csv"
    ),
    gamma_ZA_values = c(0, 0.6, 1)
) {
    coverage <- utils::read.csv(summary_file, stringsAsFactors = FALSE)
    required <- c(
        "n", "p", "gamma_ZA", "outcome_interaction", "estimator",
        "component", "coverage", "n_success"
    )
    missing <- setdiff(required, names(coverage))
    if (length(missing) > 0) {
        stop("Missing uncertainty summary column(s): ", paste(missing, collapse = ", "))
    }

    selected_gamma <- vapply(coverage$gamma_ZA, function(value) {
        any(abs(value - gamma_ZA_values) < 1e-12)
    }, logical(1))
    coverage <- coverage[
        selected_gamma &
            coverage$outcome_interaction %in% c("none", "ZA") &
            coverage$estimator %in% c("topological", "modified_direct") &
            coverage$component %in% c("delta_Z", "delta_X", "delta_A", "delta_res"),
        ,
        drop = FALSE
    ]
    if (nrow(coverage) == 0) {
        stop("No matching uncertainty coverage results were found in ", summary_file, ".")
    }

    coverage$order <- factor(
        ifelse(coverage$estimator == "topological", "Topological", "Modified"),
        levels = c("Topological", "Modified")
    )
    coverage$model <- factor(
        ifelse(coverage$outcome_interaction == "none", "main effects", "Z:A"),
        levels = c("main effects", "Z:A")
    )
    coverage$plot_estimator <- factor(
        paste0(coverage$order, ": ", coverage$model),
        levels = manuscript_uncertainty_estimator_levels()
    )
    coverage$scenario <- factor(
        paste0("n==", coverage$n, "*','~p==", coverage$p),
        levels = unique(paste0(
            "n==", c(500, 1500, 500, 1500),
            "*','~p==", c(5, 5, 15, 15)
        ))
    )
    coverage$component <- factor(
        coverage$component,
        levels = c("delta_Z", "delta_X", "delta_A", "delta_res")
    )
    coverage
}

make_nested_ZA_manuscript_coverage_plot <- function(coverage_data,
                                                    gamma_ZA,
                                                    y_limits = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The ggplot2 package is required to make manuscript figures.")
    }
    dat <- coverage_data[
        abs(coverage_data$gamma_ZA - gamma_ZA) < 1e-12,
        ,
        drop = FALSE
    ]
    if (nrow(dat) == 0) {
        stop("No coverage data found for gamma_ZA = ", gamma_ZA, ".")
    }

    plot <- ggplot2::ggplot(
        dat,
        ggplot2::aes(
            x = plot_estimator,
            y = coverage,
            color = order,
            shape = model
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0.95,
            color = "gray35",
            linetype = "dashed",
            linewidth = 0.4
        ) +
        ggplot2::geom_vline(
            xintercept = 2.5,
            color = "gray75",
            linetype = "dotted",
            linewidth = 0.3
        ) +
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
            values = c("main effects" = 16, "Z:A" = 17)
        ) +
        ggplot2::labs(
            x = NULL,
            y = "Coverage probability",
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

save_nested_ZA_manuscript_coverage_plots <- function(
    summary_file = paste0(
        "outputs/uncertainty_results_nested_ZA_interaction/",
        "nested_uncertainty_combined_summary.csv"
    ),
    out_dir = "outputs/figures/nested_ZA_interaction/manuscript/uncertainty",
    gamma_ZA_values = c(0, 0.6, 1),
    width = 8,
    height = 8
) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    coverage_data <- prepare_nested_ZA_manuscript_coverage_data(
        summary_file = summary_file,
        gamma_ZA_values = gamma_ZA_values
    )
    coverage_min <- min(coverage_data$coverage, na.rm = TRUE)
    common_y_limits <- c(max(0, coverage_min - 0.04), 1)
    saved <- character()
    for (gamma_ZA in gamma_ZA_values) {
        path <- file.path(
            out_dir,
            paste0(
                "nested_gammaZA",
                format_nested_gamma_ZA_for_file(gamma_ZA),
                "_coverage.pdf"
            )
        )
        ggplot2::ggsave(
            filename = path,
            plot = make_nested_ZA_manuscript_coverage_plot(
                coverage_data,
                gamma_ZA,
                y_limits = common_y_limits
            ),
            width = width,
            height = height
        )
        saved <- c(saved, path)
    }
    invisible(saved)
}
