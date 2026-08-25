default_simulation_plot_scenarios <- function() {
    expand.grid(
        n = c(500, 1500),
        p = c(5, 15),
        KEEP.OUT.ATTRS = FALSE
    )[c(1, 3, 2, 4), ]
}

simulation_component_labels <- function() {
    c(
        delta_Z = "Delta_Z",
        delta_X = "Delta_X",
        delta_A = "Delta_A",
        delta_res = "Delta_res"
    )
}

simulation_estimator_colors <- function() {
    c(
        Parametric = "#3B6EA8",
        Bayes = "#4D8B31",
        XGBoost = "#B24A3B"
    )
}

read_simulation_result_files <- function(result_dir = "outputs/simulation_results",
                                         pattern = "[.]rds$") {
    files <- list.files(result_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
        stop("No .rds simulation result files found in ", result_dir, ".")
    }
    out <- lapply(files, readRDS)
    is_result <- vapply(out, function(obj) {
        is.list(obj) &&
            length(obj$scenario_id) == 1 &&
            !is.null(obj$scenario) &&
            !is.null(obj$estimates) &&
            !is.null(obj$truth)
    }, logical(1))
    out <- out[is_result]
    files <- files[is_result]
    if (length(out) == 0) {
        stop("No simulation result .rds files found in ", result_dir, ".")
    }
    names(out) <- vapply(out, function(obj) obj$scenario_id, character(1))
    out
}

filter_simulation_results <- function(results, scenarios = default_simulation_plot_scenarios()) {
    keep <- vapply(results, function(obj) {
        any(obj$scenario$n == scenarios$n & obj$scenario$p == scenarios$p)
    }, logical(1))
    out <- results[keep]
    found <- do.call(rbind, lapply(out, function(obj) {
        data.frame(n = obj$scenario$n, p = obj$scenario$p)
    }))
    missing <- scenarios[!vapply(seq_len(nrow(scenarios)), function(i) {
        any(found$n == scenarios$n[i] & found$p == scenarios$p[i])
    }, logical(1)), , drop = FALSE]
    if (nrow(missing) > 0) {
        msg <- paste(sprintf("n = %s, p = %s", missing$n, missing$p), collapse = "; ")
        stop("Missing simulation result files for scenario(s): ", msg, ".")
    }
    out
}

scenario_label <- function(n, p) {
    paste0("n = ", n, ", p = ", p)
}

scenario_levels <- function(scenarios) {
    scenario_label(scenarios$n, scenarios$p)
}

matrix_to_estimate_data <- function(mat, scenario, estimator, components) {
    missing <- setdiff(components, colnames(mat))
    if (length(missing) > 0) {
        stop("Missing estimate component(s): ", paste(missing, collapse = ", "), ".")
    }
    if (!"total" %in% colnames(mat)) {
        stop("Missing estimate component: total.")
    }
    total <- as.numeric(mat[, "total"])
    mat <- mat[, components, drop = FALSE]
    data.frame(
        scenario = scenario,
        estimator = estimator,
        replicate = rep(seq_len(nrow(mat)), times = ncol(mat)),
        component = rep(colnames(mat), each = nrow(mat)),
        estimate = as.numeric(mat),
        total = rep(total, times = ncol(mat)),
        stringsAsFactors = FALSE
    )
}

truth_to_data <- function(truth, scenario, components) {
    missing <- setdiff(components, names(truth))
    if (length(missing) > 0) {
        stop("Missing truth component(s): ", paste(missing, collapse = ", "), ".")
    }
    if (!"total" %in% names(truth)) {
        stop("Missing truth component: total.")
    }
    data.frame(
        scenario = scenario,
        component = components,
        truth = as.numeric(truth[components]),
        total = as.numeric(truth["total"]),
        stringsAsFactors = FALSE
    )
}

tidy_topological_simulation_data <- function(results,
                                             scenarios = default_simulation_plot_scenarios(),
                                             components = c("delta_Z", "delta_X", "delta_A", "delta_res")) {
    results <- filter_simulation_results(results, scenarios = scenarios)
    scenario_order <- scenario_levels(scenarios)
    est_list <- list()
    truth_list <- list()
    k <- 0
    for (obj in results) {
        scen <- scenario_label(obj$scenario$n, obj$scenario$p)
        k <- k + 1
        est_list[[k]] <- matrix_to_estimate_data(
            obj$estimates$topological,
            scenario = scen,
            estimator = "Parametric",
            components = components
        )
        if (!is.null(obj$estimates$topological_xgb)) {
            k <- k + 1
            est_list[[k]] <- matrix_to_estimate_data(
                obj$estimates$topological_xgb,
                scenario = scen,
                estimator = "XGBoost",
                components = components
            )
        }
        truth_list[[length(truth_list) + 1]] <- truth_to_data(
            obj$truth$topological,
            scenario = scen,
            components = components
        )
    }
    estimates <- do.call(rbind, est_list)
    truths <- do.call(rbind, truth_list)
    estimates$scenario <- factor(estimates$scenario, levels = scenario_order)
    truths$scenario <- factor(truths$scenario, levels = scenario_order)
    estimates$component <- factor(estimates$component, levels = components)
    truths$component <- factor(truths$component, levels = components)
    estimates$estimator <- factor(estimates$estimator, levels = c("Parametric", "Bayes", "XGBoost"))
    list(estimates = estimates, truths = truths)
}

tidy_modified_simulation_data <- function(results,
                                          scenarios = default_simulation_plot_scenarios(),
                                          components = c("delta_X", "delta_Z", "delta_A", "delta_res")) {
    results <- filter_simulation_results(results, scenarios = scenarios)
    scenario_order <- scenario_levels(scenarios)
    est_list <- list()
    truth_list <- list()
    k <- 0
    for (obj in results) {
        scen <- scenario_label(obj$scenario$n, obj$scenario$p)
        k <- k + 1
        est_list[[k]] <- matrix_to_estimate_data(
            obj$estimates$modified_direct,
            scenario = scen,
            estimator = "Parametric",
            components = components
        )
        if (!is.null(obj$estimates$modified_bayes)) {
            k <- k + 1
            est_list[[k]] <- matrix_to_estimate_data(
                obj$estimates$modified_bayes,
                scenario = scen,
                estimator = "Bayes",
                components = components
            )
        }
        if (!is.null(obj$estimates$modified_direct_xgb)) {
            k <- k + 1
            est_list[[k]] <- matrix_to_estimate_data(
                obj$estimates$modified_direct_xgb,
                scenario = scen,
                estimator = "XGBoost",
                components = components
            )
        }
        truth_list[[length(truth_list) + 1]] <- truth_to_data(
            obj$truth$modified,
            scenario = scen,
            components = components
        )
    }
    estimates <- do.call(rbind, est_list)
    truths <- do.call(rbind, truth_list)
    estimates$scenario <- factor(estimates$scenario, levels = scenario_order)
    truths$scenario <- factor(truths$scenario, levels = scenario_order)
    estimates$component <- factor(estimates$component, levels = components)
    truths$component <- factor(truths$component, levels = components)
    estimates$estimator <- factor(estimates$estimator, levels = c("Parametric", "Bayes", "XGBoost"))
    list(estimates = estimates, truths = truths)
}

summarize_simulation_mean_overlay <- function(estimates) {
    groups <- split(
        estimates,
        list(estimates$scenario, estimates$component, estimates$estimator),
        drop = TRUE
    )
    out <- lapply(groups, function(dat) {
        b <- sum(!is.na(dat$estimate))
        mean_estimate <- mean(dat$estimate, na.rm = TRUE)
        mcse <- if (b > 1) stats::sd(dat$estimate, na.rm = TRUE) / sqrt(b) else NA_real_
        data.frame(
            scenario = dat$scenario[1],
            component = dat$component[1],
            estimator = dat$estimator[1],
            mean_estimate = mean_estimate,
            mcse = mcse,
            lower = mean_estimate - 1.96 * mcse,
            upper = mean_estimate + 1.96 * mcse,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    out$scenario <- factor(out$scenario, levels = levels(estimates$scenario))
    out$component <- factor(out$component, levels = levels(estimates$component))
    out$estimator <- factor(out$estimator, levels = levels(estimates$estimator))
    out
}

scale_simulation_plot_data <- function(plot_data,
                                       y_scale = c("absolute", "percent_total")) {
    y_scale <- match.arg(y_scale)
    if (identical(y_scale, "absolute")) {
        return(plot_data)
    }
    if (!"total" %in% names(plot_data$estimates) || !"total" %in% names(plot_data$truths)) {
        stop("Percent-total plots require total columns in estimates and truths.")
    }
    out <- plot_data
    if (any(!is.finite(out$estimates$total) | out$estimates$total == 0)) {
        stop("Cannot scale estimates by total variance: non-finite or zero total encountered.")
    }
    if (any(!is.finite(out$truths$total) | out$truths$total == 0)) {
        stop("Cannot scale truths by total variance: non-finite or zero total encountered.")
    }
    out$estimates$estimate <- 100 * out$estimates$estimate / out$estimates$total
    out$truths$truth <- 100 * out$truths$truth / out$truths$total
    out
}

simulation_y_limits <- function(plot_data,
                                show_mean_overlay = FALSE,
                                expand = 0.05) {
    vals <- c(plot_data$estimates$estimate, plot_data$truths$truth)
    if (show_mean_overlay) {
        mean_overlay_data <- summarize_simulation_mean_overlay(plot_data$estimates)
        vals <- c(vals, mean_overlay_data$lower, mean_overlay_data$upper)
    }
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) {
        return(NULL)
    }
    lim <- range(vals)
    span <- diff(lim)
    if (span <= 0) {
        span <- max(abs(lim), 1) * 0.1
    }
    lim + c(-1, 1) * expand * span
}

simulation_y_limits_by_component <- function(plot_data,
                                             show_mean_overlay = FALSE,
                                             expand = 0.05) {
    component_levels <- levels(plot_data$estimates$component)
    out <- lapply(component_levels, function(component) {
        comp_data <- list(
            estimates = plot_data$estimates[plot_data$estimates$component == component, , drop = FALSE],
            truths = plot_data$truths[plot_data$truths$component == component, , drop = FALSE]
        )
        simulation_y_limits(comp_data, show_mean_overlay = show_mean_overlay, expand = expand)
    })
    names(out) <- component_levels
    out
}

make_simulation_boxplot <- function(plot_data,
                                    title = NULL,
                                    component_labels = simulation_component_labels(),
                                    estimator_colors = simulation_estimator_colors(),
                                    facet_mode = c("default", "independent_y"),
                                    show_mean_overlay = FALSE,
                                    mean_overlay_side = c("right", "left"),
                                    mean_overlay_width = 0.18,
                                    mean_overlay_nudge = 0.46,
                                    mean_overlay_size = 1.2,
                                    mean_overlay_linewidth = 0.35,
                                    y_limits = NULL,
                                    y_scale = c("absolute", "percent_total")) {
    facet_mode <- match.arg(facet_mode)
    mean_overlay_side <- match.arg(mean_overlay_side)
    y_scale <- match.arg(y_scale)
    plot_data <- scale_simulation_plot_data(plot_data, y_scale = y_scale)
    y_label <- if (identical(y_scale, "percent_total")) {
        "Variance component (% of total)"
    } else {
        "Estimated variance component"
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The ggplot2 package is required to make simulation plots.")
    }
    if (identical(facet_mode, "independent_y") &&
        !requireNamespace("ggh4x", quietly = TRUE)) {
        stop("The ggh4x package is required for facet_mode = 'independent_y'.")
    }
    facet <- if (identical(facet_mode, "independent_y")) {
        ggh4x::facet_grid2(
            rows = ggplot2::vars(scenario),
            cols = ggplot2::vars(component),
            scales = "free_y",
            independent = "y",
            labeller = ggplot2::labeller(component = ggplot2::as_labeller(component_labels))
        )
    } else {
        ggplot2::facet_grid(
            rows = ggplot2::vars(scenario),
            cols = ggplot2::vars(component),
            scales = "free_y",
            labeller = ggplot2::labeller(component = ggplot2::as_labeller(component_labels))
        )
    }
    plot <- ggplot2::ggplot() +
        ggplot2::geom_hline(
            data = plot_data$truths,
            ggplot2::aes(yintercept = truth),
            linetype = "dashed",
            color = "gray20",
            linewidth = 0.35
        ) +
        ggplot2::geom_boxplot(
            data = plot_data$estimates,
            ggplot2::aes(x = estimator, y = estimate, fill = estimator),
            width = 0.62,
            outlier.size = 0.5,
            linewidth = 0.25
        ) +
        facet +
        ggplot2::scale_fill_manual(values = estimator_colors, drop = TRUE) +
        ggplot2::labs(title = title, x = NULL, y = y_label, fill = NULL) +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", hjust = 0),
            strip.background = ggplot2::element_rect(fill = "gray92", color = "gray60", linewidth = 0.3),
            strip.text = ggplot2::element_text(face = "bold", size = 8),
            axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1, size = 7),
            axis.text.y = ggplot2::element_text(size = 7),
            axis.title.y = ggplot2::element_text(size = 8),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size = 8),
            legend.key.size = grid::unit(0.35, "cm"),
            plot.margin = ggplot2::margin(5, 5, 5, 5)
        )
    if (show_mean_overlay) {
        mean_overlay_data <- summarize_simulation_mean_overlay(plot_data$estimates)
        mean_overlay_nudge <- if (identical(mean_overlay_side, "right")) {
            abs(mean_overlay_nudge)
        } else {
            -abs(mean_overlay_nudge)
        }
        plot <- plot +
            ggplot2::geom_errorbar(
                data = mean_overlay_data,
                ggplot2::aes(x = estimator, ymin = lower, ymax = upper, color = estimator),
                width = mean_overlay_width,
                linewidth = mean_overlay_linewidth,
                position = ggplot2::position_nudge(x = mean_overlay_nudge),
                show.legend = FALSE
            ) +
            ggplot2::geom_point(
                data = mean_overlay_data,
                ggplot2::aes(x = estimator, y = mean_estimate, color = estimator, fill = estimator),
                shape = 21,
                size = mean_overlay_size,
                stroke = mean_overlay_linewidth,
                position = ggplot2::position_nudge(x = mean_overlay_nudge),
                show.legend = FALSE
            ) +
            ggplot2::scale_color_manual(values = estimator_colors, drop = TRUE)
    }
    if (!is.null(y_limits)) {
        if (is.list(y_limits) && !is.null(names(y_limits))) {
            if (!requireNamespace("ggh4x", quietly = TRUE)) {
                stop("The ggh4x package is required for component-specific y limits.")
            }
            y_scales <- lapply(names(y_limits), function(component) {
                stats::as.formula(
                    paste0("component == '", component, "' ~ ggplot2::scale_y_continuous(limits = c(", y_limits[[component]][1], ", ", y_limits[[component]][2], "))")
                )
            })
            plot <- plot + ggh4x::facetted_pos_scales(y = y_scales)
        } else {
            plot <- plot + ggplot2::coord_cartesian(ylim = y_limits)
        }
    }
    plot
}

build_simulation_plots <- function(result_dir = "outputs/simulation_results",
                                   scenarios = default_simulation_plot_scenarios(),
                                   facet_mode = c("default", "independent_y"),
                                   show_mean_overlay = FALSE,
                                   mean_overlay_side = c("right", "left"),
                                   mean_overlay_width = 0.18,
                                   mean_overlay_nudge = 0.46,
                                   mean_overlay_size = 1.2,
                                   mean_overlay_linewidth = 0.35,
                                   y_limits = NULL,
                                   y_scale = c("absolute", "percent_total")) {
    facet_mode <- match.arg(facet_mode)
    mean_overlay_side <- match.arg(mean_overlay_side)
    y_scale <- match.arg(y_scale)
    results <- read_simulation_result_files(result_dir = result_dir)
    topological_data <- tidy_topological_simulation_data(results, scenarios = scenarios)
    modified_data <- tidy_modified_simulation_data(results, scenarios = scenarios)
    topological_limit_data <- scale_simulation_plot_data(topological_data, y_scale = y_scale)
    modified_limit_data <- scale_simulation_plot_data(modified_data, y_scale = y_scale)
    if (identical(y_limits, "global")) {
        y_limits <- list(
            topological = simulation_y_limits(topological_limit_data, show_mean_overlay = show_mean_overlay),
            modified = simulation_y_limits(modified_limit_data, show_mean_overlay = show_mean_overlay)
        )
    } else if (identical(y_limits, "component")) {
        y_limits <- list(
            topological = simulation_y_limits_by_component(topological_limit_data, show_mean_overlay = show_mean_overlay),
            modified = simulation_y_limits_by_component(modified_limit_data, show_mean_overlay = show_mean_overlay)
        )
    }
    topological_y_limits <- if (is.list(y_limits)) y_limits$topological else y_limits
    modified_y_limits <- if (is.list(y_limits)) y_limits$modified else y_limits
    topological_plot <- make_simulation_boxplot(
        topological_data,
        title = "Topological decomposition",
        facet_mode = facet_mode,
        show_mean_overlay = show_mean_overlay,
        mean_overlay_side = mean_overlay_side,
        mean_overlay_width = mean_overlay_width,
        mean_overlay_nudge = mean_overlay_nudge,
        mean_overlay_size = mean_overlay_size,
        mean_overlay_linewidth = mean_overlay_linewidth,
        y_limits = topological_y_limits,
        y_scale = y_scale
    )
    modified_plot <- make_simulation_boxplot(
        modified_data,
        title = "Modified-order decomposition",
        facet_mode = facet_mode,
        show_mean_overlay = show_mean_overlay,
        mean_overlay_side = mean_overlay_side,
        mean_overlay_width = mean_overlay_width,
        mean_overlay_nudge = mean_overlay_nudge,
        mean_overlay_size = mean_overlay_size,
        mean_overlay_linewidth = mean_overlay_linewidth,
        y_limits = modified_y_limits,
        y_scale = y_scale
    )
    list(
        topological_plot = topological_plot,
        modified_plot = modified_plot,
        topological_data = topological_data,
        modified_data = modified_data
    )
}

save_simulation_plots <- function(plot_objects,
                                  out_dir = "outputs/figures",
                                  topological_file = "topological_simulation_boxplots.pdf",
                                  modified_file = "modified_simulation_boxplots.pdf",
                                  width = 8.2,
                                  height = 8.2,
                                  save_tidy_data = TRUE) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The ggplot2 package is required to save simulation plots.")
    }
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    topological_path <- file.path(out_dir, topological_file)
    modified_path <- file.path(out_dir, modified_file)
    ggplot2::ggsave(topological_path, plot_objects$topological_plot, width = width, height = height)
    ggplot2::ggsave(modified_path, plot_objects$modified_plot, width = width, height = height)
    saved <- list(topological = topological_path, modified = modified_path)
    if (save_tidy_data) {
        tidy_path <- file.path(out_dir, "simulation_boxplot_tidy_data.rds")
        saveRDS(
            list(
                topological = plot_objects$topological_data,
                modified = plot_objects$modified_data
            ),
            tidy_path
        )
        saved$tidy_data <- tidy_path
    }
    invisible(saved)
}

save_simulation_plots_by_n <- function(result_dir = "outputs/simulation_results",
                                       out_dir = "outputs/figures",
                                       n_values = c(500, 1500),
                                       p_values = c(5, 15),
                                       facet_mode = c("default", "independent_y"),
                                       show_mean_overlay = TRUE,
                                       y_limits = "global",
                                       y_scale = c("absolute", "percent_total"),
                                       file_suffix = NULL,
                                       width = 8.2,
                                       height = 4.8) {
    facet_mode <- match.arg(facet_mode)
    y_scale <- match.arg(y_scale)
    suffix <- if (is.null(file_suffix)) {
        paste0(
            if (identical(facet_mode, "independent_y")) "_independent_y" else "",
            if (identical(y_scale, "percent_total")) "_percent_total" else ""
        )
    } else {
        file_suffix
    }
    shared_y_limits <- y_limits
    if (identical(y_limits, "global")) {
        results <- read_simulation_result_files(result_dir = result_dir)
        all_scenarios <- expand.grid(
            n = n_values,
            p = p_values,
            KEEP.OUT.ATTRS = FALSE
        )
        all_scenarios <- all_scenarios[order(match(all_scenarios$n, n_values), match(all_scenarios$p, p_values)), ]
        topological_data <- tidy_topological_simulation_data(results, scenarios = all_scenarios)
        modified_data <- tidy_modified_simulation_data(results, scenarios = all_scenarios)
        topological_data <- scale_simulation_plot_data(topological_data, y_scale = y_scale)
        modified_data <- scale_simulation_plot_data(modified_data, y_scale = y_scale)
        shared_y_limits <- list(
            topological = simulation_y_limits(topological_data, show_mean_overlay = show_mean_overlay),
            modified = simulation_y_limits(modified_data, show_mean_overlay = show_mean_overlay)
        )
    } else if (identical(y_limits, "component")) {
        results <- read_simulation_result_files(result_dir = result_dir)
        all_scenarios <- expand.grid(
            n = n_values,
            p = p_values,
            KEEP.OUT.ATTRS = FALSE
        )
        all_scenarios <- all_scenarios[order(match(all_scenarios$n, n_values), match(all_scenarios$p, p_values)), ]
        topological_data <- tidy_topological_simulation_data(results, scenarios = all_scenarios)
        modified_data <- tidy_modified_simulation_data(results, scenarios = all_scenarios)
        topological_data <- scale_simulation_plot_data(topological_data, y_scale = y_scale)
        modified_data <- scale_simulation_plot_data(modified_data, y_scale = y_scale)
        shared_y_limits <- list(
            topological = simulation_y_limits_by_component(topological_data, show_mean_overlay = show_mean_overlay),
            modified = simulation_y_limits_by_component(modified_data, show_mean_overlay = show_mean_overlay)
        )
    }
    saved <- list()
    for (n_value in n_values) {
        scenarios <- data.frame(
            n = rep(n_value, length(p_values)),
            p = p_values
        )
        plot_objects <- build_simulation_plots(
            result_dir = result_dir,
            scenarios = scenarios,
            facet_mode = facet_mode,
            show_mean_overlay = show_mean_overlay,
            y_limits = shared_y_limits,
            y_scale = y_scale
        )
        saved[[paste0("n", n_value)]] <- save_simulation_plots(
            plot_objects,
            out_dir = out_dir,
            topological_file = paste0("topological_simulation_boxplots_n", n_value, suffix, ".pdf"),
            modified_file = paste0("modified_simulation_boxplots_n", n_value, suffix, ".pdf"),
            width = width,
            height = height,
            save_tidy_data = FALSE
        )
    }
    invisible(saved)
}
