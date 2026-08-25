rm(list = ls())

source("R/plot_simulation_results.R")
source("R/plot_nested_ZA_interaction_results.R")
source("R/plot_nested_ZA_interaction_manuscript_results.R")
source("R/plot_nested_ZA_interaction_manuscript_uncertainty.R")
source("R/plot_nested_ZA_interaction_full_uncertainty_diagnostics.R")

point_dir <- Sys.getenv(
    "POINT_RESULTS_DIR",
    unset = "outputs/simulation_results_nested_ZA_interaction"
)
uncertainty_dir <- Sys.getenv(
    "UNCERTAINTY_RESULTS_DIR",
    unset = "outputs/uncertainty_results_nested_ZA_interaction"
)
figure_dir <- Sys.getenv(
    "FIGURE_DIR",
    unset = "outputs/figures/nested_ZA_interaction"
)
gamma_values <- c(0, 0.6, 1)
summary_file <- file.path(
    uncertainty_dir,
    "nested_uncertainty_combined_summary.csv"
)

save_nested_ZA_interaction_plots(
    result_dir = point_dir,
    out_dir = figure_dir,
    gamma_ZA_values = gamma_values,
    write_fixed_absolute = FALSE,
    write_independent_absolute = FALSE,
    write_fixed_percent = TRUE,
    write_independent_percent = FALSE,
    save_plot_data = FALSE,
    width = 10,
    height = 10
)

save_nested_ZA_manuscript_bias_plots(
    result_dir = point_dir,
    out_dir = file.path(figure_dir, "manuscript"),
    gamma_ZA_values = gamma_values,
    show_mc_intervals = FALSE,
    width = 8,
    height = 8
)

save_nested_ZA_manuscript_coverage_plots(
    summary_file = summary_file,
    out_dir = file.path(figure_dir, "manuscript", "uncertainty"),
    gamma_ZA_values = gamma_values,
    width = 8,
    height = 8
)

save_nested_ZA_full_uncertainty_diagnostic_plots(
    summary_file = summary_file,
    out_dir = file.path(figure_dir, "full_uncertainty_diagnostics"),
    gamma_ZA_values = gamma_values,
    width = 10,
    height = 10
)

message("Figure generation complete: ", figure_dir)
