log_nested_scenario_error <- function(error_log, scenario_id, err) {
    cat(
        paste0(
            "[", Sys.time(), "] ", scenario_id, "\n",
            conditionMessage(err), "\n\n"
        ),
        file = error_log,
        append = TRUE
    )
}

make_nested_scenario_pars <- function(p_scenario,
                                      gamma_ZA,
                                      calibration_mc_n = 200000,
                                      calibration_seed = 202606) {
    make_sim_params_nested(
        R = 3,
        H = 5,
        p = p_scenario,
        sparse = TRUE,
        n_active_A = ceiling(p_scenario / 2),
        n_active_Y = ceiling(p_scenario / 2),
        s_ZX = 1.1,
        s_ZA = 1,
        s_XA = 1.1,
        s_AY = 1,
        s_ZY = 0.8,
        s_XY = 1,
        intercept_Y = -0.3,
        rho_X = 0.25,
        gamma_ZA = gamma_ZA,
        calibration_mc_n = calibration_mc_n,
        calibration_seed = calibration_seed
    )
}
