softmax_mat <- function(eta) {
    eta_centered <- eta - apply(eta, 1, max)
    exp_eta <- exp(eta_centered)
    exp_eta/rowSums(exp_eta)
}

sample_multinom_rows <- function(prob_mat, labels = NULL) {
    n <- nrow(prob_mat)
    H <- ncol(prob_mat)
    idx <- apply(prob_mat, 1, function(p) sample.int(H, size = 1, prob = p))
    if (is.null(labels)) {
        idx
    }
    else {
        factor(labels[idx], levels = labels)
    }
}

default_race_labels <- function(R) {
    if (R == 3) {
        c("neutral", "advantaged", "disadvantaged")
    }
    else {
        paste0("race", seq_len(R))
    }
}

default_race_direct_scores <- function(R) {
    if (R == 3) {
        c(0, 1, -1)
    }
    else {
        seq(-1, 1, length.out = R)
    }
}

default_race_severity_scores <- function(R) {
    if (R == 3) {
        c(0, -1, 1)
    }
    else {
        -seq(-1, 1, length.out = R)
    }
}

default_hospital_labels <- function(H) {
    if (H == 5) {
        c("H1_low", "H2_low", "H3_mid", "H4_high", "H5_high")
    }
    else {
        paste0("H", seq_len(H))
    }
}

default_hospital_scores <- function(H) {
    if (H == 5) {
        c(-1, -0.5, 0, 0.5, 1)
    }
    else {
        seq(-1, 1, length.out = H)
    }
}

make_positive_pattern <- function(p, sparse = TRUE, n_active = NULL) {
    if (is.null(n_active)) {
        n_active <- if (sparse) 
            max(2, ceiling(p/2))
        else p
    }
    n_active <- min(n_active, p)
    out <- rep(0, p)
    vals <- seq(1, 0.5, length.out = n_active)
    out[seq_len(n_active)] <- vals
    s <- sqrt(sum(out^2))
    if (s > 0) 
        out <- out/s
    out
}

make_sim_params <- function(R = 3, H = 5, p = 5, race_probs = NULL, race_direct_scores = NULL, race_severity_scores = NULL, hospital_scores = NULL, rho_X = 0.25, sparse = TRUE, n_active_A = NULL, n_active_Y = NULL, s_ZX = 0.9, s_ZA = 0.7, s_XA = 1, s_AY = 0.9, s_ZY = 0.5, s_XY = 0.9, intercept_Y = -0.3, alpha_A = NULL, seed = 123) {
    set.seed(seed)
    if (is.null(race_probs)) {
        race_probs <- rep(1/R, R)
    }
    if (is.null(race_direct_scores)) {
        race_direct_scores <- default_race_direct_scores(R)
    }
    if (is.null(race_severity_scores)) {
        race_severity_scores <- default_race_severity_scores(R)
    }
    if (is.null(hospital_scores)) {
        hospital_scores <- default_hospital_scores(H)
    }
    if (is.null(alpha_A)) {
        alpha_A <- rep(0, H)
    }
    stopifnot(length(race_probs) == R)
    stopifnot(length(race_direct_scores) == R)
    stopifnot(length(race_severity_scores) == R)
    stopifnot(length(hospital_scores) == H)
    stopifnot(length(alpha_A) == H)
    race_labels <- default_race_labels(R)
    hospital_labels <- default_hospital_labels(H)
    mu_pattern <- seq(1, 0.6, length.out = p)
    mu_pattern <- mu_pattern/sqrt(sum(mu_pattern^2))
    Sigma_X <- matrix(rho_X, nrow = p, ncol = p)
    diag(Sigma_X) <- 1
    beta_A_base <- make_positive_pattern(p, sparse = sparse, n_active = n_active_A)
    beta_Y_base <- make_positive_pattern(p, sparse = sparse, n_active = n_active_Y)
    better_hospitals <- hospital_labels[hospital_scores > 0]
    worse_hospitals <- hospital_labels[hospital_scores < 0]
    middle_hospitals <- hospital_labels[which.min(abs(hospital_scores))]
    baseline_hospital <- hospital_labels[which.min(abs(hospital_scores))]
    baseline_race <- race_labels[which.min(abs(race_direct_scores))]
    list(R = R, H = H, p = p, race_probs = race_probs, race_labels = race_labels, race_direct_scores = race_direct_scores, race_severity_scores = race_severity_scores, hospital_labels = hospital_labels, hospital_scores = hospital_scores, better_hospitals = better_hospitals, worse_hospitals = worse_hospitals, middle_hospitals = middle_hospitals, baseline_hospital = baseline_hospital, baseline_race = baseline_race, rho_X = rho_X, Sigma_X = Sigma_X, mu_pattern = mu_pattern, beta_A_base = beta_A_base, 
        beta_Y_base = beta_Y_base, sparse = sparse, s_ZX = s_ZX, s_ZA = s_ZA, s_XA = s_XA, s_AY = s_AY, s_ZY = s_ZY, s_XY = s_XY, intercept_Y = intercept_Y, alpha_A = alpha_A, seed = seed)
}

simulate_one <- function(n, pars) {
    R <- pars$R
    H <- pars$H
    p <- pars$p
    Z_idx <- sample.int(R, size = n, replace = TRUE, prob = pars$race_probs)
    Z <- factor(pars$race_labels[Z_idx], levels = pars$race_labels)
    z_direct <- pars$race_direct_scores[Z_idx]
    z_severity <- pars$race_severity_scores[Z_idx]
    X <- matrix(NA_real_, nrow = n, ncol = p)
    for (r in seq_len(R)) {
        idx <- which(Z_idx == r)
        if (length(idx) > 0) {
            mu_r <- pars$s_ZX * pars$race_severity_scores[r] * pars$mu_pattern
            X[idx, ] <- MASS::mvrnorm(n = length(idx), mu = mu_r, Sigma = pars$Sigma_X)
        }
    }
    colnames(X) <- paste0("X", seq_len(p))
    x_score_A <- as.vector(X %*% pars$beta_A_base)
    x_score_Y <- as.vector(X %*% pars$beta_Y_base)
    eta_A <- matrix(0, nrow = n, ncol = H)
    for (h in seq_len(H)) {
        qh <- pars$hospital_scores[h]
        eta_A[, h] <- pars$alpha_A[h] + pars$s_ZA * z_direct * qh + pars$s_XA * x_score_A * qh
    }
    prob_A <- softmax_mat(eta_A)
    A <- sample_multinom_rows(prob_A, labels = pars$hospital_labels)
    h_score <- pars$hospital_scores[match(as.character(A), pars$hospital_labels)]
    lp_Y <- pars$intercept_Y + pars$s_AY * h_score + pars$s_ZY * z_direct + pars$s_XY * x_score_Y
    p_Y <- plogis(lp_Y)
    Y <- rbinom(n, size = 1, prob = p_Y)
    dat <- data.frame(id = seq_len(n), Z = Z, A = A, Y = Y, X, stringsAsFactors = FALSE)
    attr(dat, "true_prob_A") <- prob_A
    attr(dat, "true_prob_Y") <- p_Y
    attr(dat, "true_x_score_A") <- x_score_A
    attr(dat, "true_x_score_Y") <- x_score_Y
    dat
}

fit_oracle_models <- function(dat, pars) {
    p <- pars$p
    x_names <- paste0("X", seq_len(p))
    dat$Z <- relevel(dat$Z, ref = pars$baseline_race)
    dat$A <- relevel(dat$A, ref = pars$baseline_hospital)
    mu_hat <- sapply(split(dat[, x_names, drop = FALSE], dat$Z), colMeans)
    mu_hat <- t(mu_hat)
    Sigma_num <- matrix(0, nrow = p, ncol = p)
    denom <- 0
    for (g in levels(dat$Z)) {
        Xg <- as.matrix(dat[dat$Z == g, x_names, drop = FALSE])
        if (nrow(Xg) > 1) {
            cg <- colMeans(Xg)
            centered <- sweep(Xg, 2, cg, "-")
            Sigma_num <- Sigma_num + crossprod(centered)
            denom <- denom + (nrow(Xg) - 1)
        }
    }
    Sigma_hat <- Sigma_num/denom
    colnames(Sigma_hat) <- rownames(Sigma_hat) <- x_names
    form_A <- as.formula(paste("A ~ Z +", paste(x_names, collapse = " + ")))
    fit_A <- tryCatch(nnet::multinom(form_A, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    form_Y <- as.formula(paste("Y ~ A + Z +", paste(x_names, collapse = " + ")))
    fit_Y <- tryCatch(glm(form_Y, data = dat, family = binomial()), error = function(e) e)
    list(fit_X = list(mu_hat = mu_hat, Sigma_hat = Sigma_hat), fit_A = fit_A, fit_Y = fit_Y)
}

replicate_diagnostics <- function(dat, pars, fits) {
    better <- dat$A %in% pars$better_hospitals
    pY_by_Z <- tapply(dat$Y, dat$Z, mean)
    pBetter_by_Z <- tapply(better, dat$Z, mean)
    meanX_by_Z <- sapply(split(dat[, paste0("X", seq_len(pars$p)), drop = FALSE], dat$Z), function(dd) mean(as.matrix(dd)))
    pY_by_A <- tapply(dat$Y, dat$A, mean)
    pA_by_Z <- prop.table(table(dat$Z, dat$A), 1)
    ok_A <- inherits(fits$fit_A, "multinom")
    ok_Y <- inherits(fits$fit_Y, "glm")
    list(pY_overall = mean(dat$Y), pY_by_Z = pY_by_Z, pBetter_by_Z = pBetter_by_Z, meanX_by_Z = meanX_by_Z, pY_by_A = pY_by_A, pA_by_Z = pA_by_Z, ok_A = ok_A, ok_Y = ok_Y)
}

run_simulation <- function(B = 50, n = 3000, pars = make_sim_params(), keep_data = FALSE, keep_fits = FALSE, verbose = TRUE) {
    out <- vector("list", B)
    for (b in seq_len(B)) {
        dat_b <- simulate_one(n = n, pars = pars)
        fits_b <- fit_oracle_models(dat_b, pars)
        diag_b <- replicate_diagnostics(dat_b, pars, fits_b)
        out[[b]] <- list(diagnostics = diag_b, data = if (keep_data) dat_b else NULL, fits = if (keep_fits) fits_b else NULL)
        if (verbose && (b%%max(1, floor(B/10)) == 0)) {
            message("Completed replicate ", b, " / ", B)
        }
    }
    race_levels <- pars$race_labels
    hosp_levels <- pars$hospital_labels
    pY_overall <- sapply(out, function(obj) obj$diagnostics$pY_overall)
    pY_by_Z_mat <- do.call(rbind, lapply(out, function(obj) obj$diagnostics$pY_by_Z[race_levels]))
    pBetter_by_Z_mat <- do.call(rbind, lapply(out, function(obj) obj$diagnostics$pBetter_by_Z[race_levels]))
    meanX_by_Z_mat <- do.call(rbind, lapply(out, function(obj) obj$diagnostics$meanX_by_Z[race_levels]))
    pY_by_A_mat <- do.call(rbind, lapply(out, function(obj) obj$diagnostics$pY_by_A[hosp_levels]))
    ok_A <- sapply(out, function(obj) obj$diagnostics$ok_A)
    ok_Y <- sapply(out, function(obj) obj$diagnostics$ok_Y)
    summary_list <- list(pY_overall_mean = mean(pY_overall), pY_overall_sd = sd(pY_overall), pY_by_Z_mean = colMeans(pY_by_Z_mat, na.rm = TRUE), pBetter_by_Z_mean = colMeans(pBetter_by_Z_mat, na.rm = TRUE), meanX_by_Z_mean = colMeans(meanX_by_Z_mat, na.rm = TRUE), pY_by_A_mean = colMeans(pY_by_A_mat, na.rm = TRUE), conv_rate_A = mean(ok_A), conv_rate_Y = mean(ok_Y))
    list(params = pars, reps = out, summary = summary_list)
}

print_sim_summary <- function(sim) {
    cat("\n========================================\n")
    cat("Simulation summary\n")
    cat("========================================\n")
    cat(sprintf("Mean Pr(Y=1): %.3f\n", sim$summary$pY_overall_mean))
    cat(sprintf("SD   Pr(Y=1): %.3f\n", sim$summary$pY_overall_sd))
    cat(sprintf("Convergence rate A model: %.3f\n", sim$summary$conv_rate_A))
    cat(sprintf("Convergence rate Y model: %.3f\n", sim$summary$conv_rate_Y))
    cat("\nMean disease severity (mean of X) by race:\n")
    print(round(sim$summary$meanX_by_Z_mean, 3))
    cat("\nMean Pr(A in better hospitals | Z):\n")
    print(round(sim$summary$pBetter_by_Z_mean, 3))
    cat("\nMean Pr(Y=1 | Z):\n")
    print(round(sim$summary$pY_by_Z_mean, 3))
    cat("\nMean Pr(Y=1 | A):\n")
    print(round(sim$summary$pY_by_A_mean, 3))
    cat("========================================\n")
}

mvnorm_logdens_rows <- function(X, mu, Sigma) {
    p <- ncol(X)
    cholS <- chol(Sigma)
    logdetS <- 2 * sum(log(diag(cholS)))
    quad <- mahalanobis(X, center = mu, cov = Sigma, inverted = FALSE)
    -0.5 * (p * log(2 * pi) + logdetS + quad)
}

true_mu_X_given_Z <- function(pars) {
    R <- pars$R
    out <- vector("list", R)
    for (r in seq_len(R)) {
        out[[r]] <- pars$s_ZX * pars$race_severity_scores[r] * pars$mu_pattern
    }
    out
}

true_prob_A_given_XZ <- function(X, Z_idx, pars) {
    n <- nrow(X)
    H <- pars$H
    z_direct <- pars$race_direct_scores[Z_idx]
    x_score_A <- as.vector(X %*% pars$beta_A_base)
    eta_A <- matrix(0, nrow = n, ncol = H)
    for (h in seq_len(H)) {
        qh <- pars$hospital_scores[h]
        eta_A[, h] <- pars$alpha_A[h] + pars$s_ZA * z_direct * qh + pars$s_XA * x_score_A * qh
    }
    softmax_mat(eta_A)
}

true_mean_Y_given_XZA_allA <- function(X, Z_idx, pars) {
    n <- nrow(X)
    H <- pars$H
    z_direct <- pars$race_direct_scores[Z_idx]
    x_score_Y <- as.vector(X %*% pars$beta_Y_base)
    m_mat <- matrix(0, nrow = n, ncol = H)
    for (h in seq_len(H)) {
        qh <- pars$hospital_scores[h]
        lp <- pars$intercept_Y + pars$s_AY * qh + pars$s_ZY * z_direct + pars$s_XY * x_score_Y
        m_mat[, h] <- plogis(lp)
    }
    m_mat
}

true_mean_Y_given_XZA_allA_nested <- function(X, Z_idx, pars) {
    n <- nrow(X)
    H <- pars$H
    z_direct <- pars$race_direct_scores[Z_idx]
    x_score_Y <- as.vector(X %*% pars$beta_Y_base)
    gamma_ZA <- if (!is.null(pars$gamma_ZA)) pars$gamma_ZA else 0
    lambda_Z <- pars$lambda_Z_calibrated
    if (is.null(lambda_Z)) {
        lambda_Z <- pars$s_ZY * pars$race_direct_scores
    }
    m_mat <- matrix(0, nrow = n, ncol = H)
    for (h in seq_len(H)) {
        qh <- pars$hospital_scores[h]
        lp <- pars$intercept_Y +
            pars$s_XY * x_score_Y +
            pars$s_AY * qh +
            lambda_Z[Z_idx] +
            gamma_ZA * z_direct * qh
        m_mat[, h] <- plogis(lp)
    }
    m_mat
}

simulate_one_nested <- function(n, pars) {
    R <- pars$R
    H <- pars$H
    p <- pars$p
    Z_idx <- sample.int(R, size = n, replace = TRUE, prob = pars$race_probs)
    Z <- factor(pars$race_labels[Z_idx], levels = pars$race_labels)
    z_direct <- pars$race_direct_scores[Z_idx]
    z_severity <- pars$race_severity_scores[Z_idx]
    X <- matrix(NA_real_, nrow = n, ncol = p)
    for (r in seq_len(R)) {
        idx <- which(Z_idx == r)
        if (length(idx) > 0) {
            mu_r <- pars$s_ZX * pars$race_severity_scores[r] * pars$mu_pattern
            X[idx, ] <- MASS::mvrnorm(n = length(idx), mu = mu_r, Sigma = pars$Sigma_X)
        }
    }
    colnames(X) <- paste0("X", seq_len(p))
    x_score_A <- as.vector(X %*% pars$beta_A_base)
    x_score_Y <- as.vector(X %*% pars$beta_Y_base)
    eta_A <- matrix(0, nrow = n, ncol = H)
    for (h in seq_len(H)) {
        qh <- pars$hospital_scores[h]
        eta_A[, h] <- pars$alpha_A[h] + pars$s_ZA * z_direct * qh + pars$s_XA * x_score_A * qh
    }
    prob_A <- softmax_mat(eta_A)
    A <- sample_multinom_rows(prob_A, labels = pars$hospital_labels)
    h_score <- pars$hospital_scores[match(as.character(A), pars$hospital_labels)]
    gamma_ZA <- if (!is.null(pars$gamma_ZA)) pars$gamma_ZA else 0
    lambda_Z <- pars$lambda_Z_calibrated
    if (is.null(lambda_Z)) {
        lambda_Z <- pars$s_ZY * pars$race_direct_scores
    }
    lp_Y <- pars$intercept_Y +
        pars$s_XY * x_score_Y +
        pars$s_AY * h_score +
        lambda_Z[Z_idx] +
        gamma_ZA * z_direct * h_score
    p_Y <- plogis(lp_Y)
    Y <- rbinom(n, size = 1, prob = p_Y)
    dat <- data.frame(id = seq_len(n), Z = Z, A = A, Y = Y, X, stringsAsFactors = FALSE)
    attr(dat, "true_prob_A") <- prob_A
    attr(dat, "true_prob_Y") <- p_Y
    attr(dat, "true_x_score_A") <- x_score_A
    attr(dat, "true_x_score_Y") <- x_score_Y
    dat
}

true_prob_Z_given_X <- function(X, pars) {
    n <- nrow(X)
    R <- pars$R
    mu_list <- true_mu_X_given_Z(pars)
    log_num <- matrix(0, nrow = n, ncol = R)
    for (r in seq_len(R)) {
        log_px_given_z <- mvnorm_logdens_rows(X, mu = mu_list[[r]], Sigma = pars$Sigma_X)
        log_num[, r] <- log(pars$race_probs[r]) + log_px_given_z
    }
    row_max <- apply(log_num, 1, max)
    num <- exp(log_num - row_max)
    num/rowSums(num)
}

marginal_EY_by_Z_nested_mc <- function(pars, gamma_ZA, lambda_Z, mc_n = 200000, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    pars_eval <- pars
    pars_eval$gamma_ZA <- gamma_ZA
    pars_eval$lambda_Z_calibrated <- lambda_Z
    out <- numeric(pars$R)
    names(out) <- pars$race_labels
    for (r in seq_len(pars$R)) {
        mu_r <- pars$s_ZX * pars$race_severity_scores[r] * pars$mu_pattern
        X_r <- MASS::mvrnorm(n = mc_n, mu = mu_r, Sigma = pars$Sigma_X)
        colnames(X_r) <- paste0("X", seq_len(pars$p))
        Z_idx <- rep(r, mc_n)
        prob_A <- true_prob_A_given_XZ(X_r, Z_idx, pars)
        m_A <- true_mean_Y_given_XZA_allA_nested(X_r, Z_idx, pars_eval)
        out[r] <- mean(rowSums(prob_A * m_A))
    }
    out
}

calibrate_lambda_Z_for_ZA_interaction <- function(pars,
                                                  gamma_ZA,
                                                  target_EY_by_Z = NULL,
                                                  mc_n = 200000,
                                                  seed = 202606,
                                                  interval = c(-10, 10)) {
    lambda_Z_baseline <- pars$s_ZY * pars$race_direct_scores
    names(lambda_Z_baseline) <- pars$race_labels
    if (is.null(target_EY_by_Z)) {
        target_EY_by_Z <- marginal_EY_by_Z_nested_mc(
            pars = pars,
            gamma_ZA = 0,
            lambda_Z = lambda_Z_baseline,
            mc_n = mc_n,
            seed = seed
        )
    }
    target_EY_by_Z <- target_EY_by_Z[pars$race_labels]
    if (gamma_ZA == 0) {
        achieved_EY_by_Z <- marginal_EY_by_Z_nested_mc(
            pars = pars,
            gamma_ZA = 0,
            lambda_Z = lambda_Z_baseline,
            mc_n = mc_n,
            seed = seed
        )
        return(list(
            lambda_Z_calibrated = lambda_Z_baseline,
            target_EY_by_Z = target_EY_by_Z,
            achieved_EY_by_Z = achieved_EY_by_Z[pars$race_labels],
            gamma_ZA = gamma_ZA,
            mc_n = mc_n,
            seed = seed
        ))
    }

    set.seed(seed)
    mc_by_z <- vector("list", pars$R)
    for (r in seq_len(pars$R)) {
        mu_r <- pars$s_ZX * pars$race_severity_scores[r] * pars$mu_pattern
        X_r <- MASS::mvrnorm(n = mc_n, mu = mu_r, Sigma = pars$Sigma_X)
        colnames(X_r) <- paste0("X", seq_len(pars$p))
        Z_idx <- rep(r, mc_n)
        mc_by_z[[r]] <- list(
            X = X_r,
            Z_idx = Z_idx,
            prob_A = true_prob_A_given_XZ(X_r, Z_idx, pars)
        )
    }

    eval_group_mean <- function(r, lambda_r) {
        pars_eval <- pars
        lambda_eval <- lambda_Z_baseline
        lambda_eval[r] <- lambda_r
        pars_eval$gamma_ZA <- gamma_ZA
        pars_eval$lambda_Z_calibrated <- lambda_eval
        obj <- mc_by_z[[r]]
        m_A <- true_mean_Y_given_XZA_allA_nested(obj$X, obj$Z_idx, pars_eval)
        mean(rowSums(obj$prob_A * m_A))
    }

    lambda_Z_calibrated <- lambda_Z_baseline
    achieved_EY_by_Z <- numeric(pars$R)
    names(achieved_EY_by_Z) <- pars$race_labels
    for (r in seq_len(pars$R)) {
        target_r <- unname(target_EY_by_Z[pars$race_labels[r]])
        f_root <- function(lambda_r) eval_group_mean(r, lambda_r) - target_r
        lower_val <- f_root(interval[1])
        upper_val <- f_root(interval[2])
        if (lower_val * upper_val > 0) {
            stop(
                "Calibration interval does not bracket the target for group ",
                pars$race_labels[r],
                "."
            )
        }
        lambda_Z_calibrated[r] <- uniroot(f_root, interval = interval)$root
        achieved_EY_by_Z[r] <- eval_group_mean(r, lambda_Z_calibrated[r])
    }
    names(lambda_Z_calibrated) <- pars$race_labels
    list(
        lambda_Z_calibrated = lambda_Z_calibrated,
        target_EY_by_Z = target_EY_by_Z,
        achieved_EY_by_Z = achieved_EY_by_Z,
        gamma_ZA = gamma_ZA,
        mc_n = mc_n,
        seed = seed
    )
}

make_sim_params_nested <- function(...,
                                   gamma_ZA = 0,
                                   calibration_mc_n = 200000,
                                   calibration_seed = 202606) {
    pars <- make_sim_params(...)
    pars$gamma_ZA <- gamma_ZA
    calibration <- calibrate_lambda_Z_for_ZA_interaction(
        pars = pars,
        gamma_ZA = gamma_ZA,
        mc_n = calibration_mc_n,
        seed = calibration_seed
    )
    pars$lambda_Z_calibrated <- calibration$lambda_Z_calibrated
    pars$calibration_target_EY_by_Z <- calibration$target_EY_by_Z
    pars$calibration_achieved_EY_by_Z <- calibration$achieved_EY_by_Z
    pars$calibration_mc_n <- calibration_mc_n
    pars$calibration_seed <- calibration_seed
    pars$nested_outcome <- TRUE
    pars
}

true_conditional_objects <- function(X, pars) {
    n <- nrow(X)
    R <- pars$R
    H <- pars$H
    g_all <- matrix(0, nrow = n, ncol = R)
    vA_all <- matrix(0, nrow = n, ncol = R)
    rA_all <- matrix(0, nrow = n, ncol = R)
    probA_list <- vector("list", R)
    mA_list <- vector("list", R)
    for (r in seq_len(R)) {
        Z_idx <- rep(r, n)
        probA_r <- true_prob_A_given_XZ(X, Z_idx, pars)
        mA_r <- true_mean_Y_given_XZA_allA(X, Z_idx, pars)
        g_r <- rowSums(probA_r * mA_r)
        vA_r <- rowSums(probA_r * (mA_r - g_r)^2)
        rA_r <- rowSums(probA_r * (mA_r * (1 - mA_r)))
        g_all[, r] <- g_r
        vA_all[, r] <- vA_r
        rA_all[, r] <- rA_r
        probA_list[[r]] <- probA_r
        mA_list[[r]] <- mA_r
    }
    list(g_all = g_all, vA_all = vA_all, rA_all = rA_all, probA_list = probA_list, mA_list = mA_list)
}

true_topological_decomposition_mc <- function(pars, mc_n = 1e+05, seed = NULL) {
    if (!is.null(seed)) 
        set.seed(seed)
    dat_mc <- simulate_one(n = mc_n, pars = pars)
    X <- as.matrix(dat_mc[, paste0("X", seq_len(pars$p)), drop = FALSE])
    Z_idx <- match(dat_mc$Z, pars$race_labels)
    obj <- true_conditional_objects(X, pars)
    g_all <- obj$g_all
    vA_all <- obj$vA_all
    rA_all <- obj$rA_all
    u_z <- numeric(pars$R)
    for (r in seq_len(pars$R)) {
        idx <- which(Z_idx == r)
        u_z[r] <- mean(g_all[idx, r])
    }
    mu_u <- sum(pars$race_probs * u_z)
    delta_Z <- sum(pars$race_probs * (u_z - mu_u)^2)
    delta_X <- 0
    for (r in seq_len(pars$R)) {
        idx <- which(Z_idx == r)
        delta_X <- delta_X + pars$race_probs[r] * mean((g_all[idx, r] - u_z[r])^2)
    }
    delta_A <- 0
    for (r in seq_len(pars$R)) {
        idx <- which(Z_idx == r)
        delta_A <- delta_A + pars$race_probs[r] * mean(vA_all[idx, r])
    }
    delta_res <- 0
    for (r in seq_len(pars$R)) {
        idx <- which(Z_idx == r)
        delta_res <- delta_res + pars$race_probs[r] * mean(rA_all[idx, r])
    }
    pY_true <- attr(dat_mc, "true_prob_Y")
    varY_true <- mean(pY_true * (1 - pY_true)) + var(pY_true)
    out <- c(delta_Z = delta_Z, delta_X = delta_X, delta_A = delta_A, delta_res = delta_res, total = delta_Z + delta_X + delta_A + delta_res, varY_true = varY_true)
    list(components = out, mc_sample = dat_mc)
}

true_modified_decomposition_mc <- function(pars, mc_n = 1e+05, seed = NULL) {
    if (!is.null(seed)) 
        set.seed(seed)
    dat_mc <- simulate_one(n = mc_n, pars = pars)
    X <- as.matrix(dat_mc[, paste0("X", seq_len(pars$p)), drop = FALSE])
    probZ_given_X <- true_prob_Z_given_X(X, pars)
    obj <- true_conditional_objects(X, pars)
    g_all <- obj$g_all
    vA_all <- obj$vA_all
    rA_all <- obj$rA_all
    h_x <- rowSums(probZ_given_X * g_all)
    delta_X <- mean((h_x - mean(h_x))^2)
    delta_Z <- mean(rowSums(probZ_given_X * (g_all - h_x)^2))
    delta_A <- mean(rowSums(probZ_given_X * vA_all))
    delta_res <- mean(rowSums(probZ_given_X * rA_all))
    pY_true <- attr(dat_mc, "true_prob_Y")
    varY_true <- mean(pY_true * (1 - pY_true)) + var(pY_true)
    out <- c(delta_X = delta_X, delta_Z = delta_Z, delta_A = delta_A, delta_res = delta_res, total = delta_X + delta_Z + delta_A + delta_res, varY_true = varY_true)
    list(components = out, mc_sample = dat_mc, probZ_given_X = probZ_given_X)
}

true_decompositions_mc <- function(pars, mc_n = 1e+05, seed_top = 1001, seed_mod = 2001) {
    top <- true_topological_decomposition_mc(pars, mc_n = mc_n, seed = seed_top)
    mod <- true_modified_decomposition_mc(pars, mc_n = mc_n, seed = seed_mod)
    list(topological = top, modified = mod)
}

true_conditional_objects_nested <- function(X, pars) {
    n <- nrow(X)
    R <- pars$R
    g_all <- matrix(0, nrow = n, ncol = R)
    vA_all <- matrix(0, nrow = n, ncol = R)
    rA_all <- matrix(0, nrow = n, ncol = R)
    probA_list <- vector("list", R)
    mA_list <- vector("list", R)
    for (r in seq_len(R)) {
        Z_idx <- rep(r, n)
        probA_r <- true_prob_A_given_XZ(X, Z_idx, pars)
        mA_r <- true_mean_Y_given_XZA_allA_nested(X, Z_idx, pars)
        g_r <- rowSums(probA_r * mA_r)
        vA_r <- rowSums(probA_r * (mA_r - g_r)^2)
        rA_r <- rowSums(probA_r * (mA_r * (1 - mA_r)))
        g_all[, r] <- g_r
        vA_all[, r] <- vA_r
        rA_all[, r] <- rA_r
        probA_list[[r]] <- probA_r
        mA_list[[r]] <- mA_r
    }
    list(g_all = g_all, vA_all = vA_all, rA_all = rA_all, probA_list = probA_list, mA_list = mA_list)
}

true_topological_decomposition_mc_nested <- function(pars, mc_n = 1e+05, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    dat_mc <- simulate_one_nested(n = mc_n, pars = pars)
    X <- as.matrix(dat_mc[, paste0("X", seq_len(pars$p)), drop = FALSE])
    Z_idx <- match(dat_mc$Z, pars$race_labels)
    obj <- true_conditional_objects_nested(X, pars)
    g_all <- obj$g_all
    vA_all <- obj$vA_all
    rA_all <- obj$rA_all
    u_z <- numeric(pars$R)
    for (r in seq_len(pars$R)) {
        idx <- which(Z_idx == r)
        u_z[r] <- mean(g_all[idx, r])
    }
    mu_u <- sum(pars$race_probs * u_z)
    delta_Z <- sum(pars$race_probs * (u_z - mu_u)^2)
    delta_X <- 0
    for (r in seq_len(pars$R)) {
        idx <- which(Z_idx == r)
        delta_X <- delta_X + pars$race_probs[r] * mean((g_all[idx, r] - u_z[r])^2)
    }
    delta_A <- 0
    for (r in seq_len(pars$R)) {
        idx <- which(Z_idx == r)
        delta_A <- delta_A + pars$race_probs[r] * mean(vA_all[idx, r])
    }
    delta_res <- 0
    for (r in seq_len(pars$R)) {
        idx <- which(Z_idx == r)
        delta_res <- delta_res + pars$race_probs[r] * mean(rA_all[idx, r])
    }
    pY_true <- attr(dat_mc, "true_prob_Y")
    varY_true <- mean(pY_true * (1 - pY_true)) + var(pY_true)
    out <- c(delta_Z = delta_Z, delta_X = delta_X, delta_A = delta_A, delta_res = delta_res, total = delta_Z + delta_X + delta_A + delta_res, varY_true = varY_true)
    list(components = out, mc_sample = dat_mc)
}

true_modified_decomposition_mc_nested <- function(pars, mc_n = 1e+05, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    dat_mc <- simulate_one_nested(n = mc_n, pars = pars)
    X <- as.matrix(dat_mc[, paste0("X", seq_len(pars$p)), drop = FALSE])
    probZ_given_X <- true_prob_Z_given_X(X, pars)
    obj <- true_conditional_objects_nested(X, pars)
    g_all <- obj$g_all
    vA_all <- obj$vA_all
    rA_all <- obj$rA_all
    h_x <- rowSums(probZ_given_X * g_all)
    delta_X <- mean((h_x - mean(h_x))^2)
    delta_Z <- mean(rowSums(probZ_given_X * (g_all - h_x)^2))
    delta_A <- mean(rowSums(probZ_given_X * vA_all))
    delta_res <- mean(rowSums(probZ_given_X * rA_all))
    pY_true <- attr(dat_mc, "true_prob_Y")
    varY_true <- mean(pY_true * (1 - pY_true)) + var(pY_true)
    out <- c(delta_X = delta_X, delta_Z = delta_Z, delta_A = delta_A, delta_res = delta_res, total = delta_X + delta_Z + delta_A + delta_res, varY_true = varY_true)
    list(components = out, mc_sample = dat_mc, probZ_given_X = probZ_given_X)
}

true_decompositions_mc_nested <- function(pars, mc_n = 1e+05, seed_top = 1001, seed_mod = 2001) {
    top <- true_topological_decomposition_mc_nested(pars, mc_n = mc_n, seed = seed_top)
    mod <- true_modified_decomposition_mc_nested(pars, mc_n = mc_n, seed = seed_mod)
    list(topological = top, modified = mod)
}

predict_multinom_prob <- function(fit, newdata, levels_out = NULL) {
    pr <- predict(fit, newdata = newdata, type = "probs")
    if (is.null(dim(pr))) {
        if (!is.null(fit$lev) && length(fit$lev) == 2 && length(pr) == nrow(newdata)) {
            pr <- cbind(1 - pr, pr)
            colnames(pr) <- fit$lev
        }
        else if (is.null(levels_out)) {
            pr <- matrix(pr, nrow = 1)
        }
        else {
            tmp <- matrix(0, nrow = 1, ncol = length(levels_out))
            colnames(tmp) <- levels_out
            nm <- names(pr)
            tmp[1, nm] <- pr
            pr <- tmp
        }
    }
    else {
        pr <- as.matrix(pr)
    }
    if (!is.null(levels_out)) {
        tmp <- matrix(0, nrow = nrow(pr), ncol = length(levels_out))
        colnames(tmp) <- levels_out
        common <- intersect(colnames(pr), levels_out)
        tmp[, common] <- pr[, common, drop = FALSE]
        missing <- setdiff(levels_out, colnames(pr))
        if (length(missing) == 1) {
            tmp[, missing] <- 1 - rowSums(tmp[, common, drop = FALSE])
        }
        pr <- tmp
    }
    pr
}

fit_decomp_models <- function(dat, pars, fit_X_given_Z = TRUE) {
    x_names <- paste0("X", seq_len(pars$p))
    dat <- dat
    dat$Z <- relevel(dat$Z, ref = pars$baseline_race)
    dat$A <- relevel(dat$A, ref = pars$baseline_hospital)
    pZ_hat <- prop.table(table(dat$Z))
    pZ_hat <- pZ_hat[pars$race_labels]
    form_Z <- as.formula(paste("Z ~", paste(x_names, collapse = " + ")))
    fit_Z_given_X <- tryCatch(nnet::multinom(form_Z, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    form_A <- as.formula(paste("A ~ Z +", paste(x_names, collapse = " + ")))
    fit_A_given_XZ <- tryCatch(nnet::multinom(form_A, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    form_Y <- as.formula(paste("Y ~ A + Z +", paste(x_names, collapse = " + ")))
    fit_Y_given_XZA <- tryCatch(glm(form_Y, data = dat, family = binomial()), error = function(e) e)
    fit_X_given_Z_obj <- NULL
    if (fit_X_given_Z) {
        mu_hat <- sapply(split(dat[, x_names, drop = FALSE], dat$Z), colMeans)
        mu_hat <- t(mu_hat)
        mu_hat <- mu_hat[pars$race_labels, , drop = FALSE]
        Sigma_num <- matrix(0, nrow = pars$p, ncol = pars$p)
        denom <- 0
        for (g in levels(dat$Z)) {
            Xg <- as.matrix(dat[dat$Z == g, x_names, drop = FALSE])
            if (nrow(Xg) > 1) {
                cg <- colMeans(Xg)
                centered <- sweep(Xg, 2, cg, "-")
                Sigma_num <- Sigma_num + crossprod(centered)
                denom <- denom + (nrow(Xg) - 1)
            }
        }
        Sigma_hat <- Sigma_num/denom
        rownames(Sigma_hat) <- colnames(Sigma_hat) <- x_names
        fit_X_given_Z_obj <- list(mu_hat = mu_hat, Sigma_hat = Sigma_hat)
    }
    list(x_names = x_names, pZ_hat = pZ_hat, fit_Z_given_X = fit_Z_given_X, fit_A_given_XZ = fit_A_given_XZ, fit_Y_given_XZA = fit_Y_given_XZA, fit_X_given_Z = fit_X_given_Z_obj)
}

build_conditional_surfaces <- function(Xmat, fits, pars) {
    n <- nrow(Xmat)
    R <- pars$R
    H <- pars$H
    x_names <- fits$x_names
    stopifnot(all(colnames(Xmat) == x_names))
    m_array <- array(NA_real_, dim = c(n, R, H), dimnames = list(NULL, pars$race_labels, pars$hospital_labels))
    pA_array <- array(NA_real_, dim = c(n, R, H), dimnames = list(NULL, pars$race_labels, pars$hospital_labels))
    for (r in seq_len(R)) {
        zlab <- pars$race_labels[r]
        newA <- data.frame(Xmat, Z = factor(zlab, levels = pars$race_labels))
        pA_r <- predict_multinom_prob(fits$fit_A_given_XZ, newdata = newA, levels_out = pars$hospital_labels)
        pA_array[, r, ] <- pA_r
        for (h in seq_len(H)) {
            alab <- pars$hospital_labels[h]
            newY <- data.frame(Xmat, Z = factor(zlab, levels = pars$race_labels), A = factor(alab, levels = pars$hospital_labels))
            m_array[, r, h] <- predict(fits$fit_Y_given_XZA, newdata = newY, type = "response")
        }
    }
    g_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    vA_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    rA_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    for (r in seq_len(R)) {
        pA_r <- pA_array[, r, ]
        m_r <- m_array[, r, ]
        g_r <- rowSums(pA_r * m_r)
        v_r <- rowSums(pA_r * (m_r - g_r)^2)
        rr_r <- rowSums(pA_r * (m_r * (1 - m_r)))
        g_mat[, r] <- g_r
        vA_mat[, r] <- v_r
        rA_mat[, r] <- rr_r
    }
    list(m_array = m_array, pA_array = pA_array, g_mat = g_mat, vA_mat = vA_mat, rA_mat = rA_mat)
}

fitted_prob_Z_given_X_bayes <- function(Xmat, fits, pars) {
    stopifnot(!is.null(fits$fit_X_given_Z))
    R <- pars$R
    mu_hat <- fits$fit_X_given_Z$mu_hat
    Sigma_hat <- fits$fit_X_given_Z$Sigma_hat
    pZ_hat <- as.numeric(fits$pZ_hat)
    names(pZ_hat) <- names(fits$pZ_hat)
    log_num <- matrix(0, nrow = nrow(Xmat), ncol = R, dimnames = list(NULL, pars$race_labels))
    for (r in seq_len(R)) {
        zlab <- pars$race_labels[r]
        log_px_given_z <- mvnorm_logdens_rows(Xmat, mu = mu_hat[zlab, ], Sigma = Sigma_hat)
        log_num[, r] <- log(pZ_hat[zlab]) + log_px_given_z
    }
    row_max <- apply(log_num, 1, max)
    num <- exp(log_num - row_max)
    num/rowSums(num)
}

# Approximate uncertainty is added after model fitting.  Bayesian bootstrap
# weights propagate uncertainty in empirical integration distributions, while
# Gaussian draws from the fitted covariance matrices propagate parametric model
# uncertainty.  We intentionally do not refit models under BB weights: the
# likelihood factors are treated as independent post-estimation components.
draw_bb_weights <- function(n) {
    w <- stats::rexp(n, rate = 1)
    w <- w/sum(w)
    if (!isTRUE(all.equal(sum(w), 1, tolerance = 1e-10))) {
        stop("Bayesian bootstrap weights do not sum to 1.")
    }
    w
}

weighted_level_probs <- function(W, f, levels) {
    if (!isTRUE(all.equal(sum(W), 1, tolerance = 1e-10))) {
        stop("Bayesian bootstrap weights must sum to 1.")
    }
    out <- stats::setNames(rep(0, length(levels)), levels)
    tab <- tapply(W, factor(f, levels = levels), sum)
    tab[is.na(tab)] <- 0
    out[names(tab)] <- as.numeric(tab)
    if (!isTRUE(all.equal(sum(out), 1, tolerance = 1e-10))) {
        stop("Weighted empirical level probabilities do not sum to 1.")
    }
    out
}

draw_mvn_approx <- function(mu, Sigma) {
    ev <- eigen(Sigma, symmetric = TRUE)
    vals <- pmax(ev$values, 0)
    as.numeric(mu + ev$vectors %*% (sqrt(vals) * stats::rnorm(length(mu))))
}

draw_glm_coef <- function(fit) {
    theta <- stats::coef(fit)
    vc <- stats::vcov(fit)
    out <- draw_mvn_approx(theta, vc)
    names(out) <- names(theta)
    out
}

draw_multinom_coef <- function(fit) {
    theta <- stats::coef(fit)
    theta_mat <- if (is.null(dim(theta))) {
        matrix(theta, nrow = 1, dimnames = list(fit$lev[-1], names(theta)))
    } else {
        as.matrix(theta)
    }
    vc <- stats::vcov(fit)
    theta_vec <- as.vector(t(theta_mat))
    names(theta_vec) <- if (nrow(theta_mat) == 1) {
        colnames(theta_mat)
    } else {
        as.vector(t(outer(rownames(theta_mat), colnames(theta_mat), paste, sep = ":")))
    }
    vc_names <- rownames(vc)
    if (is.null(vc_names) || is.null(colnames(vc)) || !identical(vc_names, colnames(vc))) {
        stop("Multinomial covariance matrix must have matching row and column names.")
    }
    missing <- setdiff(vc_names, names(theta_vec))
    if (length(missing) > 0) {
        stop("Fitted multinomial coefficients are missing covariance entries: ", paste(missing, collapse = ", "), ".")
    }
    theta_vec <- theta_vec[vc_names]
    draw <- draw_mvn_approx(theta_vec, vc[vc_names, vc_names, drop = FALSE])
    draw_mat <- matrix(draw, nrow = nrow(theta_mat), ncol = ncol(theta_mat), byrow = TRUE, dimnames = dimnames(theta_mat))
    draw_mat
}

fitted_multinom_coef_matrix <- function(fit) {
    theta <- stats::coef(fit)
    if (is.null(dim(theta))) {
        matrix(theta, nrow = 1, dimnames = list(fit$lev[-1], names(theta)))
    } else {
        as.matrix(theta)
    }
}

predict_glm_prob_draw <- function(fit, newdata, theta) {
    x <- stats::model.matrix(stats::delete.response(stats::terms(fit)), newdata, xlev = fit$xlevels)
    theta <- theta[colnames(x)]
    stats::plogis(as.vector(x %*% theta))
}

predict_multinom_prob_draw <- function(fit, newdata, theta_mat, levels_out = NULL) {
    x <- stats::model.matrix(stats::delete.response(stats::terms(fit)), newdata, xlev = fit$xlevels)
    theta_mat <- theta_mat[, colnames(x), drop = FALSE]
    eta_nonbase <- x %*% t(theta_mat)
    lev <- fit$lev
    eta <- matrix(0, nrow = nrow(x), ncol = length(lev), dimnames = list(NULL, lev))
    eta[, rownames(theta_mat)] <- eta_nonbase
    pr <- softmax_mat(eta)
    if (!is.null(levels_out)) {
        tmp <- matrix(0, nrow = nrow(pr), ncol = length(levels_out), dimnames = list(NULL, levels_out))
        tmp[, intersect(colnames(pr), levels_out)] <- pr[, intersect(colnames(pr), levels_out), drop = FALSE]
        pr <- tmp
    }
    pr
}

check_multinom_draw_prediction_mean <- function(fit, newdata, levels_out, tol = 1e-8) {
    theta_mat <- fitted_multinom_coef_matrix(fit)
    pr_draw <- predict_multinom_prob_draw(fit, newdata = newdata, theta_mat = theta_mat, levels_out = levels_out)
    pr_fit <- predict_multinom_prob(fit, newdata = newdata, levels_out = levels_out)
    max_diff <- max(abs(pr_draw - pr_fit))
    if (!is.finite(max_diff) || max_diff > tol) {
        stop("Draw-based multinomial prediction does not match fitted prediction; max difference = ", max_diff, ".")
    }
    invisible(TRUE)
}

draw_X_given_Z_bb <- function(dat, fits, pars, W) {
    x_names <- fits$x_names
    Xmat <- as.matrix(dat[, x_names, drop = FALSE])
    Z_obs <- dat$Z
    pZ_hat <- weighted_level_probs(W, Z_obs, pars$race_labels)
    mu_hat <- matrix(NA_real_, nrow = pars$R, ncol = pars$p, dimnames = list(pars$race_labels, x_names))
    Sigma_hat <- matrix(0, nrow = pars$p, ncol = pars$p, dimnames = list(x_names, x_names))
    for (r in seq_len(pars$R)) {
        zlab <- pars$race_labels[r]
        idx <- which(Z_obs == zlab)
        if (length(idx) == 0 || pZ_hat[zlab] <= 0) {
            stop("Cannot compute within-stratum Bayesian bootstrap quantity: stratum has zero observations.")
        }
        wz <- W[idx]
        wz <- wz/sum(wz)
        mu_hat[zlab, ] <- colSums(Xmat[idx, , drop = FALSE] * wz)
        centered <- sweep(Xmat[idx, , drop = FALSE], 2, mu_hat[zlab, ], "-")
        Sigma_hat <- Sigma_hat + crossprod(centered, centered * W[idx])
    }
    Sigma_hat <- Sigma_hat + diag(1e-8, nrow(Sigma_hat))
    list(mu_hat = mu_hat, Sigma_hat = Sigma_hat, pZ_hat = pZ_hat)
}

build_conditional_surfaces_uncertainty <- function(Xmat, fits, pars, theta_A, theta_Y) {
    n <- nrow(Xmat)
    R <- pars$R
    H <- pars$H
    x_names <- fits$x_names
    stopifnot(all(colnames(Xmat) == x_names))
    m_array <- array(NA_real_, dim = c(n, R, H), dimnames = list(NULL, pars$race_labels, pars$hospital_labels))
    pA_array <- array(NA_real_, dim = c(n, R, H), dimnames = list(NULL, pars$race_labels, pars$hospital_labels))
    for (r in seq_len(R)) {
        zlab <- pars$race_labels[r]
        newA <- data.frame(Xmat, Z = factor(zlab, levels = pars$race_labels))
        pA_r <- predict_multinom_prob_draw(fits$fit_A_given_XZ, newdata = newA, theta_mat = theta_A, levels_out = pars$hospital_labels)
        pA_array[, r, ] <- pA_r
        for (h in seq_len(H)) {
            alab <- pars$hospital_labels[h]
            newY <- data.frame(Xmat, Z = factor(zlab, levels = pars$race_labels), A = factor(alab, levels = pars$hospital_labels))
            m_array[, r, h] <- predict_glm_prob_draw(fits$fit_Y_given_XZA, newdata = newY, theta = theta_Y)
        }
    }
    g_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    vA_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    rA_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    for (r in seq_len(R)) {
        pA_r <- pA_array[, r, ]
        m_r <- m_array[, r, ]
        g_r <- rowSums(pA_r * m_r)
        g_mat[, r] <- g_r
        vA_mat[, r] <- rowSums(pA_r * (m_r - g_r)^2)
        rA_mat[, r] <- rowSums(pA_r * (m_r * (1 - m_r)))
    }
    list(m_array = m_array, pA_array = pA_array, g_mat = g_mat, vA_mat = vA_mat, rA_mat = rA_mat)
}

fitted_prob_Z_given_X_bayes_uncertainty <- function(Xmat, x_given_z, pars) {
    R <- pars$R
    mu_hat <- x_given_z$mu_hat
    Sigma_hat <- x_given_z$Sigma_hat
    pZ_hat <- x_given_z$pZ_hat
    log_num <- matrix(0, nrow = nrow(Xmat), ncol = R, dimnames = list(NULL, pars$race_labels))
    for (r in seq_len(R)) {
        zlab <- pars$race_labels[r]
        log_px_given_z <- mvnorm_logdens_rows(Xmat, mu = mu_hat[zlab, ], Sigma = Sigma_hat)
        log_num[, r] <- log(pZ_hat[zlab]) + log_px_given_z
    }
    row_max <- apply(log_num, 1, max)
    num <- exp(log_num - row_max)
    num/rowSums(num)
}

weighted_topological_decomp <- function(dat, fits, pars, surf, W) {
    Z_obs <- dat$Z
    g_mat <- surf$g_mat
    vA_mat <- surf$vA_mat
    rA_mat <- surf$rA_mat
    pZ_hat <- weighted_level_probs(W, Z_obs, pars$race_labels)
    u_z <- numeric(pars$R)
    names(u_z) <- pars$race_labels
    for (r in seq_len(pars$R)) {
        zlab <- pars$race_labels[r]
        idx <- which(Z_obs == zlab)
        if (length(idx) == 0 || pZ_hat[zlab] <= 0) {
            stop("Cannot compute within-stratum Bayesian bootstrap quantity: stratum has zero observations.")
        }
        u_z[zlab] <- sum(W[idx] * g_mat[idx, zlab])/pZ_hat[zlab]
    }
    mu_u <- sum(pZ_hat * u_z)
    delta_Z <- sum(pZ_hat * (u_z - mu_u)^2)
    delta_X <- delta_A <- delta_res <- 0
    for (r in seq_len(pars$R)) {
        zlab <- pars$race_labels[r]
        idx <- which(Z_obs == zlab)
        delta_X <- delta_X + sum(W[idx] * (g_mat[idx, zlab] - u_z[zlab])^2)
        delta_A <- delta_A + sum(W[idx] * vA_mat[idx, zlab])
        delta_res <- delta_res + sum(W[idx] * rA_mat[idx, zlab])
    }
    c(delta_Z = delta_Z, delta_X = delta_X, delta_A = delta_A, delta_res = delta_res, total = delta_Z + delta_X + delta_A + delta_res)
}

weighted_modified_decomp <- function(pZ_given_X, surf, W) {
    g_mat <- surf$g_mat
    vA_mat <- surf$vA_mat
    rA_mat <- surf$rA_mat
    h_x <- rowSums(pZ_given_X * g_mat)
    h_bar <- sum(W * h_x)
    delta_X <- sum(W * (h_x - h_bar)^2)
    delta_Z <- sum(W * rowSums(pZ_given_X * (g_mat - h_x)^2))
    delta_A <- sum(W * rowSums(pZ_given_X * vA_mat))
    delta_res <- sum(W * rowSums(pZ_given_X * rA_mat))
    c(delta_X = delta_X, delta_Z = delta_Z, delta_A = delta_A, delta_res = delta_res, total = delta_X + delta_Z + delta_A + delta_res)
}

check_decomposition_draws <- function(draws) {
    if (anyNA(draws)) {
        stop("Posterior decomposition draws contain NA values.")
    }
    component_cols <- setdiff(colnames(draws), "total")
    total_diff <- draws[, "total"] - rowSums(draws[, component_cols, drop = FALSE])
    if (any(abs(total_diff) > 1e-8)) {
        stop("Posterior decomposition totals do not equal the sum of components.")
    }
    invisible(TRUE)
}

summarize_decomposition_uncertainty <- function(draws) {
    rbind(
        mean = colMeans(draws, na.rm = TRUE),
        sd = apply(draws, 2, stats::sd, na.rm = TRUE),
        q2.5 = apply(draws, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
        q97.5 = apply(draws, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    )
}

estimate_decomposition_uncertainty <- function(dat,
                                               pars,
                                               fits,
                                               estimator = c("topological", "modified_direct", "modified_bayes"),
                                               B_post = 1000) {
    estimator <- match.arg(estimator)
    if (!inherits(fits$fit_A_given_XZ, "multinom") || !inherits(fits$fit_Y_given_XZA, "glm")) {
        stop("Uncertainty estimation requires fitted A|X,Z multinomial and Y|X,Z,A GLM models.")
    }
    if (identical(estimator, "modified_direct") && !inherits(fits$fit_Z_given_X, "multinom")) {
        stop("Uncertainty estimation for modified_direct requires fitted Z|X multinomial model.")
    }
    if (identical(estimator, "modified_bayes") && is.null(fits$fit_X_given_Z)) {
        stop("Uncertainty estimation for modified_bayes requires fit_X_given_Z.")
    }
    x_names <- fits$x_names
    Xmat <- as.matrix(dat[, x_names, drop = FALSE])
    comp_names <- if (identical(estimator, "topological")) {
        c("delta_Z", "delta_X", "delta_A", "delta_res", "total")
    } else {
        c("delta_X", "delta_Z", "delta_A", "delta_res", "total")
    }
    draws <- matrix(NA_real_, nrow = B_post, ncol = length(comp_names), dimnames = list(NULL, comp_names))
    for (b in seq_len(B_post)) {
        W <- draw_bb_weights(nrow(dat))
        theta_A <- draw_multinom_coef(fits$fit_A_given_XZ)
        theta_Y <- draw_glm_coef(fits$fit_Y_given_XZA)
        surf <- build_conditional_surfaces_uncertainty(Xmat, fits, pars, theta_A = theta_A, theta_Y = theta_Y)
        if (identical(estimator, "topological")) {
            draws[b, ] <- weighted_topological_decomp(dat, fits, pars, surf, W)[comp_names]
        } else if (identical(estimator, "modified_direct")) {
            theta_Z <- draw_multinom_coef(fits$fit_Z_given_X)
            pZ_given_X <- predict_multinom_prob_draw(fits$fit_Z_given_X, newdata = data.frame(Xmat), theta_mat = theta_Z, levels_out = pars$race_labels)
            draws[b, ] <- weighted_modified_decomp(pZ_given_X, surf, W)[comp_names]
        } else {
            x_given_z <- draw_X_given_Z_bb(dat, fits, pars, W)
            pZ_given_X <- fitted_prob_Z_given_X_bayes_uncertainty(Xmat, x_given_z, pars)
            draws[b, ] <- weighted_modified_decomp(pZ_given_X, surf, W)[comp_names]
        }
    }
    stopifnot(identical(colnames(draws), comp_names))
    out <- summarize_decomposition_uncertainty(draws)
    check_decomposition_draws(draws)
    if (any(out["q2.5", ] > out["q97.5", ])) {
        stop("Uncertainty interval endpoints are not ordered.")
    }
    attr(out, "draws") <- draws
    out
}

make_uncertainty_interval_df <- function(replicate, estimator, point_estimate, uncertainty, truth) {
    component <- colnames(uncertainty)
    truth_component <- truth[component]
    out <- data.frame(
        replicate = replicate,
        estimator = estimator,
        component = component,
        point_estimate = as.numeric(point_estimate[component]),
        posterior_mean = as.numeric(uncertainty["mean", component]),
        posterior_sd = as.numeric(uncertainty["sd", component]),
        q2.5 = as.numeric(uncertainty["q2.5", component]),
        q97.5 = as.numeric(uncertainty["q97.5", component]),
        truth = as.numeric(truth_component),
        stringsAsFactors = FALSE
    )
    if (any(out$q2.5 > out$q97.5)) {
        stop("Uncertainty interval endpoints are not ordered.")
    }
    out$covered <- out$q2.5 <= out$truth & out$truth <= out$q97.5
    if (!is.logical(out$covered)) {
        stop("Coverage indicator must be logical.")
    }
    out$interval_length <- out$q97.5 - out$q2.5
    out
}

summarize_uncertainty_performance <- function(interval_df) {
    components <- unique(interval_df$component)
    out <- do.call(rbind, lapply(components, function(component) {
        df <- interval_df[interval_df$component == component, , drop = FALSE]
        data.frame(
            component = component,
            mean_point_estimate = mean(df$point_estimate, na.rm = TRUE),
            mean_posterior_mean = mean(df$posterior_mean, na.rm = TRUE),
            empirical_bias_point = mean(df$point_estimate - df$truth, na.rm = TRUE),
            mean_posterior_sd = mean(df$posterior_sd, na.rm = TRUE),
            empirical_sd_point = stats::sd(df$point_estimate, na.rm = TRUE),
            coverage = mean(df$covered, na.rm = TRUE),
            mean_interval_length = mean(df$interval_length, na.rm = TRUE),
            rmse_point = sqrt(mean((df$point_estimate - df$truth)^2, na.rm = TRUE)),
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    out
}

estimate_topological_decomp <- function(dat, fits, pars) {
    x_names <- fits$x_names
    Xmat <- as.matrix(dat[, x_names, drop = FALSE])
    surf <- build_conditional_surfaces(Xmat, fits, pars)
    g_mat <- surf$g_mat
    vA_mat <- surf$vA_mat
    rA_mat <- surf$rA_mat
    pZ_hat <- as.numeric(fits$pZ_hat)
    names(pZ_hat) <- pars$race_labels
    Z_obs <- dat$Z
    u_z <- numeric(pars$R)
    delta_X <- 0
    delta_A <- 0
    delta_res <- 0
    for (r in seq_len(pars$R)) {
        zlab <- pars$race_labels[r]
        idx <- which(Z_obs == zlab)
        g_r <- g_mat[idx, zlab]
        u_z[r] <- mean(g_r)
        pz_r <- unname(pZ_hat[zlab])
        delta_X <- delta_X + pz_r * mean((g_r - u_z[r])^2)
        delta_A <- delta_A + pz_r * mean(vA_mat[idx, zlab])
        delta_res <- delta_res + pz_r * mean(rA_mat[idx, zlab])
    }
    mu_u <- sum(pZ_hat * u_z)
    delta_Z <- sum(pZ_hat * (u_z - mu_u)^2)
    out <- c(delta_Z = delta_Z, delta_X = delta_X, delta_A = delta_A, delta_res = delta_res, total = delta_Z + delta_X + delta_A + delta_res)
    out
}

estimate_modified_decomp_direct <- function(dat, fits, pars) {
    x_names <- fits$x_names
    Xmat <- as.matrix(dat[, x_names, drop = FALSE])
    pZ_given_X <- predict_multinom_prob(fits$fit_Z_given_X, newdata = data.frame(Xmat), levels_out = pars$race_labels)
    surf <- build_conditional_surfaces(Xmat, fits, pars)
    g_mat <- surf$g_mat
    vA_mat <- surf$vA_mat
    rA_mat <- surf$rA_mat
    h_x <- rowSums(pZ_given_X * g_mat)
    delta_X <- mean((h_x - mean(h_x))^2)
    delta_Z <- mean(rowSums(pZ_given_X * (g_mat - h_x)^2))
    delta_A <- mean(rowSums(pZ_given_X * vA_mat))
    delta_res <- mean(rowSums(pZ_given_X * rA_mat))
    out <- c(delta_X = delta_X, delta_Z = delta_Z, delta_A = delta_A, delta_res = delta_res, total = delta_X + delta_Z + delta_A + delta_res)
    out
}

estimate_modified_decomp_bayes <- function(dat, fits, pars) {
    if (is.null(fits$fit_X_given_Z)) {
        stop("fit_X_given_Z is NULL; refit with fit_X_given_Z = TRUE.")
    }
    x_names <- fits$x_names
    Xmat <- as.matrix(dat[, x_names, drop = FALSE])
    pZ_given_X <- fitted_prob_Z_given_X_bayes(Xmat, fits, pars)
    surf <- build_conditional_surfaces(Xmat, fits, pars)
    g_mat <- surf$g_mat
    vA_mat <- surf$vA_mat
    rA_mat <- surf$rA_mat
    h_x <- rowSums(pZ_given_X * g_mat)
    delta_X <- mean((h_x - mean(h_x))^2)
    delta_Z <- mean(rowSums(pZ_given_X * (g_mat - h_x)^2))
    delta_A <- mean(rowSums(pZ_given_X * vA_mat))
    delta_res <- mean(rowSums(pZ_given_X * rA_mat))
    out <- c(delta_X = delta_X, delta_Z = delta_Z, delta_A = delta_A, delta_res = delta_res, total = delta_X + delta_Z + delta_A + delta_res)
    out
}

estimate_all_decomps <- function(dat, pars, fit_X_given_Z = TRUE) {
    fits <- fit_decomp_models(dat, pars, fit_X_given_Z = fit_X_given_Z)
    top_hat <- estimate_topological_decomp(dat, fits, pars)
    mod_hat_direct <- estimate_modified_decomp_direct(dat, fits, pars)
    mod_hat_bayes <- NULL
    if (fit_X_given_Z) {
        mod_hat_bayes <- estimate_modified_decomp_bayes(dat, fits, pars)
    }
    list(fits = fits, topological_hat = top_hat, modified_hat_direct = mod_hat_direct, modified_hat_bayes = mod_hat_bayes)
}

summarize_estimator_performance <- function(est_mat, truth_vec) {
    est_mean <- colMeans(est_mat, na.rm = TRUE)
    est_sd <- apply(est_mat, 2, sd, na.rm = TRUE)
    bias <- est_mean - truth_vec[colnames(est_mat)]
    rmse <- sqrt(colMeans((sweep(est_mat, 2, truth_vec[colnames(est_mat)], "-"))^2, na.rm = TRUE))
    out <- rbind(truth = truth_vec[colnames(est_mat)], mean = est_mean, bias = bias, sd = est_sd, rmse = rmse)
    out
}

run_estimator_simulation <- function(B = 200,
                                     n = 3000,
                                     pars,
                                     true_vals = NULL,
                                     fit_X_given_Z = TRUE,
                                     seed = 1,
                                     verbose = TRUE) {
    set.seed(seed)
    if (is.null(true_vals)) {
        true_vals <- true_decompositions_mc(pars = pars, mc_n = 1e+05)
    }
    truth_top <- true_vals$topological$components[c("delta_Z", "delta_X", "delta_A", "delta_res", "total")]
    truth_mod <- true_vals$modified$components[c("delta_X", "delta_Z", "delta_A", "delta_res", "total")]
    top_store <- matrix(NA_real_, nrow = B, ncol = length(truth_top))
    colnames(top_store) <- names(truth_top)
    mod_dir_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
    colnames(mod_dir_store) <- names(truth_mod)
    mod_bayes_store <- NULL
    if (fit_X_given_Z) {
        mod_bayes_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
        colnames(mod_bayes_store) <- names(truth_mod)
    }
    conv_Z <- logical(B)
    conv_A <- logical(B)
    conv_Y <- logical(B)
    for (b in seq_len(B)) {
        dat_b <- simulate_one(n = n, pars = pars)
        est_b <- estimate_all_decomps(dat = dat_b, pars = pars, fit_X_given_Z = fit_X_given_Z)
        top_store[b, ] <- est_b$topological_hat[colnames(top_store)]
        mod_dir_store[b, ] <- est_b$modified_hat_direct[colnames(mod_dir_store)]
        if (fit_X_given_Z) {
            mod_bayes_store[b, ] <- est_b$modified_hat_bayes[colnames(mod_bayes_store)]
        }
        conv_Z[b] <- inherits(est_b$fits$fit_Z_given_X, "multinom")
        conv_A[b] <- inherits(est_b$fits$fit_A_given_XZ, "multinom")
        conv_Y[b] <- inherits(est_b$fits$fit_Y_given_XZA, "glm")
        if (verbose && (b%%max(1, floor(B/10)) == 0)) {
            message("Completed estimator replicate ", b, " / ", B)
        }
    }
    summ_top <- summarize_estimator_performance(top_store, truth_top)
    summ_mod_dir <- summarize_estimator_performance(mod_dir_store, truth_mod)
    summ_mod_bayes <- NULL
    if (fit_X_given_Z) {
        summ_mod_bayes <- summarize_estimator_performance(mod_bayes_store, truth_mod)
    }
    list(params = c(pars, list(n = n)), truth = list(topological = truth_top, modified = truth_mod), estimates = list(topological = top_store, modified_direct = mod_dir_store, modified_bayes = mod_bayes_store), summary = list(topological = summ_top, modified_direct = summ_mod_dir, modified_bayes = summ_mod_bayes), convergence = data.frame(conv_Z = conv_Z, conv_A = conv_A, conv_Y = conv_Y))
}

fit_decomp_models_nested <- function(dat,
                                     pars,
                                     fit_X_given_Z = TRUE,
                                     outcome_interaction = c("none", "ZA")) {
    outcome_interaction <- match.arg(outcome_interaction)
    x_names <- paste0("X", seq_len(pars$p))
    dat <- dat
    dat$Z <- relevel(dat$Z, ref = pars$baseline_race)
    dat$A <- relevel(dat$A, ref = pars$baseline_hospital)
    pZ_hat <- prop.table(table(dat$Z))
    pZ_hat <- pZ_hat[pars$race_labels]
    form_Z <- as.formula(paste("Z ~", paste(x_names, collapse = " + ")))
    fit_Z_given_X <- tryCatch(nnet::multinom(form_Z, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    form_A <- as.formula(paste("A ~ Z +", paste(x_names, collapse = " + ")))
    fit_A_given_XZ <- tryCatch(nnet::multinom(form_A, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    if (identical(outcome_interaction, "ZA")) {
        form_Y <- as.formula(paste("Y ~ A * Z +", paste(x_names, collapse = " + ")))
    } else {
        form_Y <- as.formula(paste("Y ~ A + Z +", paste(x_names, collapse = " + ")))
    }
    fit_Y_warnings <- character()
    fit_Y_given_XZA <- tryCatch(
        withCallingHandlers(
            glm(form_Y, data = dat, family = binomial()),
            warning = function(w) {
                fit_Y_warnings <<- c(fit_Y_warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) e
    )
    fit_X_given_Z_obj <- NULL
    if (fit_X_given_Z) {
        mu_hat <- sapply(split(dat[, x_names, drop = FALSE], dat$Z), colMeans)
        mu_hat <- t(mu_hat)
        mu_hat <- mu_hat[pars$race_labels, , drop = FALSE]
        Sigma_num <- matrix(0, nrow = pars$p, ncol = pars$p)
        denom <- 0
        for (g in levels(dat$Z)) {
            Xg <- as.matrix(dat[dat$Z == g, x_names, drop = FALSE])
            if (nrow(Xg) > 1) {
                cg <- colMeans(Xg)
                centered <- sweep(Xg, 2, cg, "-")
                Sigma_num <- Sigma_num + crossprod(centered)
                denom <- denom + (nrow(Xg) - 1)
            }
        }
        Sigma_hat <- Sigma_num/denom
        rownames(Sigma_hat) <- colnames(Sigma_hat) <- x_names
        fit_X_given_Z_obj <- list(mu_hat = mu_hat, Sigma_hat = Sigma_hat)
    }
    list(
        x_names = x_names,
        pZ_hat = pZ_hat,
        fit_Z_given_X = fit_Z_given_X,
        fit_A_given_XZ = fit_A_given_XZ,
        fit_Y_given_XZA = fit_Y_given_XZA,
        fit_Y_warnings = fit_Y_warnings,
        fit_X_given_Z = fit_X_given_Z_obj,
        outcome_interaction = outcome_interaction
    )
}

fit_decomp_models_nested_firth <- function(dat, pars, fit_X_given_Z = TRUE) {
    if (!requireNamespace("logistf", quietly = TRUE)) {
        stop("The logistf package is required for fit_decomp_models_nested_firth().")
    }
    x_names <- paste0("X", seq_len(pars$p))
    dat <- dat
    dat$Z <- relevel(dat$Z, ref = pars$baseline_race)
    dat$A <- relevel(dat$A, ref = pars$baseline_hospital)
    pZ_hat <- prop.table(table(dat$Z))
    pZ_hat <- pZ_hat[pars$race_labels]
    form_Z <- as.formula(paste("Z ~", paste(x_names, collapse = " + ")))
    fit_Z_given_X <- tryCatch(nnet::multinom(form_Z, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    form_A <- as.formula(paste("A ~ Z +", paste(x_names, collapse = " + ")))
    fit_A_given_XZ <- tryCatch(nnet::multinom(form_A, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    form_Y <- as.formula(paste("Y ~ A * Z +", paste(x_names, collapse = " + ")))
    fit_Y_warnings <- character()
    fit_Y_given_XZA <- tryCatch(
        withCallingHandlers(
            logistf::logistf(form_Y, data = dat),
            warning = function(w) {
                fit_Y_warnings <<- c(fit_Y_warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) e
    )
    fit_X_given_Z_obj <- NULL
    if (fit_X_given_Z) {
        mu_hat <- sapply(split(dat[, x_names, drop = FALSE], dat$Z), colMeans)
        mu_hat <- t(mu_hat)
        mu_hat <- mu_hat[pars$race_labels, , drop = FALSE]
        Sigma_num <- matrix(0, nrow = pars$p, ncol = pars$p)
        denom <- 0
        for (g in levels(dat$Z)) {
            Xg <- as.matrix(dat[dat$Z == g, x_names, drop = FALSE])
            if (nrow(Xg) > 1) {
                cg <- colMeans(Xg)
                centered <- sweep(Xg, 2, cg, "-")
                Sigma_num <- Sigma_num + crossprod(centered)
                denom <- denom + (nrow(Xg) - 1)
            }
        }
        Sigma_hat <- Sigma_num/denom
        rownames(Sigma_hat) <- colnames(Sigma_hat) <- x_names
        fit_X_given_Z_obj <- list(mu_hat = mu_hat, Sigma_hat = Sigma_hat)
    }
    list(
        x_names = x_names,
        pZ_hat = pZ_hat,
        fit_Z_given_X = fit_Z_given_X,
        fit_A_given_XZ = fit_A_given_XZ,
        fit_Y_given_XZA = fit_Y_given_XZA,
        fit_Y_warnings = fit_Y_warnings,
        fit_X_given_Z = fit_X_given_Z_obj,
        outcome_interaction = "ZA_firth"
    )
}

outcome_ZA_cell_diagnostics_nested <- function(dat, pars) {
    cell_n <- table(
        factor(dat$Z, levels = pars$race_labels),
        factor(dat$A, levels = pars$hospital_labels)
    )
    cell_events <- tapply(
        dat$Y,
        list(
            factor(dat$Z, levels = pars$race_labels),
            factor(dat$A, levels = pars$hospital_labels)
        ),
        sum
    )
    cell_events[is.na(cell_events)] <- 0
    cell_rate <- cell_events / cell_n
    cell_rate[cell_n == 0] <- NA_real_
    c(
        min_ZA_cell_n = min(as.numeric(cell_n)),
        median_ZA_cell_n = stats::median(as.numeric(cell_n)),
        zero_event_ZA_cells = sum(cell_n > 0 & cell_events == 0),
        all_event_ZA_cells = sum(cell_n > 0 & cell_events == cell_n),
        min_ZA_cell_event_rate = min(cell_rate, na.rm = TRUE),
        max_ZA_cell_event_rate = max(cell_rate, na.rm = TRUE)
    )
}

glm_outcome_diagnostics_nested <- function(fit, warnings = character(), extreme_tol = 1e-6, big_coef_threshold = 10) {
    if (!inherits(fit, "glm")) {
        return(c(
            fit_Y_ok = FALSE,
            fit_Y_converged = FALSE,
            fit_Y_boundary = NA,
            fit_Y_iterations = NA,
            fit_Y_warning_count = length(warnings),
            fit_Y_separation_warning = any(grepl("fitted probabilities numerically 0 or 1|algorithm did not converge", warnings)),
            fit_Y_max_abs_coef = NA,
            fit_Y_min_fitted = NA,
            fit_Y_max_fitted = NA,
            fit_Y_extreme_fitted = NA,
            fit_Y_big_coef = NA
        ))
    }
    fit_vals <- fitted(fit)
    coefs <- coef(fit)
    max_abs_coef <- suppressWarnings(max(abs(coefs), na.rm = TRUE))
    if (!is.finite(max_abs_coef)) {
        max_abs_coef <- NA_real_
    }
    c(
        fit_Y_ok = TRUE,
        fit_Y_converged = isTRUE(fit$converged),
        fit_Y_boundary = isTRUE(fit$boundary),
        fit_Y_iterations = fit$iter,
        fit_Y_warning_count = length(warnings),
        fit_Y_separation_warning = any(grepl("fitted probabilities numerically 0 or 1|algorithm did not converge", warnings)),
        fit_Y_max_abs_coef = max_abs_coef,
        fit_Y_min_fitted = min(fit_vals, na.rm = TRUE),
        fit_Y_max_fitted = max(fit_vals, na.rm = TRUE),
        fit_Y_extreme_fitted = any(fit_vals < extreme_tol | fit_vals > 1 - extreme_tol),
        fit_Y_big_coef = any(abs(coefs) > big_coef_threshold, na.rm = TRUE)
    )
}

logistf_outcome_diagnostics_nested <- function(fit, warnings = character(), extreme_tol = 1e-6, big_coef_threshold = 10) {
    if (!inherits(fit, "logistf")) {
        return(c(
            fit_Y_ok = FALSE,
            fit_Y_converged = FALSE,
            fit_Y_boundary = NA,
            fit_Y_iterations = NA,
            fit_Y_warning_count = length(warnings),
            fit_Y_separation_warning = FALSE,
            fit_Y_max_abs_coef = NA,
            fit_Y_min_fitted = NA,
            fit_Y_max_fitted = NA,
            fit_Y_extreme_fitted = NA,
            fit_Y_big_coef = NA
        ))
    }
    fit_vals <- fit$predict
    coefs <- coef(fit)
    max_abs_coef <- suppressWarnings(max(abs(coefs), na.rm = TRUE))
    if (!is.finite(max_abs_coef)) {
        max_abs_coef <- NA_real_
    }
    converged <- if (!is.null(fit$conv)) {
        all(is.finite(fit$conv)) && max(abs(fit$conv), na.rm = TRUE) < 1e-4
    } else {
        TRUE
    }
    c(
        fit_Y_ok = TRUE,
        fit_Y_converged = converged,
        fit_Y_boundary = FALSE,
        fit_Y_iterations = if (!is.null(fit$iter)) fit$iter else NA,
        fit_Y_warning_count = length(warnings),
        fit_Y_separation_warning = FALSE,
        fit_Y_max_abs_coef = max_abs_coef,
        fit_Y_min_fitted = min(fit_vals, na.rm = TRUE),
        fit_Y_max_fitted = max(fit_vals, na.rm = TRUE),
        fit_Y_extreme_fitted = any(fit_vals < extreme_tol | fit_vals > 1 - extreme_tol),
        fit_Y_big_coef = any(abs(coefs) > big_coef_threshold, na.rm = TRUE)
    )
}

outcome_model_diagnostics_nested <- function(fit, warnings = character()) {
    if (inherits(fit, "logistf")) {
        return(logistf_outcome_diagnostics_nested(fit, warnings = warnings))
    }
    glm_outcome_diagnostics_nested(fit, warnings = warnings)
}

model_diagnostics_nested <- function(dat, fits, pars, replicate, outcome_interaction, simulation_seed = NA_integer_) {
    fit_diag <- outcome_model_diagnostics_nested(
        fit = fits$fit_Y_given_XZA,
        warnings = fits$fit_Y_warnings
    )
    cell_diag <- outcome_ZA_cell_diagnostics_nested(dat, pars)
    out <- data.frame(
        replicate = replicate,
        simulation_seed = simulation_seed,
        outcome_interaction = outcome_interaction,
        conv_Z = inherits(fits$fit_Z_given_X, "multinom"),
        conv_A = inherits(fits$fit_A_given_XZ, "multinom"),
        conv_Y = inherits(fits$fit_Y_given_XZA, "glm") || inherits(fits$fit_Y_given_XZA, "logistf"),
        as.list(fit_diag),
        as.list(cell_diag),
        fit_Y_warnings = paste(fits$fit_Y_warnings, collapse = " | "),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    out
}

estimate_all_decomps_nested <- function(dat,
                                        pars,
                                        fit_X_given_Z = TRUE,
                                        outcome_interaction = c("none", "ZA")) {
    outcome_interaction <- match.arg(outcome_interaction)
    fits <- fit_decomp_models_nested(
        dat = dat,
        pars = pars,
        fit_X_given_Z = fit_X_given_Z,
        outcome_interaction = outcome_interaction
    )
    top_hat <- estimate_topological_decomp(dat, fits, pars)
    mod_hat_direct <- estimate_modified_decomp_direct(dat, fits, pars)
    mod_hat_bayes <- NULL
    if (fit_X_given_Z) {
        mod_hat_bayes <- estimate_modified_decomp_bayes(dat, fits, pars)
    }
    list(
        fits = fits,
        topological_hat = top_hat,
        modified_hat_direct = mod_hat_direct,
        modified_hat_bayes = mod_hat_bayes,
        outcome_interaction = outcome_interaction
    )
}

estimate_all_decomps_nested_firth <- function(dat, pars, fit_X_given_Z = TRUE) {
    fits <- fit_decomp_models_nested_firth(
        dat = dat,
        pars = pars,
        fit_X_given_Z = fit_X_given_Z
    )
    top_names <- c("delta_Z", "delta_X", "delta_A", "delta_res", "total")
    mod_names <- c("delta_X", "delta_Z", "delta_A", "delta_res", "total")
    if (!inherits(fits$fit_Z_given_X, "multinom") ||
        !inherits(fits$fit_A_given_XZ, "multinom") ||
        !inherits(fits$fit_Y_given_XZA, "logistf")) {
        top_hat <- stats::setNames(rep(NA_real_, length(top_names)), top_names)
        mod_hat_direct <- stats::setNames(rep(NA_real_, length(mod_names)), mod_names)
        mod_hat_bayes <- if (fit_X_given_Z) stats::setNames(rep(NA_real_, length(mod_names)), mod_names) else NULL
    } else {
        top_hat <- estimate_topological_decomp(dat, fits, pars)
        mod_hat_direct <- estimate_modified_decomp_direct(dat, fits, pars)
        mod_hat_bayes <- NULL
        if (fit_X_given_Z) {
            mod_hat_bayes <- estimate_modified_decomp_bayes(dat, fits, pars)
        }
    }
    list(
        fits = fits,
        topological_hat = top_hat,
        modified_hat_direct = mod_hat_direct,
        modified_hat_bayes = mod_hat_bayes,
        outcome_interaction = "ZA_firth"
    )
}

run_estimator_simulation_nested <- function(B = 200,
                                            n = 3000,
                                            pars,
                                            true_vals = NULL,
                                            fit_X_given_Z = TRUE,
                                            outcome_interaction = c("none", "ZA"),
                                            seed = 1,
                                            verbose = TRUE) {
    outcome_interaction <- match.arg(outcome_interaction)
    set.seed(seed)
    if (is.null(true_vals)) {
        true_vals <- true_decompositions_mc_nested(pars = pars, mc_n = 1e+05)
    }
    truth_top <- true_vals$topological$components[c("delta_Z", "delta_X", "delta_A", "delta_res", "total")]
    truth_mod <- true_vals$modified$components[c("delta_X", "delta_Z", "delta_A", "delta_res", "total")]
    top_store <- matrix(NA_real_, nrow = B, ncol = length(truth_top))
    colnames(top_store) <- names(truth_top)
    mod_dir_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
    colnames(mod_dir_store) <- names(truth_mod)
    mod_bayes_store <- NULL
    if (fit_X_given_Z) {
        mod_bayes_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
        colnames(mod_bayes_store) <- names(truth_mod)
    }
    conv_Z <- logical(B)
    conv_A <- logical(B)
    conv_Y <- logical(B)
    simulation_seeds <- seed + seq_len(B)
    model_diagnostics <- vector("list", B)
    for (b in seq_len(B)) {
        set.seed(simulation_seeds[b])
        dat_b <- simulate_one_nested(n = n, pars = pars)
        est_b <- estimate_all_decomps_nested(
            dat = dat_b,
            pars = pars,
            fit_X_given_Z = fit_X_given_Z,
            outcome_interaction = outcome_interaction
        )
        top_store[b, ] <- est_b$topological_hat[colnames(top_store)]
        mod_dir_store[b, ] <- est_b$modified_hat_direct[colnames(mod_dir_store)]
        if (fit_X_given_Z) {
            mod_bayes_store[b, ] <- est_b$modified_hat_bayes[colnames(mod_bayes_store)]
        }
        conv_Z[b] <- inherits(est_b$fits$fit_Z_given_X, "multinom")
        conv_A[b] <- inherits(est_b$fits$fit_A_given_XZ, "multinom")
        conv_Y[b] <- inherits(est_b$fits$fit_Y_given_XZA, "glm")
        model_diagnostics[[b]] <- model_diagnostics_nested(
            dat = dat_b,
            fits = est_b$fits,
            pars = pars,
            replicate = b,
            outcome_interaction = outcome_interaction,
            simulation_seed = simulation_seeds[b]
        )
        if (verbose && (b%%max(1, floor(B/10)) == 0)) {
            message(
                "Completed nested estimator replicate ",
                b,
                " / ",
                B,
                " (outcome_interaction = ",
                outcome_interaction,
                ")"
            )
        }
    }
    summ_top <- summarize_estimator_performance(top_store, truth_top)
    summ_mod_dir <- summarize_estimator_performance(mod_dir_store, truth_mod)
    summ_mod_bayes <- NULL
    if (fit_X_given_Z) {
        summ_mod_bayes <- summarize_estimator_performance(mod_bayes_store, truth_mod)
    }
    list(
        params = c(pars, list(n = n, outcome_interaction = outcome_interaction)),
        truth = list(topological = truth_top, modified = truth_mod),
        estimates = list(topological = top_store, modified_direct = mod_dir_store, modified_bayes = mod_bayes_store),
        summary = list(topological = summ_top, modified_direct = summ_mod_dir, modified_bayes = summ_mod_bayes),
        convergence = data.frame(conv_Z = conv_Z, conv_A = conv_A, conv_Y = conv_Y),
        simulation_seeds = data.frame(replicate = seq_len(B), simulation_seed = simulation_seeds),
        model_diagnostics = do.call(rbind, model_diagnostics),
        nested_metadata = nested_calibration_metadata(pars, n = n, outcome_interaction = outcome_interaction)
    )
}

run_estimator_simulation_nested_firth <- function(B = 200,
                                                  n = 3000,
                                                  pars,
                                                  true_vals = NULL,
                                                  fit_X_given_Z = TRUE,
                                                  seed = 1,
                                                  verbose = TRUE) {
    set.seed(seed)
    if (is.null(true_vals)) {
        true_vals <- true_decompositions_mc_nested(pars = pars, mc_n = 1e+05)
    }
    truth_top <- true_vals$topological$components[c("delta_Z", "delta_X", "delta_A", "delta_res", "total")]
    truth_mod <- true_vals$modified$components[c("delta_X", "delta_Z", "delta_A", "delta_res", "total")]
    top_store <- matrix(NA_real_, nrow = B, ncol = length(truth_top))
    colnames(top_store) <- names(truth_top)
    mod_dir_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
    colnames(mod_dir_store) <- names(truth_mod)
    mod_bayes_store <- NULL
    if (fit_X_given_Z) {
        mod_bayes_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
        colnames(mod_bayes_store) <- names(truth_mod)
    }
    conv_Z <- logical(B)
    conv_A <- logical(B)
    conv_Y <- logical(B)
    simulation_seeds <- seed + seq_len(B)
    model_diagnostics <- vector("list", B)
    for (b in seq_len(B)) {
        set.seed(simulation_seeds[b])
        dat_b <- simulate_one_nested(n = n, pars = pars)
        est_b <- estimate_all_decomps_nested_firth(
            dat = dat_b,
            pars = pars,
            fit_X_given_Z = fit_X_given_Z
        )
        top_store[b, ] <- est_b$topological_hat[colnames(top_store)]
        mod_dir_store[b, ] <- est_b$modified_hat_direct[colnames(mod_dir_store)]
        if (fit_X_given_Z) {
            mod_bayes_store[b, ] <- est_b$modified_hat_bayes[colnames(mod_bayes_store)]
        }
        conv_Z[b] <- inherits(est_b$fits$fit_Z_given_X, "multinom")
        conv_A[b] <- inherits(est_b$fits$fit_A_given_XZ, "multinom")
        conv_Y[b] <- inherits(est_b$fits$fit_Y_given_XZA, "logistf")
        model_diagnostics[[b]] <- model_diagnostics_nested(
            dat = dat_b,
            fits = est_b$fits,
            pars = pars,
            replicate = b,
            outcome_interaction = "ZA_firth",
            simulation_seed = simulation_seeds[b]
        )
        if (verbose && (b%%max(1, floor(B/10)) == 0)) {
            message("Completed nested Firth estimator replicate ", b, " / ", B)
        }
    }
    summ_top <- summarize_estimator_performance(top_store, truth_top)
    summ_mod_dir <- summarize_estimator_performance(mod_dir_store, truth_mod)
    summ_mod_bayes <- NULL
    if (fit_X_given_Z) {
        summ_mod_bayes <- summarize_estimator_performance(mod_bayes_store, truth_mod)
    }
    list(
        params = c(pars, list(n = n, outcome_interaction = "ZA_firth")),
        truth = list(topological = truth_top, modified = truth_mod),
        estimates = list(topological = top_store, modified_direct = mod_dir_store, modified_bayes = mod_bayes_store),
        summary = list(topological = summ_top, modified_direct = summ_mod_dir, modified_bayes = summ_mod_bayes),
        convergence = data.frame(conv_Z = conv_Z, conv_A = conv_A, conv_Y = conv_Y),
        simulation_seeds = data.frame(replicate = seq_len(B), simulation_seed = simulation_seeds),
        model_diagnostics = do.call(rbind, model_diagnostics),
        nested_metadata = nested_calibration_metadata(pars, n = n, outcome_interaction = "ZA_firth")
    )
}

print_estimator_simulation_summary <- function(sim_res) {
    cat("\n========================================\n")
    cat("Estimator simulation summary\n")
    cat("========================================\n")
    cat("\nConvergence rates:\n")
    print(round(colMeans(sim_res$convergence), 3))
    cat("\nTopological-order estimator:\n")
    print(round(sim_res$summary$topological, 5))
    cat("\nModified-order estimator (direct Z|X model):\n")
    print(round(sim_res$summary$modified_direct, 5))
    if (!is.null(sim_res$summary$modified_bayes)) {
        cat("\nModified-order estimator (Bayes inversion from X|Z model):\n")
        print(round(sim_res$summary$modified_bayes, 5))
    }
    cat("========================================\n")
}

default_xgb_grid <- function() {
    expand.grid(eta = c(0.05, 0.1), max_depth = c(2, 3, 4), min_child_weight = c(1, 3, 5), stringsAsFactors = FALSE)
}

default_xgb_grid_expanded <- function() {
    expand.grid(
        eta = c(0.03, 0.05, 0.1, 0.2),
        max_depth = c(2, 3, 4, 5, 6),
        min_child_weight = c(0.5, 1, 2, 3),
        stringsAsFactors = FALSE
    )
}

default_xgb_grid_aggressive <- function() {
    expand.grid(
        eta = c(0.03, 0.05, 0.1, 0.2, 0.3),
        max_depth = c(3, 4, 5, 6, 8, 10),
        min_child_weight = c(0, 0.1, 0.5, 1, 2),
        subsample = c(0.8, 1.0),
        colsample_bytree = c(0.8, 1.0),
        lambda = c(0, 0.01, 0.1, 1),
        stringsAsFactors = FALSE
    )
}

default_xgb_grid_aggressive_small <- function() {
    expand.grid(
        eta = c(0.05, 0.1, 0.2, 0.3),
        max_depth = c(4, 6, 8, 10),
        min_child_weight = c(0, 0.1, 0.5, 1),
        subsample = c(1.0),
        colsample_bytree = c(1.0),
        lambda = c(0, 0.1, 1),
        stringsAsFactors = FALSE
    )
}

make_xgb_rhs_formula <- function(x_names) {
    as.formula(paste("~ A + Z +", paste(x_names, collapse = " + ")))
}

make_xgb_matrix <- function(dat, rhs_formula) {
    mm <- Matrix::sparse.model.matrix(rhs_formula, data = dat)
    if ("(Intercept)" %in% colnames(mm)) {
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    }
    mm
}

xgb_grid_value <- function(xgb_grid, name, row, default) {
    if (name %in% names(xgb_grid)) {
        xgb_grid[[name]][row]
    } else {
        default
    }
}

fit_xgb_outcome_cv <- function(dat, pars, xgb_grid = default_xgb_grid(), nfold = 5, max_nrounds = 500, early_stopping_rounds = 20, seed = 123) {
    set.seed(seed)
    x_names <- paste0("X", seq_len(pars$p))
    rhs_formula <- make_xgb_rhs_formula(x_names)
    dat_xgb <- dat
    dat_xgb$Z <- factor(dat_xgb$Z, levels = pars$race_labels)
    dat_xgb$A <- factor(dat_xgb$A, levels = pars$hospital_labels)
    Xmat <- make_xgb_matrix(dat_xgb, rhs_formula)
    yvec <- dat_xgb$Y
    dtrain <- xgboost::xgb.DMatrix(data = Xmat, label = yvec)
    best_score <- Inf
    best_iter <- NA_integer_
    best_params <- NULL
    cv_store <- vector("list", nrow(xgb_grid))
    for (g in seq_len(nrow(xgb_grid))) {
        pars_g <- list(
            booster = "gbtree",
            objective = "binary:logistic",
            eval_metric = "logloss",
            eta = xgb_grid$eta[g],
            max_depth = xgb_grid$max_depth[g],
            min_child_weight = xgb_grid$min_child_weight[g],
            subsample = xgb_grid_value(xgb_grid, "subsample", g, 0.8),
            colsample_bytree = xgb_grid_value(xgb_grid, "colsample_bytree", g, 0.8),
            lambda = xgb_grid_value(xgb_grid, "lambda", g, 1),
            alpha = xgb_grid_value(xgb_grid, "alpha", g, 0),
            gamma = xgb_grid_value(xgb_grid, "gamma", g, 0),
            nthread = xgb_grid_value(xgb_grid, "nthread", g, 4)
        )
        cv_fit <- xgboost::xgb.cv(params = pars_g, data = dtrain, nrounds = max_nrounds, nfold = nfold, early_stopping_rounds = early_stopping_rounds, verbose = 0)
        iter_g <- NULL
        if (!is.null(cv_fit$early_stop) && !is.null(cv_fit$early_stop$best_iteration)) {
            iter_g <- cv_fit$early_stop$best_iteration
        }
        else if (!is.null(cv_fit$best_iteration)) {
            iter_g <- cv_fit$best_iteration
        }
        else {
            iter_g <- nrow(cv_fit$evaluation_log)
        }
        eval_log <- as.data.frame(cv_fit$evaluation_log)
        test_cols <- grep("^test.*logloss.*mean$|^test_logloss_mean$", names(eval_log), value = TRUE)
        if (length(test_cols) == 0) {
            stop("Could not find test log-loss column in xgb.cv evaluation_log.")
        }
        score_g <- eval_log[[test_cols[1]]][iter_g]
        if (length(score_g) != 1 || !is.finite(score_g)) {
            stop("Failed to extract a valid CV log-loss value from xgb.cv output.")
        }
        cv_store[[g]] <- list(params = pars_g, best_iteration = iter_g, best_score = score_g, evaluation_log = eval_log)
        if (score_g < best_score) {
            best_score <- score_g
            best_iter <- iter_g
            best_params <- pars_g
        }
    }
    final_model <- xgboost::xgb.train(params = best_params, data = dtrain, nrounds = best_iter, verbose = 0)
    list(model = final_model, rhs_formula = rhs_formula, x_names = x_names, best_params = best_params, best_nrounds = best_iter, best_cv_logloss = best_score, cv_store = cv_store)
}

predict_xgb_prob <- function(fit_obj, newdata, pars) {
    nd <- newdata
    nd$Z <- factor(nd$Z, levels = pars$race_labels)
    nd$A <- factor(nd$A, levels = pars$hospital_labels)
    Xnew <- make_xgb_matrix(nd, fit_obj$rhs_formula)
    as.numeric(predict(fit_obj$model, newdata = Xnew))
}

fit_decomp_models_xgb <- function(dat, pars, xgb_grid = default_xgb_grid(), nfold = 5, max_nrounds = 500, early_stopping_rounds = 20, seed = 123) {
    x_names <- paste0("X", seq_len(pars$p))
    dat <- dat
    dat$Z <- relevel(factor(dat$Z, levels = pars$race_labels), ref = pars$baseline_race)
    dat$A <- relevel(factor(dat$A, levels = pars$hospital_labels), ref = pars$baseline_hospital)
    pZ_hat <- prop.table(table(dat$Z))
    pZ_hat <- pZ_hat[pars$race_labels]
    form_Z <- as.formula(paste("Z ~", paste(x_names, collapse = " + ")))
    fit_Z_given_X <- tryCatch(nnet::multinom(form_Z, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    form_A <- as.formula(paste("A ~ Z +", paste(x_names, collapse = " + ")))
    fit_A_given_XZ <- tryCatch(nnet::multinom(form_A, data = dat, trace = FALSE, MaxNWts = 5000), error = function(e) e)
    fit_Y_xgb <- fit_xgb_outcome_cv(dat = dat, pars = pars, xgb_grid = xgb_grid, nfold = nfold, max_nrounds = max_nrounds, early_stopping_rounds = early_stopping_rounds, seed = seed)
    list(x_names = x_names, pZ_hat = pZ_hat, fit_Z_given_X = fit_Z_given_X, fit_A_given_XZ = fit_A_given_XZ, fit_Y_xgb = fit_Y_xgb)
}

build_conditional_surfaces_xgb <- function(Xmat, fits, pars) {
    n <- nrow(Xmat)
    R <- pars$R
    H <- pars$H
    x_names <- fits$x_names
    stopifnot(all(colnames(Xmat) == x_names))
    m_array <- array(NA_real_, dim = c(n, R, H), dimnames = list(NULL, pars$race_labels, pars$hospital_labels))
    pA_array <- array(NA_real_, dim = c(n, R, H), dimnames = list(NULL, pars$race_labels, pars$hospital_labels))
    for (r in seq_len(R)) {
        zlab <- pars$race_labels[r]
        newA <- data.frame(Xmat, Z = factor(zlab, levels = pars$race_labels))
        pA_r <- predict_multinom_prob(fits$fit_A_given_XZ, newdata = newA, levels_out = pars$hospital_labels)
        pA_array[, r, ] <- pA_r
        for (h in seq_len(H)) {
            alab <- pars$hospital_labels[h]
            newY <- data.frame(Xmat, Z = factor(zlab, levels = pars$race_labels), A = factor(alab, levels = pars$hospital_labels))
            m_array[, r, h] <- predict_xgb_prob(fits$fit_Y_xgb, newdata = newY, pars = pars)
        }
    }
    g_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    vA_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    rA_mat <- matrix(0, nrow = n, ncol = R, dimnames = list(NULL, pars$race_labels))
    for (r in seq_len(R)) {
        pA_r <- pA_array[, r, ]
        m_r <- m_array[, r, ]
        g_r <- rowSums(pA_r * m_r)
        v_r <- rowSums(pA_r * (m_r - g_r)^2)
        rr_r <- rowSums(pA_r * (m_r * (1 - m_r)))
        g_mat[, r] <- g_r
        vA_mat[, r] <- v_r
        rA_mat[, r] <- rr_r
    }
    list(m_array = m_array, pA_array = pA_array, g_mat = g_mat, vA_mat = vA_mat, rA_mat = rA_mat)
}

estimate_topological_decomp_xgb <- function(dat, fits, pars) {
    x_names <- fits$x_names
    Xmat <- as.matrix(dat[, x_names, drop = FALSE])
    surf <- build_conditional_surfaces_xgb(Xmat, fits, pars)
    g_mat <- surf$g_mat
    vA_mat <- surf$vA_mat
    rA_mat <- surf$rA_mat
    pZ_hat <- as.numeric(fits$pZ_hat)
    names(pZ_hat) <- pars$race_labels
    Z_obs <- dat$Z
    u_z <- numeric(pars$R)
    delta_X <- 0
    delta_A <- 0
    delta_res <- 0
    for (r in seq_len(pars$R)) {
        zlab <- pars$race_labels[r]
        idx <- which(Z_obs == zlab)
        g_r <- g_mat[idx, zlab]
        u_z[r] <- mean(g_r)
        pz_r <- unname(pZ_hat[zlab])
        delta_X <- delta_X + pz_r * mean((g_r - u_z[r])^2)
        delta_A <- delta_A + pz_r * mean(vA_mat[idx, zlab])
        delta_res <- delta_res + pz_r * mean(rA_mat[idx, zlab])
    }
    mu_u <- sum(pZ_hat * u_z)
    delta_Z <- sum(pZ_hat * (u_z - mu_u)^2)
    c(delta_Z = delta_Z, delta_X = delta_X, delta_A = delta_A, delta_res = delta_res, total = delta_Z + delta_X + delta_A + delta_res)
}

estimate_modified_decomp_direct_xgb <- function(dat, fits, pars) {
    x_names <- fits$x_names
    Xmat <- as.matrix(dat[, x_names, drop = FALSE])
    pZ_given_X <- predict_multinom_prob(fits$fit_Z_given_X, newdata = data.frame(Xmat), levels_out = pars$race_labels)
    surf <- build_conditional_surfaces_xgb(Xmat, fits, pars)
    g_mat <- surf$g_mat
    vA_mat <- surf$vA_mat
    rA_mat <- surf$rA_mat
    h_x <- rowSums(pZ_given_X * g_mat)
    delta_X <- mean((h_x - mean(h_x))^2)
    delta_Z <- mean(rowSums(pZ_given_X * (g_mat - h_x)^2))
    delta_A <- mean(rowSums(pZ_given_X * vA_mat))
    delta_res <- mean(rowSums(pZ_given_X * rA_mat))
    c(delta_X = delta_X, delta_Z = delta_Z, delta_A = delta_A, delta_res = delta_res, total = delta_X + delta_Z + delta_A + delta_res)
}

estimate_all_decomps_xgb <- function(dat, pars, xgb_grid = default_xgb_grid(), nfold = 5, max_nrounds = 500, early_stopping_rounds = 20, seed = 123) {
    fits <- fit_decomp_models_xgb(dat = dat, pars = pars, xgb_grid = xgb_grid, nfold = nfold, max_nrounds = max_nrounds, early_stopping_rounds = early_stopping_rounds, seed = seed)
    top_hat <- estimate_topological_decomp_xgb(dat, fits, pars)
    mod_hat <- estimate_modified_decomp_direct_xgb(dat, fits, pars)
    list(fits = fits, topological_hat_xgb = top_hat, modified_hat_direct_xgb = mod_hat)
}

run_estimator_simulation_with_xgb <- function(B = 100,
                                             n = 1000,
                                             pars,
                                             true_vals = NULL,
                                             fit_X_given_Z = TRUE,
                                             xgb_grid = default_xgb_grid(),
                                             nfold = 5,
                                             max_nrounds = 500,
                                             early_stopping_rounds = 20,
                                             seed = 1,
                                             verbose = TRUE) {
    set.seed(seed)
    if (is.null(true_vals)) {
        true_vals <- true_decompositions_mc(pars = pars, mc_n = 1e+05)
    }
    truth_top <- true_vals$topological$components[c("delta_Z", "delta_X", "delta_A", "delta_res", "total")]
    truth_mod <- true_vals$modified$components[c("delta_X", "delta_Z", "delta_A", "delta_res", "total")]
    top_store <- matrix(NA_real_, nrow = B, ncol = length(truth_top))
    colnames(top_store) <- names(truth_top)
    mod_dir_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
    colnames(mod_dir_store) <- names(truth_mod)
    mod_bayes_store <- NULL
    if (fit_X_given_Z) {
        mod_bayes_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
        colnames(mod_bayes_store) <- names(truth_mod)
    }
    top_xgb_store <- matrix(NA_real_, nrow = B, ncol = length(truth_top))
    colnames(top_xgb_store) <- names(truth_top)
    mod_xgb_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
    colnames(mod_xgb_store) <- names(truth_mod)
    conv_Z <- logical(B)
    conv_A <- logical(B)
    conv_Y <- logical(B)
    best_nrounds_xgb <- numeric(B)
    best_logloss_xgb <- numeric(B)
    for (b in seq_len(B)) {
        dat_b <- simulate_one(n = n, pars = pars)
        est_b <- estimate_all_decomps(dat = dat_b, pars = pars, fit_X_given_Z = fit_X_given_Z)
        top_store[b, ] <- est_b$topological_hat[colnames(top_store)]
        mod_dir_store[b, ] <- est_b$modified_hat_direct[colnames(mod_dir_store)]
        if (fit_X_given_Z) {
            mod_bayes_store[b, ] <- est_b$modified_hat_bayes[colnames(mod_bayes_store)]
        }
        conv_Z[b] <- inherits(est_b$fits$fit_Z_given_X, "multinom")
        conv_A[b] <- inherits(est_b$fits$fit_A_given_XZ, "multinom")
        conv_Y[b] <- inherits(est_b$fits$fit_Y_given_XZA, "glm")
        est_xgb_b <- estimate_all_decomps_xgb(dat = dat_b, pars = pars, xgb_grid = xgb_grid, nfold = nfold, max_nrounds = max_nrounds, early_stopping_rounds = early_stopping_rounds, seed = seed + b)
        top_xgb_store[b, ] <- est_xgb_b$topological_hat_xgb[colnames(top_xgb_store)]
        mod_xgb_store[b, ] <- est_xgb_b$modified_hat_direct_xgb[colnames(mod_xgb_store)]
        best_nrounds_xgb[b] <- est_xgb_b$fits$fit_Y_xgb$best_nrounds
        best_logloss_xgb[b] <- est_xgb_b$fits$fit_Y_xgb$best_cv_logloss
        if (verbose && (b%%max(1, floor(B/10)) == 0)) {
            message("Completed full estimator replicate ", b, " / ", B)
        }
    }
    summ_top <- summarize_estimator_performance(top_store, truth_top)
    summ_mod_dir <- summarize_estimator_performance(mod_dir_store, truth_mod)
    summ_mod_bayes <- NULL
    if (fit_X_given_Z) {
        summ_mod_bayes <- summarize_estimator_performance(mod_bayes_store, truth_mod)
    }
    summ_top_xgb <- summarize_estimator_performance(top_xgb_store, truth_top)
    summ_mod_xgb <- summarize_estimator_performance(mod_xgb_store, truth_mod)
    list(params = c(pars, list(n = n)), truth = list(topological = truth_top, modified = truth_mod), estimates = list(topological = top_store, modified_direct = mod_dir_store, modified_bayes = mod_bayes_store, topological_xgb = top_xgb_store, modified_direct_xgb = mod_xgb_store), summary = list(topological = summ_top, modified_direct = summ_mod_dir, modified_bayes = summ_mod_bayes, topological_xgb = summ_top_xgb, modified_direct_xgb = summ_mod_xgb), convergence = data.frame(conv_Z = conv_Z, conv_A = conv_A, conv_Y = conv_Y), 
        xgb_tuning = data.frame(best_nrounds = best_nrounds_xgb, best_cv_logloss = best_logloss_xgb))
}

estimate_all_decomps_xgb_nested <- function(dat,
                                            pars,
                                            xgb_grid = default_xgb_grid(),
                                            nfold = 5,
                                            max_nrounds = 500,
                                            early_stopping_rounds = 20,
                                            seed = 123) {
    estimate_all_decomps_xgb(
        dat = dat,
        pars = pars,
        xgb_grid = xgb_grid,
        nfold = nfold,
        max_nrounds = max_nrounds,
        early_stopping_rounds = early_stopping_rounds,
        seed = seed
    )
}

run_estimator_simulation_with_xgb_nested <- function(B = 100,
                                                    n = 1000,
                                                    pars,
                                                    true_vals = NULL,
                                                    xgb_grid = default_xgb_grid(),
                                                    nfold = 5,
                                                    max_nrounds = 500,
                                                    early_stopping_rounds = 20,
                                                    seed = 1,
                                                    verbose = TRUE,
                                                    xgb_grid_name = "default") {
    set.seed(seed)
    if (is.null(true_vals)) {
        true_vals <- true_decompositions_mc_nested(pars = pars, mc_n = 1e+05)
    }
    truth_top <- true_vals$topological$components[c("delta_Z", "delta_X", "delta_A", "delta_res", "total")]
    truth_mod <- true_vals$modified$components[c("delta_X", "delta_Z", "delta_A", "delta_res", "total")]
    top_xgb_store <- matrix(NA_real_, nrow = B, ncol = length(truth_top))
    colnames(top_xgb_store) <- names(truth_top)
    mod_xgb_store <- matrix(NA_real_, nrow = B, ncol = length(truth_mod))
    colnames(mod_xgb_store) <- names(truth_mod)
    conv_Z <- logical(B)
    conv_A <- logical(B)
    best_nrounds_xgb <- numeric(B)
    best_logloss_xgb <- numeric(B)
    simulation_seeds <- seed + seq_len(B)
    for (b in seq_len(B)) {
        set.seed(simulation_seeds[b])
        dat_b <- simulate_one_nested(n = n, pars = pars)
        est_xgb_b <- estimate_all_decomps_xgb_nested(
            dat = dat_b,
            pars = pars,
            xgb_grid = xgb_grid,
            nfold = nfold,
            max_nrounds = max_nrounds,
            early_stopping_rounds = early_stopping_rounds,
            seed = seed + b
        )
        top_xgb_store[b, ] <- est_xgb_b$topological_hat_xgb[colnames(top_xgb_store)]
        mod_xgb_store[b, ] <- est_xgb_b$modified_hat_direct_xgb[colnames(mod_xgb_store)]
        conv_Z[b] <- inherits(est_xgb_b$fits$fit_Z_given_X, "multinom")
        conv_A[b] <- inherits(est_xgb_b$fits$fit_A_given_XZ, "multinom")
        best_nrounds_xgb[b] <- est_xgb_b$fits$fit_Y_xgb$best_nrounds
        best_logloss_xgb[b] <- est_xgb_b$fits$fit_Y_xgb$best_cv_logloss
        if (verbose && (b%%max(1, floor(B/10)) == 0)) {
            message("Completed nested XGBoost replicate ", b, " / ", B)
        }
    }
    summ_top_xgb <- summarize_estimator_performance(top_xgb_store, truth_top)
    summ_mod_xgb <- summarize_estimator_performance(mod_xgb_store, truth_mod)
    list(
        params = c(pars, list(n = n, outcome_interaction = "xgb", xgb_grid_name = xgb_grid_name)),
        truth = list(topological = truth_top, modified = truth_mod),
        estimates = list(topological_xgb = top_xgb_store, modified_direct_xgb = mod_xgb_store),
        summary = list(topological_xgb = summ_top_xgb, modified_direct_xgb = summ_mod_xgb),
        convergence = data.frame(conv_Z = conv_Z, conv_A = conv_A),
        simulation_seeds = data.frame(replicate = seq_len(B), simulation_seed = simulation_seeds),
        xgb_tuning = data.frame(best_nrounds = best_nrounds_xgb, best_cv_logloss = best_logloss_xgb),
        nested_metadata = nested_calibration_metadata(pars, n = n, outcome_interaction = "xgb", xgb_grid_name = xgb_grid_name)
    )
}

print_estimator_simulation_summary_with_xgb <- function(sim_res) {
    cat("\n========================================\n")
    cat("Estimator simulation summary\n")
    cat("========================================\n")
    cat("\nConvergence rates (parametric Z|X, A|X,Z, Y|X,Z,A):\n")
    print(round(colMeans(sim_res$convergence), 3))
    cat("\nXGBoost tuning summary:\n")
    print(round(c(mean_best_nrounds = mean(sim_res$xgb_tuning$best_nrounds), mean_best_cv_logloss = mean(sim_res$xgb_tuning$best_cv_logloss)), 4))
    cat("\nTopological-order estimator (parametric):\n")
    print(round(sim_res$summary$topological, 5))
    cat("\nModified-order estimator (direct Z|X, parametric):\n")
    print(round(sim_res$summary$modified_direct, 5))
    if (!is.null(sim_res$summary$modified_bayes)) {
        cat("\nModified-order estimator (Bayes inversion from X|Z model):\n")
        print(round(sim_res$summary$modified_bayes, 5))
    }
    cat("\nTopological-order estimator (XGBoost outcome model):\n")
    print(round(sim_res$summary$topological_xgb, 5))
    cat("\nModified-order estimator (direct Z|X + XGBoost outcome model):\n")
    print(round(sim_res$summary$modified_direct_xgb, 5))
    cat("========================================\n")
}

make_scenario_id <- function(prefix = "sim", n, p, H = NULL, extra = NULL) {
    parts <- c(prefix, paste0("n", n), paste0("p", p))
    if (!is.null(H)) {
        parts <- c(parts, paste0("H", H))
    }
    if (!is.null(extra)) {
        parts <- c(parts, extra)
    }
    paste(parts, collapse = "_")
}

format_gamma_ZA_for_id <- function(gamma_ZA) {
    paste0("gammaZA", gsub("-", "m", gsub("[.]", "p", format(gamma_ZA, scientific = FALSE, trim = TRUE))))
}

nested_calibration_metadata <- function(pars,
                                        n = NULL,
                                        outcome_interaction = NULL,
                                        xgb_grid_name = NULL) {
    list(
        n = n,
        p = pars$p,
        H = pars$H,
        R = pars$R,
        gamma_ZA = pars$gamma_ZA,
        outcome_interaction = outcome_interaction,
        xgb_grid_name = xgb_grid_name,
        lambda_Z_calibrated = pars$lambda_Z_calibrated,
        calibration_target_EY_by_Z = pars$calibration_target_EY_by_Z,
        calibration_achieved_EY_by_Z = pars$calibration_achieved_EY_by_Z,
        calibration_mc_n = pars$calibration_mc_n,
        calibration_seed = pars$calibration_seed
    )
}

check_nested_calibration <- function(pars, tol = 0.002) {
    baseline_lambda <- pars$s_ZY * pars$race_direct_scores
    names(baseline_lambda) <- pars$race_labels
    if (isTRUE(all.equal(pars$gamma_ZA, 0))) {
        ok <- isTRUE(all.equal(pars$lambda_Z_calibrated, baseline_lambda, tolerance = tol))
        if (!ok) {
            stop("For gamma_ZA = 0, lambda_Z_calibrated differs from the baseline direct group term.")
        }
        return(invisible(TRUE))
    }
    diff <- max(abs(pars$calibration_achieved_EY_by_Z - pars$calibration_target_EY_by_Z))
    if (!is.finite(diff) || diff >= tol) {
        stop("Nested calibration did not meet tolerance; max absolute difference = ", diff, ".")
    }
    invisible(TRUE)
}

check_decomposition_total <- function(components, tol = 1e-8) {
    component_names <- setdiff(names(components), c("total", "varY_true"))
    if (length(component_names) == 0 || !"total" %in% names(components)) {
        stop("components must include component entries and total.")
    }
    diff <- abs(sum(components[component_names]) - unname(components["total"]))
    if (!is.finite(diff) || diff > tol) {
        stop("Decomposition total check failed; absolute difference = ", diff, ".")
    }
    invisible(TRUE)
}

check_nested_true_decompositions <- function(true_vals, tol = 1e-8) {
    check_decomposition_total(true_vals$topological$components, tol = tol)
    check_decomposition_total(true_vals$modified$components, tol = tol)
    invisible(TRUE)
}

check_nested_simulated_prob <- function(dat) {
    p <- attr(dat, "true_prob_Y")
    if (is.null(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) {
        stop("attr(dat, \"true_prob_Y\") must be finite and in [0, 1].")
    }
    invisible(TRUE)
}

check_nested_gamma0_matches_original <- function(pars_nested,
                                                 mc_n = 1e+05,
                                                 seed_top = 1001,
                                                 seed_mod = 2001,
                                                 tol = 0.003) {
    if (!isTRUE(all.equal(pars_nested$gamma_ZA, 0))) {
        stop("check_nested_gamma0_matches_original() requires gamma_ZA = 0.")
    }
    pars_original <- pars_nested
    pars_original$gamma_ZA <- NULL
    pars_original$lambda_Z_calibrated <- NULL
    pars_original$calibration_target_EY_by_Z <- NULL
    pars_original$calibration_achieved_EY_by_Z <- NULL
    pars_original$calibration_mc_n <- NULL
    pars_original$calibration_seed <- NULL
    pars_original$nested_outcome <- NULL
    original <- true_decompositions_mc(pars_original, mc_n = mc_n, seed_top = seed_top, seed_mod = seed_mod)
    nested <- true_decompositions_mc_nested(pars_nested, mc_n = mc_n, seed_top = seed_top, seed_mod = seed_mod)
    original_top <- original$topological$components[names(nested$topological$components)]
    original_mod <- original$modified$components[names(nested$modified$components)]
    diff_top <- max(abs(nested$topological$components - original_top), na.rm = TRUE)
    diff_mod <- max(abs(nested$modified$components - original_mod), na.rm = TRUE)
    if (diff_top > tol || diff_mod > tol) {
        stop(
            "Nested gamma_ZA = 0 true values differ from original beyond tolerance; ",
            "topological max diff = ",
            diff_top,
            ", modified max diff = ",
            diff_mod,
            "."
        )
    }
    invisible(list(original = original, nested = nested, max_abs_diff_topological = diff_top, max_abs_diff_modified = diff_mod))
}

save_simulation_result <- function(sim_res, pars, scenario_id, out_dir = "outputs/simulation_results", save_raw = TRUE, save_csv = TRUE) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    obj_to_save <- list(scenario_id = scenario_id, scenario = list(n = if (!is.null(sim_res$params$n)) sim_res$params$n else NA, p = pars$p, H = pars$H, R = pars$R, sparse = pars$sparse), pars = pars, truth = sim_res$truth, summary = sim_res$summary, convergence = sim_res$convergence)
    if (!is.null(sim_res$xgb_tuning)) {
        obj_to_save$xgb_tuning <- sim_res$xgb_tuning
    }
    if (save_raw) {
        obj_to_save$estimates <- sim_res$estimates
    }
    rds_file <- file.path(out_dir, paste0(scenario_id, ".rds"))
    saveRDS(obj_to_save, file = rds_file)
    saved_files <- list(rds = rds_file)
    if (save_csv) {
        if (!dir.exists(file.path(out_dir, "csv"))) {
            dir.create(file.path(out_dir, "csv"), recursive = TRUE)
        }
        write_summary_csv <- function(mat, name_stub) {
            df <- data.frame(metric = rownames(mat), mat, row.names = NULL, check.names = FALSE)
            f <- file.path(out_dir, "csv", paste0(scenario_id, "_", name_stub, ".csv"))
            write.csv(df, file = f, row.names = FALSE)
            f
        }
        csv_files <- list()
        if (!is.null(sim_res$summary$topological)) {
            csv_files$topological <- write_summary_csv(sim_res$summary$topological, "topological")
        }
        if (!is.null(sim_res$summary$modified_direct)) {
            csv_files$modified_direct <- write_summary_csv(sim_res$summary$modified_direct, "modified_direct")
        }
        if (!is.null(sim_res$summary$modified_bayes)) {
            csv_files$modified_bayes <- write_summary_csv(sim_res$summary$modified_bayes, "modified_bayes")
        }
        if (!is.null(sim_res$summary$topological_xgb)) {
            csv_files$topological_xgb <- write_summary_csv(sim_res$summary$topological_xgb, "topological_xgb")
        }
        if (!is.null(sim_res$summary$modified_direct_xgb)) {
            csv_files$modified_direct_xgb <- write_summary_csv(sim_res$summary$modified_direct_xgb, "modified_direct_xgb")
        }
        conv_file <- file.path(out_dir, "csv", paste0(scenario_id, "_convergence.csv"))
        write.csv(sim_res$convergence, file = conv_file, row.names = FALSE)
        csv_files$convergence <- conv_file
        if (!is.null(sim_res$xgb_tuning)) {
            xgb_file <- file.path(out_dir, "csv", paste0(scenario_id, "_xgb_tuning.csv"))
            write.csv(sim_res$xgb_tuning, file = xgb_file, row.names = FALSE)
            csv_files$xgb_tuning <- xgb_file
        }
        saved_files$csv <- csv_files
    }
    invisible(saved_files)
}

save_simulation_result_nested <- function(sim_res,
                                          pars,
                                          scenario_id,
                                          out_dir = "outputs/simulation_results_nested_ZA_interaction",
                                          save_raw = TRUE,
                                          save_csv = TRUE) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    metadata <- if (!is.null(sim_res$nested_metadata)) {
        sim_res$nested_metadata
    } else {
        nested_calibration_metadata(pars, n = if (!is.null(sim_res$params$n)) sim_res$params$n else NA)
    }
    obj_to_save <- list(
        scenario_id = scenario_id,
        scenario = list(
            n = if (!is.null(sim_res$params$n)) sim_res$params$n else NA,
            p = pars$p,
            H = pars$H,
            R = pars$R,
            sparse = pars$sparse,
            gamma_ZA = pars$gamma_ZA,
            outcome_interaction = metadata$outcome_interaction,
            xgb_grid_name = metadata$xgb_grid_name
        ),
        pars = pars,
        nested_metadata = metadata,
        truth = sim_res$truth,
        summary = sim_res$summary,
        convergence = sim_res$convergence
    )
    if (!is.null(sim_res$xgb_tuning)) {
        obj_to_save$xgb_tuning <- sim_res$xgb_tuning
    }
    if (!is.null(sim_res$model_diagnostics)) {
        obj_to_save$model_diagnostics <- sim_res$model_diagnostics
    }
    if (!is.null(sim_res$simulation_seeds)) {
        obj_to_save$simulation_seeds <- sim_res$simulation_seeds
    }
    if (save_raw) {
        obj_to_save$estimates <- sim_res$estimates
    }
    rds_file <- file.path(out_dir, paste0(scenario_id, ".rds"))
    saveRDS(obj_to_save, file = rds_file)
    saved_files <- list(rds = rds_file)
    if (save_csv) {
        csv_dir <- file.path(out_dir, "csv")
        if (!dir.exists(csv_dir)) {
            dir.create(csv_dir, recursive = TRUE)
        }
        write_summary_csv <- function(mat, name_stub) {
            df <- data.frame(metric = rownames(mat), mat, row.names = NULL, check.names = FALSE)
            f <- file.path(csv_dir, paste0(scenario_id, "_", name_stub, ".csv"))
            write.csv(df, file = f, row.names = FALSE)
            f
        }
        csv_files <- list()
        for (name_stub in names(sim_res$summary)) {
            if (!is.null(sim_res$summary[[name_stub]])) {
                csv_files[[name_stub]] <- write_summary_csv(sim_res$summary[[name_stub]], name_stub)
            }
        }
        conv_file <- file.path(csv_dir, paste0(scenario_id, "_convergence.csv"))
        write.csv(sim_res$convergence, file = conv_file, row.names = FALSE)
        csv_files$convergence <- conv_file
        if (!is.null(sim_res$xgb_tuning)) {
            xgb_file <- file.path(csv_dir, paste0(scenario_id, "_xgb_tuning.csv"))
            write.csv(sim_res$xgb_tuning, file = xgb_file, row.names = FALSE)
            csv_files$xgb_tuning <- xgb_file
        }
        if (!is.null(sim_res$model_diagnostics)) {
            diagnostics_file <- file.path(csv_dir, paste0(scenario_id, "_model_diagnostics.csv"))
            write.csv(sim_res$model_diagnostics, file = diagnostics_file, row.names = FALSE)
            csv_files$model_diagnostics <- diagnostics_file
        }
        if (!is.null(sim_res$simulation_seeds)) {
            seeds_file <- file.path(csv_dir, paste0(scenario_id, "_simulation_seeds.csv"))
            write.csv(sim_res$simulation_seeds, file = seeds_file, row.names = FALSE)
            csv_files$simulation_seeds <- seeds_file
        }
        calibration_df <- data.frame(
            race = pars$race_labels,
            gamma_ZA = pars$gamma_ZA,
            lambda_Z_calibrated = as.numeric(pars$lambda_Z_calibrated[pars$race_labels]),
            target_EY_by_Z = as.numeric(pars$calibration_target_EY_by_Z[pars$race_labels]),
            achieved_EY_by_Z = as.numeric(pars$calibration_achieved_EY_by_Z[pars$race_labels]),
            calibration_mc_n = pars$calibration_mc_n,
            calibration_seed = pars$calibration_seed,
            stringsAsFactors = FALSE
        )
        calibration_file <- file.path(csv_dir, paste0(scenario_id, "_calibration.csv"))
        write.csv(calibration_df, file = calibration_file, row.names = FALSE)
        csv_files$calibration <- calibration_file
        truth_df <- do.call(rbind, lapply(names(sim_res$truth), function(decomp) {
            data.frame(
                decomposition = decomp,
                component = names(sim_res$truth[[decomp]]),
                value = as.numeric(sim_res$truth[[decomp]]),
                stringsAsFactors = FALSE
            )
        }))
        truth_file <- file.path(csv_dir, paste0(scenario_id, "_truth.csv"))
        write.csv(truth_df, file = truth_file, row.names = FALSE)
        csv_files$truth <- truth_file
        saved_files$csv <- csv_files
    }
    invisible(saved_files)
}

load_simulation_result <- function(scenario_id, out_dir = "outputs/simulation_results") {
    rds_file <- file.path(out_dir, paste0(scenario_id, ".rds"))
    readRDS(rds_file)
}

extract_summary_tidy <- function(saved_obj, estimator_name) {
    mat <- saved_obj$summary[[estimator_name]]
    if (is.null(mat)) {
        return(NULL)
    }
    df <- as.data.frame(mat, check.names = FALSE)
    df$metric <- rownames(df)
    rownames(df) <- NULL
    out <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
        data.frame(scenario_id = saved_obj$scenario_id, estimator = estimator_name, metric = df$metric[i], component = setdiff(names(df), "metric"), value = as.numeric(df[i, setdiff(names(df), "metric")]), n = saved_obj$scenario$n, p = saved_obj$scenario$p, H = saved_obj$scenario$H, stringsAsFactors = FALSE)
    }))
    out
}

