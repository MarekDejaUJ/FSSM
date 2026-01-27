# ==============================================================================
# FSSM Package - Weighting Functions
# ==============================================================================
# These functions compute case weights for causal inference adjustment.
# Weights are applied AFTER the cIPMA baseline (importance, performance, NCA)
# and only affect the fuzzy/crisp causal claims calculations.
#
# Two methods available:
# 1. GPS (Generalized Propensity Score) - for continuous treatment (fuzzy membership)
# 2. Binary IPW (Inverse Probability Weighting) - for dichotomized treatment
# ==============================================================================


#' Calculate GPS Weights for Continuous Treatment
#'
#' @description Estimates Generalized Propensity Score (GPS) weights for continuous
#' fuzzy sufficiency membership. Uses a linear model to estimate the conditional
#' density of the treatment given confounders, then computes stabilized weights.
#'
#' @param mu_S Numeric vector. Fuzzy sufficiency membership values (0-1 scale).
#' @param confounders Data frame. Confounder variables to adjust for.
#' @param trim_quantiles Numeric vector of length 2. Quantiles for weight trimming.
#'   Default c(0.01, 0.99) trims at 1st and 99th percentiles.
#' @param min_complete Integer. Minimum complete cases required. Default 10.
#'
#' @return A list containing:
#'   \item{weights}{Numeric vector of normalized weights (sum to 1)}
#'   \item{diagnostics}{List with ESS, ESS ratio, model info, and balance stats}
#'
#' @details
#' The GPS weighting procedure:
#' 1. Fit linear model: mu_S ~ confounders
#' 2. Compute conditional density: f(mu_S | Z) using residual SD
#' 3. Compute marginal density: f(mu_S) using overall mean/SD
#' 4. Calculate stabilized weight: w = f(mu_S) / f(mu_S | Z)
#' 5. Trim extreme weights and normalize
#'
#' Stabilized weights reduce variance compared to unstabilized weights.
#'
#' @references
#' Hirano, K., & Imbens, G. W. (2004). The propensity score with continuous treatments.
#' Applied Bayesian Modeling and Causal Inference from Incomplete-Data Perspectives.
#'
#' @importFrom stats lm predict sd dnorm quantile residuals coef cor
#' @export
calculate_gps_weights <- function(mu_S,
                                  confounders,
                                  trim_quantiles = c(0.01, 0.99),
                                  min_complete = 10) {


  n <- length(mu_S)


  # ============================================================================
  # Input Validation
  # ============================================================================

  # Check for empty/invalid confounders

  if (is.null(confounders) || ncol(confounders) == 0 || nrow(confounders) != n) {
    return(list(
      weights = rep(1/n, n),
      diagnostics = list(
        method = "uniform",
        note = "No valid confounders provided",
        ess = n,
        ess_ratio = 1.0
      )
    ))
  }

  # Check mu_S validity
  if (any(is.na(mu_S)) || any(mu_S < 0) || any(mu_S > 1)) {
    warning("mu_S contains NA or out-of-range values. These cases will receive mean weight.")
  }

  # ============================================================================
  # Prepare Data & Handle Missingness
  # ============================================================================

  gps_data <- data.frame(mu_S = mu_S, confounders, check.names = FALSE)
  complete_idx <- complete.cases(gps_data)
  n_complete <- sum(complete_idx)
  gps_data_complete <- gps_data[complete_idx, , drop = FALSE]

  if (n_complete < min_complete) {
    warning(paste0(
      "Only ", n_complete, " complete cases (minimum ", min_complete, " required). ",
      "Returning uniform weights."
    ))
    return(list(
      weights = rep(1/n, n),
      diagnostics = list(
        method = "uniform",
        note = paste("Insufficient complete cases:", n_complete),
        ess = n,
        ess_ratio = 1.0
      )
    ))
  }

  # ============================================================================
  # Fit GPS Model (Linear Regression)
  # ============================================================================

  confounder_names <- names(confounders)
  formula_str <- paste("mu_S ~", paste(confounder_names, collapse = " + "))

  gps_model <- tryCatch({
    stats::lm(stats::as.formula(formula_str), data = gps_data_complete)
  }, error = function(e) {
    warning(paste("GPS model fitting failed:", e$message))
    NULL
  })

  if (is.null(gps_model)) {
    return(list(
      weights = rep(1/n, n),
      diagnostics = list(
        method = "uniform",
        note = "GPS model fitting failed",
        ess = n,
        ess_ratio = 1.0
      )
    ))
  }

  # ============================================================================
  # Calculate Densities
  # ============================================================================

  # Conditional density: f(mu_S | Z)
  mu_hat <- stats::predict(gps_model)
  sigma_hat <- stats::sd(stats::residuals(gps_model))

  # Protect against zero/tiny sigma
  if (is.na(sigma_hat) || sigma_hat < 1e-10) {
    sigma_hat <- stats::sd(gps_data_complete$mu_S) / 2
    if (sigma_hat < 1e-10) sigma_hat <- 0.1
  }

  conditional_density <- pmax(
    stats::dnorm(gps_data_complete$mu_S, mean = mu_hat, sd = sigma_hat),
    1e-10
  )

  # Marginal density: f(mu_S) - stabilizer
  marginal_mean <- mean(gps_data_complete$mu_S, na.rm = TRUE)
  marginal_sd <- stats::sd(gps_data_complete$mu_S, na.rm = TRUE)

  if (is.na(marginal_sd) || marginal_sd < 1e-10) {
    marginal_sd <- 0.1
  }

  marginal_density <- pmax(
    stats::dnorm(gps_data_complete$mu_S, mean = marginal_mean, sd = marginal_sd),
    1e-10
  )

  # ============================================================================
  # Calculate Stabilized Weights
  # ============================================================================

  raw_weights <- marginal_density / conditional_density

  # Trim extreme weights
  w_quantiles <- stats::quantile(raw_weights, trim_quantiles, na.rm = TRUE)
  trimmed_weights <- pmax(pmin(raw_weights, w_quantiles[2]), w_quantiles[1])

  # Normalize (sum to 1 for complete cases)
  normalized_weights <- trimmed_weights / sum(trimmed_weights)

  # ============================================================================
  # Reconstruct Full Weight Vector
  # ============================================================================

  full_weights <- rep(NA_real_, n)
  full_weights[complete_idx] <- normalized_weights

  # Impute missing cases with mean weight (neutral contribution)
  if (any(is.na(full_weights))) {
    mean_weight <- mean(normalized_weights, na.rm = TRUE)
    full_weights[is.na(full_weights)] <- mean_weight
  }

  # Final renormalization
  full_weights <- full_weights / sum(full_weights)

  # ============================================================================
  # Diagnostics
  # ============================================================================

  # Effective Sample Size
  ess <- (sum(full_weights)^2) / sum(full_weights^2)
  ess_ratio <- ess / n

  # Model R-squared (how much variation in mu_S explained by confounders)
  r_squared <- tryCatch({
    summary(gps_model)$r.squared
  }, error = function(e) NA_real_)

  # Weight statistics
  weight_stats <- list(
    min = min(full_weights),
    max = max(full_weights),
    cv = stats::sd(full_weights) / mean(full_weights)
  )

  diagnostics <- list(
    method = "gps",
    n_complete = n_complete,
    n_missing = n - n_complete,
    ess = ess,
    ess_ratio = ess_ratio,
    sigma_hat = sigma_hat,
    r_squared = r_squared,
    trim_quantiles = trim_quantiles,
    weight_stats = weight_stats
  )

  list(weights = full_weights, diagnostics = diagnostics)
}


#' Calculate Binary IPW Weights
#'
#' @description Estimates Inverse Probability of Treatment Weights (IPTW) for
#' binary treatment indicator. Uses logistic regression to estimate propensity
#' scores, then computes stabilized weights.
#'
#' @param X_star Numeric vector. Latent scores (0-100 scale).
#' @param T_X Numeric. Threshold for dichotomization (bottleneck level).
#' @param confounders Data frame. Confounder variables to adjust for.
#' @param trim_quantiles Numeric vector of length 2. Quantiles for weight trimming.
#'   Default c(0.01, 0.99).
#' @param ps_bounds Numeric vector of length 2. Bounds for propensity scores to
#'   avoid extreme weights. Default c(0.01, 0.99).
#' @param min_complete Integer. Minimum complete cases required. Default 10.
#'
#' @return A list containing:
#'   \item{weights}{Numeric vector of normalized weights (sum to 1)}
#'   \item{diagnostics}{List with ESS, ESS ratio, model info, and balance stats}
#'
#' @details
#' The binary IPW procedure:
#' 1. Dichotomize treatment: T = 1 if X >= T_X, else 0
#' 2. Fit logistic model: P(T=1 | Z) = logit^{-1}(Z * beta)
#' 3. Calculate propensity score: e(Z) = predicted probability
#' 4. Calculate stabilized weight:
#'    - Treated (T=1): w = P(T=1) / e(Z)
#'    - Control (T=0): w = P(T=0) / (1 - e(Z))
#' 5. Trim extreme weights and normalize
#'
#' @references
#' Austin, P. C. (2011). An introduction to propensity score methods for reducing
#' the effects of confounding in observational studies.
#'
#' @importFrom stats glm binomial predict quantile
#' @export
calculate_binary_weights <- function(X_star,
                                     T_X,
                                     confounders,
                                     trim_quantiles = c(0.01, 0.99),
                                     ps_bounds = c(0.01, 0.99),
                                     min_complete = 10) {

  n <- length(X_star)

  # ============================================================================
  # Create Binary Treatment Indicator
  # ============================================================================

  T_binary <- as.integer(X_star >= T_X)

  # Check for separation issues
  n_treated <- sum(T_binary, na.rm = TRUE)
  n_control <- sum(!T_binary, na.rm = TRUE)

  if (n_treated < 5 || n_control < 5) {
    warning(paste0(
      "Extreme imbalance: ", n_treated, " treated, ", n_control, " control. ",
      "Binary IPW may be unstable. Consider GPS weighting instead."
    ))
  }

  # ============================================================================
  # Input Validation
  # ============================================================================

  if (is.null(confounders) || ncol(confounders) == 0 || nrow(confounders) != n) {
    return(list(
      weights = rep(1/n, n),
      diagnostics = list(
        method = "uniform",
        note = "No valid confounders provided",
        n_treated = n_treated,
        n_control = n_control,
        ess = n,
        ess_ratio = 1.0
      )
    ))
  }

  # ============================================================================
  # Prepare Data & Handle Missingness
  # ============================================================================

  ipw_data <- data.frame(T_bin = T_binary, confounders, check.names = FALSE)
  complete_idx <- complete.cases(ipw_data)
  n_complete <- sum(complete_idx)
  ipw_data_complete <- ipw_data[complete_idx, , drop = FALSE]

  if (n_complete < min_complete) {
    warning(paste0(
      "Only ", n_complete, " complete cases. Returning uniform weights."
    ))
    return(list(
      weights = rep(1/n, n),
      diagnostics = list(
        method = "uniform",
        note = paste("Insufficient complete cases:", n_complete),
        ess = n,
        ess_ratio = 1.0
      )
    ))
  }

  # ============================================================================
  # Fit Propensity Score Model (Logistic Regression)
  # ============================================================================

  confounder_names <- names(confounders)
  formula_str <- paste("T_bin ~", paste(confounder_names, collapse = " + "))

  ps_model <- tryCatch({
    stats::glm(
      stats::as.formula(formula_str),
      data = ipw_data_complete,
      family = stats::binomial(link = "logit")
    )
  }, error = function(e) {
    warning(paste("Propensity score model failed:", e$message))
    NULL
  })

  if (is.null(ps_model)) {
    return(list(
      weights = rep(1/n, n),
      diagnostics = list(
        method = "uniform",
        note = "Propensity score model failed",
        ess = n,
        ess_ratio = 1.0
      )
    ))
  }

  # Check convergence
  if (!ps_model$converged) {
    warning("Propensity score model did not converge. Results may be unreliable.")
  }

  # ============================================================================
  # Calculate Propensity Scores
  # ============================================================================

  propensity_scores <- stats::predict(ps_model, type = "response")

  # Bound propensity scores to avoid extreme weights
  propensity_scores <- pmax(pmin(propensity_scores, ps_bounds[2]), ps_bounds[1])

  # ============================================================================
  # Calculate Stabilized Weights
  # ============================================================================

  # Marginal probability of treatment (stabilizer)
  p_treated <- mean(ipw_data_complete$T_bin, na.rm = TRUE)
  p_control <- 1 - p_treated

  # Stabilized weights
  raw_weights <- ifelse(
    ipw_data_complete$T_bin == 1,
    p_treated / propensity_scores,
    p_control / (1 - propensity_scores)
  )

  # Trim extreme weights
  w_quantiles <- stats::quantile(raw_weights, trim_quantiles, na.rm = TRUE)
  trimmed_weights <- pmax(pmin(raw_weights, w_quantiles[2]), w_quantiles[1])

  # Normalize
  normalized_weights <- trimmed_weights / sum(trimmed_weights)

  # ============================================================================
  # Reconstruct Full Weight Vector
  # ============================================================================

  full_weights <- rep(NA_real_, n)
  full_weights[complete_idx] <- normalized_weights

  # Impute missing with mean weight
  if (any(is.na(full_weights))) {
    mean_weight <- mean(normalized_weights, na.rm = TRUE)
    full_weights[is.na(full_weights)] <- mean_weight
  }

  # Final renormalization
  full_weights <- full_weights / sum(full_weights)

  # ============================================================================
  # Diagnostics
  # ============================================================================

  # Effective Sample Size
  ess <- (sum(full_weights)^2) / sum(full_weights^2)
  ess_ratio <- ess / n

  # Model AUC (if possible)
  auc <- tryCatch({
    # Simple AUC calculation without extra packages
    pred <- propensity_scores
    actual <- ipw_data_complete$T_bin
    n1 <- sum(actual == 1)
    n0 <- sum(actual == 0)
    if (n1 > 0 && n0 > 0) {
      sum(outer(pred[actual == 1], pred[actual == 0], ">")) / (n1 * n0)
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  # Weight statistics
  weight_stats <- list(
    min = min(full_weights),
    max = max(full_weights),
    cv = stats::sd(full_weights) / mean(full_weights)
  )

  diagnostics <- list(
    method = "binary_ipw",
    n_complete = n_complete,
    n_missing = n - n_complete,
    n_treated = sum(ipw_data_complete$T_bin == 1),
    n_control = sum(ipw_data_complete$T_bin == 0),
    p_treated = p_treated,
    ess = ess,
    ess_ratio = ess_ratio,
    model_converged = ps_model$converged,
    auc = auc,
    ps_bounds = ps_bounds,
    trim_quantiles = trim_quantiles,
    weight_stats = weight_stats
  )

  list(weights = full_weights, diagnostics = diagnostics)
}


#' Create Uniform Weights (No Adjustment)
#'
#' @description Returns uniform weights when no confounder adjustment is needed.
#' Useful as explicit "no weighting" option.
#'
#' @param n Integer. Number of observations.
#'
#' @return A list with uniform weights and diagnostics.
#' @export
calculate_uniform_weights <- function(n) {

  list(
    weights = rep(1/n, n),
    diagnostics = list(
      method = "uniform",
      note = "No confounder adjustment applied",
      ess = n,
      ess_ratio = 1.0
    )
  )
}


#' Summarize Weight Diagnostics
#'
#' @description Prints a summary of weighting diagnostics for all predictors.
#'
#' @param weight_diagnostics List of diagnostics from fssm object.
#' @param predictors Character vector of predictor names.
#'
#' @return Invisibly returns a summary data frame.
#' @export
summarize_weights <- function(weight_diagnostics, predictors) {

  summary_df <- data.frame(
    Construct = predictors,
    Method = NA_character_,
    ESS = NA_real_,
    ESS_Ratio = NA_real_,
    Note = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(predictors)) {
    pred <- predictors[i]
    diag <- weight_diagnostics[[pred]]

    if (!is.null(diag)) {
      summary_df$Method[i] <- diag$method %||% "unknown"
      summary_df$ESS[i] <- diag$ess %||% NA_real_
      summary_df$ESS_Ratio[i] <- diag$ess_ratio %||% NA_real_
      summary_df$Note[i] <- diag$note %||% ""
    }
  }

  cat("\n=== Weighting Diagnostics ===\n\n")
  print(summary_df, row.names = FALSE)

  cat("\nInterpretation:\n")
  cat("- ESS Ratio > 0.8: Good (minimal information loss)\
")
  cat("- ESS Ratio 0.5-0.8: Acceptable (moderate adjustment)\n")
  cat("- ESS Ratio < 0.5: Caution (substantial information loss)\n\n")

  invisible(summary_df)
}


# Null coalescing operator (if not already defined)
`%||%` <- function(a, b) if (is.null(a)) b else a
