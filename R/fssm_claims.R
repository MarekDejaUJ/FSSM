# ==============================================================================
# FSSM Package - Causal Claims Calculations
# ==============================================================================
# These functions implement the fuzzy set-theoretic causal inference framework.
# All calculations are performed AFTER the cIPMA baseline (importance, performance,
# NCA bottlenecks) has been established.
#
# Four causal claims are computed in both crisp and fuzzy formulations:
# 1. Typical Sufficiency (alpha): P(Y|X) - "X is sufficient for Y"
# 2. Typical Necessity (epsilon): Exception rate - "X is necessary for Y"
# 3. Probabilistic Sufficiency (delta): P(Y|X) - P(Y|~X) - causal uplift
# 4. Probabilistic Necessity (beta): P(Y|~X) - counterfactual baseline
# ==============================================================================


#' Calculate Fuzzy Membership Functions
#'
#' @description Computes fuzzy membership values for sufficiency, necessity,
#' and outcome achievement based on latent scores and NCA bottleneck thresholds.
#'
#' @param X_star Numeric vector. Rescaled latent scores for condition X (0-100).
#' @param Y_star Numeric vector. Rescaled latent scores for outcome Y (0-100).
#' @param T_X Numeric. Bottleneck threshold for X at target Y level (from NCA).
#' @param T_Y Numeric. Target outcome level (e.g., 85).
#'
#' @return A list containing:
#'   \item{mu_S}{Sufficiency membership: degree to which X exceeds threshold}
#'   \item{mu_not_S}{Complement of sufficiency: 1 - mu_S}
#'   \item{mu_N}{Necessity membership: degree to which X meets case-specific bottleneck}
#'   \item{mu_Y}{Outcome membership: degree to which Y achieves target}
#'   \item{necessity_deficit}{How far below bottleneck (normalized)}
#'   \item{bottleneck_values}{Case-specific bottleneck values b(Y)}
#'   \item{T_X, T_Y}{Threshold values}
#'   \item{bottleneck_slope}{Slope of NCA ceiling line}
#'
#' @details
#' ## Membership Function Formulas
#'
#' **Sufficiency Membership (mu_S):**
#' \deqn{\mu_S(x) = \frac{x - T_X}{100 - T_X}}
#' Measures how much X exceeds the bottleneck threshold T_X.
#' Range: [0,1] where 0 = at/below threshold, 1 = at maximum.
#'
#' **Necessity Membership (mu_N):**
#' \deqn{\mu_N(x,y) = \frac{x}{b(y)} = \frac{x}{(T_X/T_Y) \cdot y}}
#' Measures how well X meets the case-specific bottleneck for its own Y value.
#' Range: [0,1] where 0 = far below ceiling, 1 = at/above ceiling.
#'
#' **Outcome Membership (mu_Y):**
#' \deqn{\mu_Y(y) = \min(1, y/T_Y)}
#' Measures achievement toward target outcome level.
#' Range: [0,1] where 0 = no achievement, 1 = target reached/exceeded.
#'
#' @export
calculate_membership <- function(X_star, Y_star, T_X, T_Y) {

  # ============================================================================
  # Input Validation
  # ============================================================================

  n <- length(X_star)

  if (length(Y_star) != n) {
    stop("X_star and Y_star must have the same length.")
  }

  if (is.na(T_X)) {
    warning("T_X is NA. Using T_Y as fallback threshold.")
    T_X <- T_Y
  }

  if (T_Y <= 0) {
    stop("T_Y must be positive.")
  }

  # ============================================================================
  # Constants
  # ============================================================================

  U_X <- 100 # Upper bound of X scale
  L_X <- 0 # Lower bound of X scale

  # NCA ceiling line slope: b(Y) = slope * Y
  # At Y = T_Y, bottleneck = T_X, so slope = T_X / T_Y
  bottleneck_slope <- T_X / T_Y

  # Case-specific bottleneck values
  b_X_Y <- bottleneck_slope * Y_star

  # ============================================================================
  # Sufficiency Membership (mu_S)
  # ============================================================================
  # How much does X exceed the threshold T_X?
  # mu_S = (X - T_X) / (U_X - T_X)
  # Interpretation: 0 = at/below threshold, 1 = at maximum

  denom_S <- U_X - T_X

  if (denom_S > 0) {
    mu_S <- (X_star - T_X) / denom_S
    mu_S <- pmax(0, pmin(1, mu_S))
  } else {
    # Edge case: threshold at maximum (degenerate)
    mu_S <- ifelse(X_star >= T_X, 1, 0)
  }

  mu_not_S <- 1 - mu_S

  # ============================================================================
  # Necessity Membership (mu_N)
  # ============================================================================
  # How well does X meet the bottleneck for ITS OWN outcome level?
  # mu_N = (X - L_X) / (b(Y) - L_X) = X / b(Y) when L_X = 0
  # Interpretation: 0 = far below ceiling, 1 = at/above ceiling

  # b(Y) - L_X = bottleneck_slope * Y - 0 = bottleneck_slope * Y
  denom_N <- b_X_Y - L_X

  # Handle cases where denominator is zero or negative
  # This happens when Y is very low (near zero)
  denom_N_safe <- denom_N
  denom_N_safe[denom_N <= 0] <- 1e-10

  mu_N <- (X_star - L_X) / denom_N_safe
  mu_N <- pmax(0, pmin(1, mu_N))

  # Cases where b(Y) <= L_X: bottleneck is at or below minimum
  # These cases trivially satisfy necessity (mu_N = 1)
  mu_N[b_X_Y <= L_X] <- 1

  # ============================================================================
  # Outcome Membership (mu_Y)
  # ============================================================================
  # How close is Y to the target level T_Y?
  # mu_Y = min(1, Y / T_Y)
  # Interpretation: 0 = no outcome, 1 = target achieved

  mu_Y <- Y_star / T_Y
  mu_Y <- pmax(0, pmin(1, mu_Y))

  # ============================================================================
  # Necessity Deficit
  # ============================================================================
  # How far below the bottleneck is X? (normalized by scale range)
  # Used for violation mass calculation

  necessity_deficit <- pmax(0, b_X_Y - X_star) / U_X

  # ============================================================================
  # Return Results
  # ============================================================================

  list(
    mu_S = mu_S,
    mu_not_S = mu_not_S,
    mu_N = mu_N,
    mu_Y = mu_Y,
    necessity_deficit = necessity_deficit,
    bottleneck_values = b_X_Y,
    T_X = T_X,
    T_Y = T_Y,
    bottleneck_slope = bottleneck_slope,
    n = n
  )
}


#' Calculate All Four Causal Claims
#'
#' @description Master function that computes crisp claims, fuzzy claims,
#' and dose-response estimates for a single predictor.
#'
#' @param membership List from calculate_membership().
#' @param X_star Numeric vector. Latent scores for X (0-100).
#' @param Y_star Numeric vector. Latent scores for Y (0-100).
#' @param T_X Numeric. Bottleneck threshold for X.
#' @param T_Y Numeric. Target outcome level.
#' @param weights Numeric vector. Case weights (will be normalized).
#'
#' @return A list containing:
#'   \item{crisp}{Crisp (binary) formulation of all four claims}
#'   \item{fuzzy}{Fuzzy formulation of all four claims plus additional metrics}
#'   \item{dose_response}{Dose-response curve estimates}
#'
#' @export
calculate_all_claims <- function(membership, X_star, Y_star, T_X, T_Y, weights) {

  # Ensure weights are normalized (sum to 1)
  weights <- weights / sum(weights, na.rm = TRUE)

  # Handle any NA weights
  if (any(is.na(weights))) {
    warning("NA values in weights. Replacing with mean weight.")
    mean_w <- mean(weights, na.rm = TRUE)
    weights[is.na(weights)] <- mean_w
    weights <- weights / sum(weights)
  }

  # Calculate all claim types
  crisp <- calculate_crisp_claims(X_star, Y_star, T_X, T_Y, weights, membership$bottleneck_slope)
  fuzzy <- calculate_fuzzy_claims(membership, weights)
  dose_response <- estimate_dose_response(membership$mu_S, membership$mu_Y, weights)

  list(
    crisp = crisp,
    fuzzy = fuzzy,
    dose_response = dose_response
  )
}


#' Calculate Crisp (Binary) Causal Claims
#'
#' @description Computes the four causal claims using binary (0/1) indicators
#' rather than fuzzy membership. This is the traditional QCA-style formulation.
#'
#' @keywords internal
calculate_crisp_claims <- function(X_star, Y_star, T_X, T_Y, weights, bottleneck_slope) {

  # ============================================================================
  # Binary Indicators
  # ============================================================================

  high_X <- X_star >= T_X # X meets/exceeds threshold
  high_Y <- Y_star >= T_Y # Y meets/exceeds target
  low_X <- !high_X # X below threshold

  # Case-specific bottleneck check for necessity
  b_X_Y <- bottleneck_slope * Y_star
  below_bottleneck <- X_star < b_X_Y # X fails to meet ceiling

  # ============================================================================
  # Typical Sufficiency (alpha)
  # ============================================================================
  # P(Y >= T_Y | X >= T_X)
  # "Among cases with high X, what proportion achieve high Y?"

  num_alpha <- sum(weights[high_X & high_Y], na.rm = TRUE)
  den_alpha <- sum(weights[high_X], na.rm = TRUE)
  alpha <- if (den_alpha > 0) num_alpha / den_alpha else NA_real_

  # ============================================================================
  # Typical Necessity (epsilon = exception rate)
  # ============================================================================
  # P(X < b(Y) | Y >= T_Y)
  # "Among cases with high Y, what proportion fall below the ceiling?"
  # Lower epsilon = stronger necessity

  num_epsilon <- sum(weights[high_Y & below_bottleneck], na.rm = TRUE)
  den_epsilon <- sum(weights[high_Y], na.rm = TRUE)
  epsilon <- if (den_epsilon > 0) num_epsilon / den_epsilon else NA_real_

  # ============================================================================
  # Probabilistic Necessity (beta)
  # ============================================================================
  # P(Y >= T_Y | X < T_X)
  # "Among cases with low X, what proportion achieve high Y?"
  # This is the counterfactual baseline

  num_beta <- sum(weights[low_X & high_Y], na.rm = TRUE)
  den_beta <- sum(weights[low_X], na.rm = TRUE)
  beta <- if (den_beta > 0) num_beta / den_beta else NA_real_

  # ============================================================================
  # Probabilistic Sufficiency (delta = causal uplift)
  # ============================================================================
  # alpha - beta = P(Y|X) - P(Y|~X)
  # "How much does having high X increase probability of high Y?"

  delta <- if (!is.na(alpha) && !is.na(beta)) alpha - beta else NA_real_

  # ============================================================================
  # Counts for Diagnostics
  # ============================================================================

  list(
    alpha = alpha,
    epsilon = epsilon,
    delta = delta,
    beta = beta,
    n_high_X = sum(high_X, na.rm = TRUE),
    n_low_X = sum(low_X, na.rm = TRUE),
    n_high_Y = sum(high_Y, na.rm = TRUE),
    n_low_Y = sum(!high_Y, na.rm = TRUE),
    n_exceptions = sum(high_Y & below_bottleneck, na.rm = TRUE),
    n_high_X_high_Y = sum(high_X & high_Y, na.rm = TRUE),
    n_low_X_high_Y = sum(low_X & high_Y, na.rm = TRUE)
  )
}


#' Calculate Fuzzy Causal Claims
#'
#' @description Computes the four causal claims using fuzzy membership values.
#' This provides a more nuanced assessment than binary crisp claims.
#'
#' @details
#' ## Fuzzy Claim Formulas
#'
#' **Fuzzy Typical Sufficiency (alpha):**
#' \deqn{\alpha = \frac{\sum w_i \cdot \mu_S(x_i) \cdot \mu_Y(y_i)}{\sum w_i \cdot \mu_S(x_i)}}
#'
#' **Fuzzy Typical Necessity (eta, epsilon = 1 - eta):**
#' \deqn{\eta = \frac{\sum w_i \cdot \mu_Y(y_i) \cdot \mu_N(x_i,y_i)}{\sum w_i \cdot \mu_Y(y_i)}}
#'
#' **Fuzzy Probabilistic Necessity (beta):**
#' \deqn{\beta = \frac{\sum w_i \cdot (1-\mu_S(x_i)) \cdot \mu_Y(y_i)}{\sum w_i \cdot (1-\mu_S(x_i))}}
#'
#' **Fuzzy Probabilistic Sufficiency (delta):**
#' \deqn{\delta = \alpha - \beta}
#'
#' @keywords internal
calculate_fuzzy_claims <- function(membership, weights) {

  # Extract membership values
  mu_S <- membership$mu_S
  mu_not_S <- membership$mu_not_S
  mu_N <- membership$mu_N
  mu_Y <- membership$mu_Y
  necessity_deficit <- membership$necessity_deficit

  # Small epsilon to prevent division by zero
  eps <- 1e-10

  # ============================================================================
  # Fuzzy Typical Sufficiency (alpha)
  # ============================================================================
  # Weighted average of mu_Y among cases weighted by mu_S
  # "On average, how much do sufficient cases achieve the outcome?"

  num_alpha <- sum(weights * mu_S * mu_Y, na.rm = TRUE)
  den_alpha <- sum(weights * mu_S, na.rm = TRUE)
  alpha <- if (den_alpha > eps) num_alpha / den_alpha else NA_real_

  # ============================================================================
  # Fuzzy Typical Necessity (eta) and Exception Rate (epsilon)
  # ============================================================================
  # eta: Weighted average of mu_N among cases weighted by mu_Y
  # epsilon = 1 - eta: Exception rate (cases achieving Y without meeting necessity)

  num_eta <- sum(weights * mu_Y * mu_N, na.rm = TRUE)
  den_eta <- sum(weights * mu_Y, na.rm = TRUE)
  eta <- if (den_eta > eps) num_eta / den_eta else NA_real_
  epsilon <- if (!is.na(eta)) 1 - eta else NA_real_

  # ============================================================================
  # Fuzzy Probabilistic Necessity (beta)
  # ============================================================================
  # Weighted average of mu_Y among cases weighted by (1 - mu_S)
  # "On average, how much do insufficient cases achieve the outcome?"

  num_beta <- sum(weights * mu_not_S * mu_Y, na.rm = TRUE)
  den_beta <- sum(weights * mu_not_S, na.rm = TRUE)
  beta <- if (den_beta > eps) num_beta / den_beta else NA_real_

  # ============================================================================
  # Fuzzy Probabilistic Sufficiency (delta = causal uplift)
  # ============================================================================
  # alpha - beta: Difference in outcome achievement between sufficient and insufficient

  delta <- if (!is.na(alpha) && !is.na(beta)) alpha - beta else NA_real_

  # ============================================================================
  # Additional Fuzzy Metrics
  # ============================================================================

  # Violation Mass: Weighted sum of necessity violations
  # High values indicate many cases achieving Y without meeting necessity
  violation_mass <- sum(weights * mu_Y * (1 - mu_N) * necessity_deficit, na.rm = TRUE)

  # Discriminative Power: Ratio of true positives to false positives (fuzzy)
  # High values indicate X strongly discriminates between high and low Y
  num_D <- sum(weights * mu_S * mu_Y, na.rm = TRUE)
  den_D <- sum(weights * mu_S * (1 - mu_Y), na.rm = TRUE) + eps
  discriminative_power <- num_D / den_D

  # Coverage: Proportion of outcome cases that are also sufficient
  # How well does sufficiency "cover" the outcome?
  coverage <- if (den_eta > eps) num_alpha / den_eta else NA_real_

  # Consistency: Same as alpha (proportion of sufficient cases achieving outcome)
  consistency <- alpha

  # ============================================================================
  # Return Results
  # ============================================================================

  list(
    alpha = alpha,
    eta = eta,
    epsilon = epsilon,
    delta = delta,
    beta = beta,
    violation_mass = violation_mass,
    discriminative_power = discriminative_power,
    coverage = coverage,
    consistency = consistency,
    # Summary statistics
    mean_mu_S = mean(mu_S, na.rm = TRUE),
    mean_mu_N = mean(mu_N, na.rm = TRUE),
    mean_mu_Y = mean(mu_Y, na.rm = TRUE)
  )
}


#' Estimate Dose-Response Function
#'
#' @description Estimates the relationship between sufficiency membership (mu_S)
#' and expected outcome membership (mu_Y) using weighted GAM (preferred) or
#' linear regression (fallback).
#'
#' @details
#' The dose-response function models:
#' \deqn{E[\mu_Y | \mu_S = s] = f(s)}
#'
#' Primary method uses GAM with cubic regression splines:
#' \code{gam(mu_Y ~ s(mu_S, bs = "cr", k = 10), weights = w_i)}
#'
#' Fallback to weighted linear regression if GAM fails.
#'
#' Key outputs:
#' - E[mu_Y | mu_S = 0]: Expected outcome for fully insufficient cases
#' - E[mu_Y | mu_S = 1]: Expected outcome for fully sufficient cases
#' - Uplift = E[1] - E[0]: Causal effect of moving from insufficient to sufficient
#'
#' @param mu_S Numeric vector. Sufficiency membership values (0-1).
#' @param mu_Y Numeric vector. Outcome membership values (0-1).
#' @param weights Numeric vector. Case weights (should sum to 1).
#' @param use_gam Logical. Whether to try GAM first. Default TRUE.
#' @param k Integer. Number of knots for GAM spline. Default 10.
#'
#' @return A list containing dose-response estimates and predictions.
#'
#' @keywords internal
estimate_dose_response <- function(mu_S, mu_Y, weights, use_gam = TRUE, k = 10) {

  n <- length(mu_S)

  # Prepare data frame for modeling
  dr_data <- data.frame(
    mu_S = mu_S,
    mu_Y = mu_Y,
    w = weights * n
  )

  # Remove any NA cases
  complete_idx <- complete.cases(dr_data)
  dr_data <- dr_data[complete_idx, ]

  if (nrow(dr_data) < 10) {
    warning("Too few complete cases for dose-response estimation.")
    return(.dose_response_fallback(mu_S, mu_Y, weights))
  }

  # ============================================================================
  # Try GAM first (preferred method per paper)
  # ============================================================================

  gam_model <- NULL
  method_used <- "none"

  if (use_gam && requireNamespace("mgcv", quietly = TRUE)) {

    # Adjust k based on sample size (k cannot exceed unique values)
    n_unique <- length(unique(dr_data$mu_S))
    k_adj <- min(k, n_unique - 1, floor(nrow(dr_data) / 5))
    k_adj <- max(k_adj, 3)

    gam_model <- tryCatch({
      mgcv::gam(
        mu_Y ~ s(mu_S, bs = "cr", k = k_adj),
        weights = w,
        data = dr_data,
        method = "REML"
      )
    }, error = function(e) {
      warning(paste("GAM fitting failed:", e$message, "- falling back to linear."))
      NULL
    })

    if (!is.null(gam_model)) {
      method_used <- "gam"
    }
  }

  # ============================================================================
  # Fallback to weighted linear regression
  # ============================================================================

  lm_model <- NULL

  if (is.null(gam_model)) {
    lm_model <- tryCatch({
      stats::lm(mu_Y ~ mu_S, weights = w, data = dr_data)
    }, error = function(e) {
      warning(paste("Linear regression also failed:", e$message))
      NULL
    })

    if (!is.null(lm_model)) {
      method_used <- "linear"
    }
  }

  # ============================================================================
  # Extract predictions
  # ============================================================================

  eval_points <- c(0, 0.25, 0.5, 0.75, 1.0)
  newdata <- data.frame(mu_S = eval_points)

  if (!is.null(gam_model)) {
    # GAM predictions with standard errors
    pred <- mgcv::predict.gam(gam_model, newdata = newdata, se.fit = TRUE)
    predictions <- pred$fit
    pred_se <- pred$se.fit

    # Model diagnostics
    r_squared <- summary(gam_model)$r.sq
    dev_explained <- summary(gam_model)$dev.expl
    edf <- sum(summary(gam_model)$edf)

    # Linear approximation for backwards compatibility
    beta_0 <- predictions[1]
    beta_1 <- predictions[5] - predictions[1]

  } else if (!is.null(lm_model)) {
    # Linear predictions
    predictions <- stats::predict(lm_model, newdata = newdata)
    pred_se <- tryCatch({
      stats::predict(lm_model, newdata = newdata, se.fit = TRUE)$se.fit
    }, error = function(e) rep(NA_real_, length(eval_points)))

    # Model diagnostics
    r_squared <- summary(lm_model)$r.squared
    dev_explained <- r_squared
    edf <- 2

    beta_0 <- stats::coef(lm_model)[1]
    beta_1 <- stats::coef(lm_model)[2]

  } else {
    return(.dose_response_fallback(mu_S, mu_Y, weights))
  }

  # ============================================================================
  # Compute key estimates
  # ============================================================================

  predictions <- pmax(0, pmin(1, predictions))

  at_0 <- predictions[1]
  at_1 <- predictions[5]
  uplift <- at_1 - at_0

  # Gradient analysis (per paper Step 8b)
  delta_low_med <- predictions[3] - predictions[2]
  delta_med_high <- predictions[4] - predictions[3]

  # ============================================================================
  # Build predictions data frame
  # ============================================================================

  predictions_df <- data.frame(
    mu_S = eval_points,
    E_mu_Y = predictions,
    SE = if (length(pred_se) == length(eval_points)) pred_se else NA_real_,
    stringsAsFactors = FALSE
  )

  if (!all(is.na(predictions_df$SE))) {
    predictions_df$CI_lower <- pmax(0, predictions_df$E_mu_Y - 1.96 * predictions_df$SE)
    predictions_df$CI_upper <- pmin(1, predictions_df$E_mu_Y + 1.96 * predictions_df$SE)
  }

  # ============================================================================
  # Return results
  # ============================================================================

  list(
    E_mu_Y_at_0 = at_0,
    E_mu_Y_at_1 = at_1,
    uplift = uplift,
    intercept = beta_0,
    slope = beta_1,
    method = method_used,
    r_squared = r_squared,
    dev_explained = dev_explained,
    edf = edf,
    delta_low_med = delta_low_med,
    delta_med_high = delta_med_high,
    returns_type = if (delta_med_high > delta_low_med) "increasing" else "diminishing",
    predictions = predictions_df,
    model = if (!is.null(gam_model)) gam_model else lm_model
  )
}


#' Fallback dose-response when models fail
#' @keywords internal
.dose_response_fallback <- function(mu_S, mu_Y, weights) {

  low_S <- mu_S < 0.5
  high_S <- mu_S >= 0.5

  w_low <- sum(weights[low_S], na.rm = TRUE)
  w_high <- sum(weights[high_S], na.rm = TRUE)

  at_0 <- if (w_low > 1e-10) {
    sum(weights[low_S] * mu_Y[low_S], na.rm = TRUE) / w_low
  } else 0.5

  at_1 <- if (w_high > 1e-10) {
    sum(weights[high_S] * mu_Y[high_S], na.rm = TRUE) / w_high
  } else 0.5

  at_0 <- max(0, min(1, at_0))
  at_1 <- max(0, min(1, at_1))

  eval_points <- c(0, 0.25, 0.5, 0.75, 1.0)
  predictions <- at_0 + (at_1 - at_0) * eval_points

  list(
    E_mu_Y_at_0 = at_0,
    E_mu_Y_at_1 = at_1,
    uplift = at_1 - at_0,
    intercept = at_0,
    slope = at_1 - at_0,
    method = "fallback",
    r_squared = NA_real_,
    dev_explained = NA_real_,
    edf = NA_real_,
    delta_low_med = NA_real_,
    delta_med_high = NA_real_,
    returns_type = NA_character_,
    predictions = data.frame(
      mu_S = eval_points,
      E_mu_Y = predictions,
      SE = NA_real_,
      stringsAsFactors = FALSE
    ),
    model = NULL
  )
}


#' Interpret Causal Claims
#'
#' @description Provides human-readable interpretation of causal claim values.
#'
#' @param claims List containing fuzzy or crisp claims.
#' @param thresholds List with alpha, epsilon, delta, beta thresholds.
#'
#' @return Character vector with interpretations.
#' @export
interpret_claims <- function(claims, thresholds = list(alpha = 0.85, epsilon = 0.05,
                                                       delta = 0.15, beta = 0.20)) {

  interpretations <- character(0)

  # Typical Sufficiency
  if (!is.na(claims$alpha)) {
    if (claims$alpha >= thresholds$alpha) {
      interpretations <- c(interpretations,
                           sprintf("Typical Sufficiency SUPPORTED (alpha = %.3f >= %.2f): X is sufficient for Y",
                                   claims$alpha, thresholds$alpha))
    } else {
      interpretations <- c(interpretations,
                           sprintf("Typical Sufficiency NOT supported (alpha = %.3f < %.2f)",
                                   claims$alpha, thresholds$alpha))
    }
  }

  # Typical Necessity
  if (!is.na(claims$epsilon)) {
    if (claims$epsilon <= thresholds$epsilon) {
      interpretations <- c(interpretations,
                           sprintf("Typical Necessity SUPPORTED (epsilon = %.3f <= %.2f): X is necessary for Y",
                                   claims$epsilon, thresholds$epsilon))
    } else {
      interpretations <- c(interpretations,
                           sprintf("Typical Necessity NOT supported (epsilon = %.3f > %.2f)",
                                   claims$epsilon, thresholds$epsilon))
    }
  }

  # Probabilistic Sufficiency
  if (!is.na(claims$delta)) {
    if (claims$delta >= thresholds$delta) {
      interpretations <- c(interpretations,
                           sprintf("Probabilistic Sufficiency SUPPORTED (delta = %.3f >= %.2f): X causally increases Y",
                                   claims$delta, thresholds$delta))
    } else {
      interpretations <- c(interpretations,
                           sprintf("Probabilistic Sufficiency NOT supported (delta = %.3f < %.2f)",
                                   claims$delta, thresholds$delta))
    }
  }

  # Probabilistic Necessity
  if (!is.na(claims$beta)) {
    if (claims$beta <= thresholds$beta) {
      interpretations <- c(interpretations,
                           sprintf("Probabilistic Necessity SUPPORTED (beta = %.3f <= %.2f): ~X prevents Y",
                                   claims$beta, thresholds$beta))
    } else {
      interpretations <- c(interpretations,
                           sprintf("Probabilistic Necessity NOT supported (beta = %.3f > %.2f)",
                                   claims$beta, thresholds$beta))
    }
  }

  interpretations
}
