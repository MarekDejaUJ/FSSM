# ==============================================================================
# FSSM Package - Print and Summary Methods
# ==============================================================================
# Output functions that clearly separate:
# - cIPMA Baseline (Importance, Performance, NCA) - should match cIPMA exactly
# - FSSM Extensions (Fuzzy claims, weighting, causal inference)
# ==============================================================================


#' Print FSSM Results
#'
#' @description Displays a concise overview of FSSM analysis results.
#' Shows both cIPMA baseline metrics and FSSM causal claims.
#'
#' @param x An fssm object.
#' @param digits Integer. Number of decimal places. Default 3.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.fssm <- function(x, digits = 3, ...) {

  cat("\n")
  cat("=======================================================================\n")
  cat("         Fuzzy Subjective Structure Model (FSSM) Analysis             \n")
  cat("=======================================================================\n\n")

  # ============================================================================
  # Model Overview
  # ============================================================================

  cat("MODEL SPECIFICATION\n")
  cat("-------------------\n")
  cat("  Target Construct:    ", x$target_construct, "\n")
  cat("  Target Level:        ", x$thresholds$target_level, "\n")
  cat("  Number of Predictors:", length(x$predictors), "\n")
  cat("  Sample Size:         ", nrow(x$data_rescaled), "\n")
  cat("  Weighting Method:    ", toupper(x$weighting_method), "\n")

  if (length(x$confounders) > 0) {
    cat("  Confounders:         ", paste(x$confounders, collapse = ", "), "\n")
  } else {
    cat("  Confounders:          None (uniform weights)\n")
  }
  cat("\n")

  # ============================================================================
  # Thresholds
  # ============================================================================

  cat("CLAIM THRESHOLDS\n")
  cat("----------------\n")
  cat("  Typical Sufficiency (alpha):    >= ", x$thresholds$alpha, "\n")
  cat("  Typical Necessity (epsilon):    <= ", x$thresholds$epsilon, "\n")
  cat("  Prob. Sufficiency (delta):      >= ", x$thresholds$delta, "\n")
  cat("  Prob. Necessity (beta):         <= ", x$thresholds$beta, "\n\n")

  # ============================================================================
  # cIPMA Baseline Results
  # ============================================================================

  cat("cIPMA BASELINE (Importance-Performance-Necessity)\n")
  cat("-------------------------------------------------\n")

  baseline_cols <- c("Construct", "Importance", "Performance", "NCA_d", "Bottleneck_Level", "Fail_Percentage")
  baseline_cols <- baseline_cols[baseline_cols %in% names(x$results)]

  baseline_df <- x$results[, baseline_cols, drop = FALSE]

  # Round numeric columns
  num_cols <- sapply(baseline_df, is.numeric)
  baseline_df[, num_cols] <- round(baseline_df[, num_cols], digits)

  # Add necessity flag
  if ("NCA_d" %in% names(x$results) && "NCA_p_value" %in% names(x$results)) {
    baseline_df$Necessary <- ifelse(
      x$results$NCA_d >= 0.10 & x$results$NCA_p_value < 0.05,
      "Yes", "No"
    )
  }

  print(baseline_df, row.names = FALSE)
  cat("\n")

  # ============================================================================
  # FSSM Causal Claims (Fuzzy)
  # ============================================================================

  cat("FSSM CAUSAL CLAIMS (Fuzzy Formulation)\n")
  cat("--------------------------------------\n")

  claims_cols <- c("Construct", "Fuzzy_Alpha", "Fuzzy_Epsilon", "Fuzzy_Delta", "Fuzzy_Beta")
  claims_cols <- claims_cols[claims_cols %in% names(x$results)]

  claims_df <- x$results[, claims_cols, drop = FALSE]
  names(claims_df) <- c("Construct", "Alpha", "Epsilon", "Delta", "Beta")

  num_cols <- sapply(claims_df, is.numeric)
  claims_df[, num_cols] <- round(claims_df[, num_cols], digits)

  print(claims_df, row.names = FALSE)
  cat("\n")

  # ============================================================================
  # Claims Support Summary
  # ============================================================================

  cat("CLAIMS SUPPORT SUMMARY\n")
  cat("----------------------\n")

  support_cols <- c("Construct", "Typical_Suff_Supported", "Typical_Nec_Supported",
                    "Prob_Suff_Supported", "Prob_Nec_Supported")
  support_cols <- support_cols[support_cols %in% names(x$results)]

  support_df <- x$results[, support_cols, drop = FALSE]
  names(support_df) <- c("Construct", "Typ.Suff", "Typ.Nec", "Prob.Suff", "Prob.Nec")

  # Convert to Yes/No for readability
  for (col in names(support_df)[-1]) {
    support_df[[col]] <- ifelse(support_df[[col]], "YES", "no")
  }

  print(support_df, row.names = FALSE)

  cat("\n")
  cat("-----------------------------------------------------------------------\n")
  cat("Use summary(x) for detailed output | plot(x) for visualization\n")
  cat("Use extract_fssm(x, 'results') to get full results table\n")
  cat("=======================================================================\n\n")

  invisible(x)
}


#' Summary of FSSM Results
#'
#' @description Provides comprehensive summary including cIPMA baseline,
#' crisp vs fuzzy comparison, and dose-response analysis.
#'
#' @param object An fssm object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns a list with summary data frames.
#' @export
summary.fssm <- function(object, ...) {

  x <- object

  cat("\n")
  cat("========================================================================\n")
  cat("                      FSSM DETAILED SUMMARY                            \n")
  cat("========================================================================\n\n")

  cat("Target:", x$target_construct, "| Level:", x$thresholds$target_level,
      "| N:", nrow(x$data_rescaled), "| Weighting:", toupper(x$weighting_method), "\n\n")

  # ============================================================================
  # 1. PLS-SEM Importance (cIPMA-identical)
  # ============================================================================

  cat("1. PLS-SEM IMPORTANCE (cIPMA Baseline)\n")
  cat("--------------------------------------\n")

  pls_cols <- c("Construct", "Importance", "PLS_p_value", "Performance")
  pls_cols <- pls_cols[pls_cols %in% names(x$results)]

  pls_df <- x$results[, pls_cols, drop = FALSE]

  # Round and add significance
  if ("Importance" %in% names(pls_df)) pls_df$Importance <- round(pls_df$Importance, 4)
  if ("PLS_p_value" %in% names(pls_df)) pls_df$PLS_p_value <- round(pls_df$PLS_p_value, 4)
  if ("Performance" %in% names(pls_df)) pls_df$Performance <- round(pls_df$Performance, 2)

  if ("PLS_p_value" %in% names(pls_df)) {
    pls_df$Sig <- ifelse(pls_df$PLS_p_value < 0.001, "***",
                         ifelse(pls_df$PLS_p_value < 0.01, "**",
                                ifelse(pls_df$PLS_p_value < 0.05, "*", "")))
  }

  print(pls_df, row.names = FALSE)
  cat("Significance: *** p<0.001, ** p<0.01, * p<0.05\n\n")

  # ============================================================================
  # 2. NCA Necessity Analysis (cIPMA-identical)
  # ============================================================================

  cat("2. NCA NECESSITY ANALYSIS (cIPMA Baseline)\n")
  cat("------------------------------------------\n")

  nca_cols <- c("Construct", "NCA_d", "NCA_p_value", "Bottleneck_Level",
                "Fail_Percentage", "Fail_Count")
  nca_cols <- nca_cols[nca_cols %in% names(x$results)]

  nca_df <- x$results[, nca_cols, drop = FALSE]

  # Round
  if ("NCA_d" %in% names(nca_df)) nca_df$NCA_d <- round(nca_df$NCA_d, 3)
  if ("NCA_p_value" %in% names(nca_df)) nca_df$NCA_p_value <- round(nca_df$NCA_p_value, 4)
  if ("Bottleneck_Level" %in% names(nca_df)) nca_df$Bottleneck_Level <- round(nca_df$Bottleneck_Level, 2)
  if ("Fail_Percentage" %in% names(nca_df)) nca_df$Fail_Percentage <- round(nca_df$Fail_Percentage, 2)

  # Add necessity flag
  if ("NCA_d" %in% names(x$results) && "NCA_p_value" %in% names(x$results)) {
    nca_df$Necessary <- ifelse(
      x$results$NCA_d >= 0.10 & x$results$NCA_p_value < 0.05,
      "Yes", "No"
    )
  }

  print(nca_df, row.names = FALSE)
  cat("Necessity criterion: d >= 0.10 AND p < 0.05\n\n")

  # ============================================================================
  # 3. Crisp vs Fuzzy Claims Comparison
  # ============================================================================

  cat("3. CRISP vs FUZZY CLAIMS COMPARISON\n")
  cat("-----------------------------------\n")

  for (pred in x$predictors) {
    cat("\n  ", pred, ":\n", sep = "")

    pred_row <- x$results$Construct == pred

    comp <- data.frame(
      Claim = c("Alpha (Typ.Suff)", "Epsilon (Typ.Nec)",
                "Delta (Prob.Suff)", "Beta (Prob.Nec)"),
      Crisp = c(
        x$results$Crisp_Alpha[pred_row],
        x$results$Crisp_Epsilon[pred_row],
        x$results$Crisp_Delta[pred_row],
        x$results$Crisp_Beta[pred_row]
      ),
      Fuzzy = c(
        x$results$Fuzzy_Alpha[pred_row],
        x$results$Fuzzy_Epsilon[pred_row],
        x$results$Fuzzy_Delta[pred_row],
        x$results$Fuzzy_Beta[pred_row]
      ),
      stringsAsFactors = FALSE
    )

    comp$Crisp <- round(comp$Crisp, 3)
    comp$Fuzzy <- round(comp$Fuzzy, 3)
    comp$Diff <- round(comp$Fuzzy - comp$Crisp, 3)

    # Add threshold and support
    thresholds <- c(x$thresholds$alpha, x$thresholds$epsilon,
                    x$thresholds$delta, x$thresholds$beta)
    directions <- c(">=", "<=", ">=", "<=")

    comp$Threshold <- paste(directions, thresholds)
    comp$Supported <- c(
      comp$Fuzzy[1] >= x$thresholds$alpha,
      comp$Fuzzy[2] <= x$thresholds$epsilon,
      comp$Fuzzy[3] >= x$thresholds$delta,
      comp$Fuzzy[4] <= x$thresholds$beta
    )
    comp$Supported <- ifelse(comp$Supported, "YES", "no")

    print(comp, row.names = FALSE)
  }
  cat("\n")

  # ============================================================================
  # 4. Dose-Response Summary
  # ============================================================================

  cat("4. DOSE-RESPONSE SUMMARY\n")
  cat("------------------------\n")

  dr_df <- do.call(rbind, lapply(x$predictors, function(pred) {
    dr <- x$claims[[pred]]$dose_response

    data.frame(
      Construct = pred,
      Intercept = round(dr$intercept, 3),
      Slope = round(dr$slope, 3),
      R2 = round(dr$r_squared, 3),
      `E[Y|S=0]` = round(dr$E_mu_Y_at_0, 3),
      `E[Y|S=1]` = round(dr$E_mu_Y_at_1, 3),
      Uplift = round(dr$uplift, 3),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  print(dr_df, row.names = FALSE)
  cat("\nE[Y|S=0]: Expected outcome without condition\n")
  cat("E[Y|S=1]: Expected outcome with condition\n")
  cat("Uplift: Causal effect of condition\n")

  # ============================================================================
  # 5. Additional Metrics
  # ============================================================================

  cat("\n5. ADDITIONAL FSSM METRICS\n")
  cat("--------------------------\n")

  add_cols <- c("Construct", "Violation_Mass", "Discriminative_Power")
  add_cols <- add_cols[add_cols %in% names(x$results)]

  if (length(add_cols) > 1) {
    add_df <- x$results[, add_cols, drop = FALSE]
    add_df$Violation_Mass <- round(add_df$Violation_Mass, 4)
    add_df$Discriminative_Power <- round(add_df$Discriminative_Power, 2)
    print(add_df, row.names = FALSE)
  }

  cat("\n========================================================================\n\n")

  # Return summary data invisibly
  invisible(list(
    pls = pls_df,
    nca = nca_df,
    dose_response = dr_df
  ))
}


#' Check FSSM/PLS-SEM Assumptions
#'
#' @description Evaluates model assumptions including measurement quality,
#' discriminant validity, multicollinearity, and weighting diagnostics.
#'
#' @param x An fssm object.
#' @param vif_threshold Numeric. VIF threshold for multicollinearity. Default 5.
#' @param htmt_threshold Numeric. HTMT threshold for discriminant validity. Default 0.90.
#'
#' @return An object of class `fssm_assumptions` with diagnostic information.
#' @export
check_assumptions <- function(x, vif_threshold = 5, htmt_threshold = 0.90) {

  if (!inherits(x, "fssm")) {
    stop("Input must be an fssm object.")
  }

  # Get the model (prefer pls_model if available)
  model <- if (!is.null(x$pls_model)) x$pls_model else x$boot_model

  if (is.null(model)) {
    warning("No PLS model available in fssm object. Limited diagnostics.")
    return(structure(
      list(
        positive_weights = list(passed = NA, note = "No model available"),
        reliability = NULL,
        htmt = NULL,
        vif = NULL,
        sample_size = list(n = nrow(x$data_rescaled), passed = NA),
        weighting = x$weight_diagnostics,
        weighting_method = x$weighting_method,
        thresholds = list(vif = vif_threshold, htmt = htmt_threshold)
      ),
      class = "fssm_assumptions"
    ))
  }

  model_smry <- tryCatch(summary(model), error = function(e) NULL)

  # ============================================================================
  # 1. Outer Weights Check
  # ============================================================================

  outer_weights <- model$outer_weights

  if (!is.null(outer_weights)) {
    weights_df <- as.data.frame(as.table(outer_weights))
    colnames(weights_df) <- c("Indicator", "Construct", "Weight")
    weights_df <- weights_df[weights_df$Weight != 0, ]
    negative_weights <- weights_df[weights_df$Weight < 0, ]

    weights_check <- list(
      passed = nrow(negative_weights) == 0,
      negative_indicators = if (nrow(negative_weights) > 0) negative_weights else NULL,
      n_indicators = nrow(weights_df)
    )
  } else {
    weights_check <- list(passed = NA, note = "Outer weights not available")
  }

  # ============================================================================
  # 2. Reliability
  # ============================================================================

  reliability <- tryCatch({
    if (!is.null(model_smry)) model_smry$reliability else NULL
  }, error = function(e) NULL)

  # ============================================================================
  # 3. HTMT (Discriminant Validity)
  # ============================================================================

  htmt <- tryCatch({
    if (!is.null(model_smry) && !is.null(model_smry$validity)) {
      model_smry$validity$htmt
    } else NULL
  }, error = function(e) NULL)

  htmt_check <- NULL
  if (!is.null(htmt)) {
    htmt_numeric <- htmt
    htmt_numeric[htmt_numeric == "."] <- NA
    htmt_numeric <- suppressWarnings(apply(htmt_numeric, 2, as.numeric))
    if (is.matrix(htmt_numeric)) {
      rownames(htmt_numeric) <- rownames(htmt)
      htmt_values <- htmt_numeric[lower.tri(htmt_numeric)]
      htmt_exceeded <- any(htmt_values > htmt_threshold, na.rm = TRUE)

      # Find problematic pairs
      problematic_pairs <- character(0)
      if (htmt_exceeded && nrow(htmt_numeric) > 1) {
        for (i in 2:nrow(htmt_numeric)) {
          for (j in 1:(i-1)) {
            val <- htmt_numeric[i, j]
            if (!is.na(val) && val > htmt_threshold) {
              problematic_pairs <- c(problematic_pairs,
                                     paste0(rownames(htmt_numeric)[i], " <-> ", colnames(htmt_numeric)[j],
                                            " (", round(val, 3), ")"))
            }
          }
        }
      }

      htmt_check <- list(
        passed = !htmt_exceeded,
        max_value = max(htmt_values, na.rm = TRUE),
        problematic_pairs = if (length(problematic_pairs) > 0) problematic_pairs else NULL
      )
    }
  }

  # ============================================================================
  # 4. VIF (Multicollinearity)
  # ============================================================================

  vif <- tryCatch({
    if (!is.null(model_smry)) model_smry$vif_antecedents else NULL
  }, error = function(e) NULL)

  vif_check <- NULL
  if (!is.null(vif)) {
    all_vifs <- unlist(vif)
    vif_exceeded <- any(all_vifs > vif_threshold, na.rm = TRUE)
    problematic_vif <- names(all_vifs)[!is.na(all_vifs) & all_vifs > vif_threshold]

    vif_check <- list(
      passed = !vif_exceeded,
      max_value = max(all_vifs, na.rm = TRUE),
      problematic = if (length(problematic_vif) > 0) problematic_vif else NULL
    )
  }

  # ============================================================================
  # 5. Sample Size
  # ============================================================================

  sample_size <- nrow(model$data)
  n_constructs <- if (!is.null(model$construct_scores)) ncol(model$construct_scores) else length(x$predictors) + 1
  min_recommended <- max(50, n_constructs * 10)

  sample_check <- list(
    n = sample_size,
    n_constructs = n_constructs,
    min_recommended = min_recommended,
    passed = sample_size >= min_recommended
  )

  # ============================================================================
  # 6. Weighting Diagnostics
  # ============================================================================

  weighting_diag <- x$weight_diagnostics

  # ============================================================================
  # Build Output
  # ============================================================================

  output <- list(
    positive_weights = weights_check,
    reliability = reliability,
    htmt = htmt,
    htmt_check = htmt_check,
    vif = vif,
    vif_check = vif_check,
    sample_size = sample_check,
    weighting = weighting_diag,
    weighting_method = x$weighting_method,
    thresholds = list(vif = vif_threshold, htmt = htmt_threshold)
  )

  class(output) <- "fssm_assumptions"
  output
}


#' Print FSSM Assumptions
#'
#' @param x An fssm_assumptions object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.fssm_assumptions <- function(x, ...) {

  cat("\n")
  cat("=======================================================================\n")
  cat("                    FSSM MODEL DIAGNOSTICS                            \n")
  cat("=======================================================================\n\n")

  # ============================================================================
  # 1. Outer Weights
  # ============================================================================

  cat("1. OUTER WEIGHTS\n")
  cat("----------------\n")

  if (is.na(x$positive_weights$passed)) {
    cat("   [INFO] ", x$positive_weights$note, "\n")
  } else if (x$positive_weights$passed) {
    cat("   [PASS] All ", x$positive_weights$n_indicators, " outer weights are positive.\n")
  } else {
    cat("   [FAIL] Negative weights detected:\n")
    print(x$positive_weights$negative_indicators)
    cat("   Consider checking for reverse-coded items.\n")
  }

  # ============================================================================
  # 2. Reliability
  # ============================================================================

  cat("\n2. RELIABILITY (Cronbach's alpha >= 0.70, rhoC >= 0.70, AVE >= 0.50)\n")
  cat("--------------------------------------------------------------------\n")

  if (!is.null(x$reliability)) {
    rel_display <- as.data.frame(round(x$reliability, 3))

    # Check thresholds
    issues <- character(0)
    if ("alpha" %in% colnames(x$reliability)) {
      low_alpha <- rownames(x$reliability)[x$reliability[, "alpha"] < 0.70]
      if (length(low_alpha) > 0) issues <- c(issues, paste("Low alpha:", paste(low_alpha, collapse = ", ")))
    }
    if ("rhoC" %in% colnames(x$reliability)) {
      low_rhoC <- rownames(x$reliability)[x$reliability[, "rhoC"] < 0.70]
      if (length(low_rhoC) > 0) issues <- c(issues, paste("Low rhoC:", paste(low_rhoC, collapse = ", ")))
    }
    if ("AVE" %in% colnames(x$reliability)) {
      low_AVE <- rownames(x$reliability)[x$reliability[, "AVE"] < 0.50]
      if (length(low_AVE) > 0) issues <- c(issues, paste("Low AVE:", paste(low_AVE, collapse = ", ")))
    }

    print(rel_display)

    if (length(issues) > 0) {
      cat("   [WARN] ", paste(issues, collapse = "; "), "\n")
    } else {
      cat("   [PASS] All reliability thresholds met.\n")
    }
  } else {
    cat("   [INFO] Reliability statistics not available.\n")
  }

  # ============================================================================
  # 3. HTMT
  # ============================================================================

  cat("\n3. DISCRIMINANT VALIDITY (HTMT < ", x$thresholds$htmt, ")\n", sep = "")
  cat("----------------------------------------------\n")

  if (!is.null(x$htmt)) {
    print(x$htmt, quote = FALSE)

    if (!is.null(x$htmt_check)) {
      if (x$htmt_check$passed) {
        cat("   [PASS] All HTMT values below threshold (max: ",
            round(x$htmt_check$max_value, 3), ")\n")
      } else {
        cat("   [FAIL] HTMT exceeded for:\n")
        for (pair in x$htmt_check$problematic_pairs) {
          cat("          ", pair, "\n")
        }
      }
    }
  } else {
    cat("   [INFO] HTMT not available.\n")
  }

  # ============================================================================
  # 4. VIF
  # ============================================================================

  cat("\n4. MULTICOLLINEARITY (VIF < ", x$thresholds$vif, ")\n", sep = "")
  cat("-----------------------------------------\n")

  if (!is.null(x$vif)) {
    for (dv in names(x$vif)) {
      cat("   ", dv, ": ")
      vifs <- x$vif[[dv]]
      cat(paste(names(vifs), "=", round(vifs, 2), collapse = ", "), "\n")
    }

    if (!is.null(x$vif_check)) {
      if (x$vif_check$passed) {
        cat("   [PASS] No multicollinearity issues (max VIF: ",
            round(x$vif_check$max_value, 2), ")\n")
      } else {
        cat("   [FAIL] High VIF for: ", paste(x$vif_check$problematic, collapse = ", "), "\n")
      }
    }
  } else {
    cat("   [INFO] VIF not available.\n")
  }

  # ============================================================================
  # 5. Sample Size
  # ============================================================================

  cat("\n5. SAMPLE SIZE\n")
  cat("--------------\n")
  cat("   N = ", x$sample_size$n, " (constructs: ", x$sample_size$n_constructs,
      ", recommended >= ", x$sample_size$min_recommended, ")\n", sep = "")

  if (x$sample_size$passed) {
    cat("   [PASS] Sample size adequate.\n")
  } else {
    cat("   [WARN] Sample size may be insufficient for stable estimates.\n")
  }

  # ============================================================================
  # 6. Weighting Diagnostics
  # ============================================================================

  cat("\n6. ", toupper(x$weighting_method), " WEIGHTING DIAGNOSTICS\n", sep = "")
  cat("----------------------------------\n")

  if (!is.null(x$weighting) && length(x$weighting) > 0) {
    for (pred in names(x$weighting)) {
      diag <- x$weighting[[pred]]

      if (!is.null(diag$ess)) {
        status <- if (diag$ess_ratio >= 0.8) "[GOOD]" else if (diag$ess_ratio >= 0.5) "[OK]" else "[WARN]"
        cat("   ", pred, ": ESS = ", round(diag$ess, 1),
            " (ratio: ", round(diag$ess_ratio, 2), ") ", status, "\n", sep = "")
      } else if (!is.null(diag$note)) {
        cat("   ", pred, ": ", diag$note, "\n", sep = "")
      }
    }
    cat("\n   ESS Ratio: >= 0.8 good, 0.5-0.8 acceptable, < 0.5 caution\n")
  } else {
    cat("   [INFO] No weighting applied (uniform weights).\n")
  }

  cat("\n=======================================================================\n\n")

  invisible(x)
}


#' Extract Components from FSSM Object
#'
#' @description Convenience function to extract specific components from an fssm object.
#'
#' @param x An fssm object.
#' @param what Character. Component to extract: "results", "claims", "membership",
#'   "weights", "bottlenecks", "dose_response", "data".
#' @param construct Character. Optional specific construct name.
#' @param crisp Logical. For claims, return crisp (TRUE) or fuzzy (FALSE). Default FALSE.
#'
#' @return The requested component.
#'
#' @examples
#' \dontrun{
#' # Get full results table
#' extract_fssm(x, "results")
#'
#' # Get fuzzy claims for all constructs
#' extract_fssm(x, "claims")
#'
#' # Get crisp claims for specific construct
#' extract_fssm(x, "claims", construct = "Adoption_Intention", crisp = TRUE)
#'
#' # Get rescaled data
#' extract_fssm(x, "data")
#' }
#'
#' @export
extract_fssm <- function(x,
                         what = c("results", "claims", "membership", "weights",
                                  "bottlenecks", "dose_response", "data"),
                         construct = NULL,
                         crisp = FALSE) {

  what <- match.arg(what)

  # Full object extractions
  if (what == "results") return(x$results)
  if (what == "data") return(x$data_rescaled)

  if (what == "bottlenecks") {
    return(list(
      levels = x$bottleneck_levels,
      full = x$bottleneck_full
    ))
  }

  # Construct-specific extractions
  if (is.null(construct)) {
    # Return for all constructs
    if (what == "claims") {
      return(lapply(x$claims, function(c) if (crisp) c$crisp else c$fuzzy))
    }
    if (what == "membership") return(x$membership)
    if (what == "weights") return(x$weights)
    if (what == "dose_response") {
      return(lapply(x$claims, function(c) c$dose_response))
    }

  } else {
    # Return for specific construct
    if (!construct %in% x$predictors) {
      stop(paste("Construct", construct, "not found. Available:",
                 paste(x$predictors, collapse = ", ")))
    }

    if (what == "claims") {
      return(if (crisp) x$claims[[construct]]$crisp else x$claims[[construct]]$fuzzy)
    }
    if (what == "membership") return(x$membership[[construct]])
    if (what == "weights") return(x$weights[[construct]])
    if (what == "dose_response") return(x$claims[[construct]]$dose_response)
  }
}


#' Export FSSM Results to Data Frame
#'
#' @description Exports the main results table in a format suitable for publication
#' or further analysis.
#'
#' @param x An fssm object.
#' @param include Character vector. Which sections to include:
#'   "baseline" (cIPMA metrics), "claims" (fuzzy claims), "support" (claim support flags).
#' @param digits Integer. Decimal places for rounding.
#'
#' @return A data frame with selected results.
#' @export
export_results <- function(x,
                           include = c("baseline", "claims", "support"),
                           digits = 3) {

  if (!inherits(x, "fssm")) stop("Input must be an fssm object.")

  include <- match.arg(include, several.ok = TRUE)

  cols <- "Construct"

  if ("baseline" %in% include) {
    cols <- c(cols, "Importance", "Performance", "NCA_d", "NCA_p_value",
              "Bottleneck_Level", "Fail_Percentage")
  }

  if ("claims" %in% include) {
    cols <- c(cols, "Fuzzy_Alpha", "Fuzzy_Epsilon", "Fuzzy_Delta", "Fuzzy_Beta")
  }

  if ("support" %in% include) {
    cols <- c(cols, "Typical_Suff_Supported", "Typical_Nec_Supported",
              "Prob_Suff_Supported", "Prob_Nec_Supported")
  }

  # Keep only columns that exist
  cols <- cols[cols %in% names(x$results)]

  result <- x$results[, cols, drop = FALSE]

  # Round numeric columns
  num_cols <- sapply(result, is.numeric)
  result[, num_cols] <- round(result[, num_cols], digits)

  result
}
