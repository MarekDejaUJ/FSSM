#' Run Fuzzy Subjective Structure Model (FSSM) Analysis
#'
#' @description The main function to perform FSSM analysis. It accepts either a raw dataset (via fssm_extract)
#' or a PLS-SEM model. It orchestrates the NCA bottleneck calculation, fuzzy membership estimation,
#' weighting (GPS or Binary), and causal claim validation.
#'
#' IMPORTANT: The pre-FSSM baseline (Importance, Performance, NCA_d, Bottleneck levels)
#' is computed IDENTICALLY to cIPMA. FSSM-specific steps (fuzzy membership, GPS/binary weighting,
#' causal claims) are applied on top of this baseline.
#'
#' @param input Either a bootstrapped seminr model OR an object from `fssm_extract()`.
#' @param target_construct String. The name of the dependent variable.
#' @param data Optional. Raw data frame (required if input is a model).
#' @param scales Optional. List of scales (required if input is a model).
#' @param target_level Numeric. The desired level of the target construct (0-100 scale). Default is 85.
#' @param confounders Character vector. Names of confounding variables for GPS weighting.
#' @param nca_rep Integer. Number of permutations for NCA. Default 10000.
#' @param weighting Character. Either "gps" (default) or "binary" for IPW method.
#' @param alpha_threshold Numeric. Threshold for typical sufficiency. Default 0.85.
#' @param epsilon_threshold Numeric. Threshold for typical necessity exception rate. Default 0.05.
#' @param delta_threshold Numeric. Threshold for probabilistic sufficiency uplift. Default 0.15.
#' @param beta_threshold Numeric. Threshold for probabilistic necessity. Default 0.20.
#' @param pls_model Optional. The base (non-bootstrapped) seminr model for full diagnostics.
#'
#' @return An object of class `fssm` containing results, claims, weights, and diagnostics.
#'
#' @importFrom NCA nca_analysis
#' @importFrom stats pnorm sd lm glm binomial dnorm predict residuals coef quantile
#' @export
fssm <- function(input,
                 target_construct,
                 data = NULL,
                 scales = NULL,
                 target_level = 85,
                 confounders = NULL,
                 nca_rep = 10000,
                 weighting = c("gps", "binary"),
                 alpha_threshold = 0.85,
                 epsilon_threshold = 0.05,
                 delta_threshold = 0.15,
                 beta_threshold = 0.20,
                 pls_model = NULL) {

  weighting <- match.arg(weighting)

  # ========================================
  # STEP 0: Data Preparation (Unified Route)
  # ========================================
  # Route A: fssm_object (from fssm_extract - raw data with entropy weighting)
  # Route B: seminr model (use cIPMA-identical rescaling)

  if (inherits(input, "fssm_object")) {
    # --- ROUTE A: fssm_object (Raw Data / Entropy Mode) ---
    lv_df <- input$scores
    raw_data <- input$raw_data
    predictors <- setdiff(names(lv_df), target_construct)

    if (!target_construct %in% names(lv_df)) {
      stop(paste("Target construct", target_construct, "not found in fssm object scores."))
    }

    # No PLS results available for raw data route
    pls_res <- data.frame(
      Construct = predictors,
      Importance = NA_real_,
      Performance = NA_real_,
      PLS_p_value = NA_real_,
      stringsAsFactors = FALSE
    )

    # Calculate performance from entropy-weighted scores
    for (i in seq_along(predictors)) {
      pred_name <- predictors[i]
      if (pred_name %in% colnames(lv_df)) {
        pls_res$Performance[i] <- mean(lv_df[[pred_name]], na.rm = TRUE)
      }
    }

  } else if (inherits(input, "boot_seminr_model") || inherits(input, "seminr_model")) {
    # --- ROUTE B: PLS-SEM Model (cIPMA-identical baseline) ---

    if (is.null(data) || is.null(scales)) {
      stop("Must provide 'data' and 'scales' for SEM input.")
    }

    # Validate bootstrapped model
    if (!inherits(input, "boot_seminr_model") || is.null(input$boot_paths)) {
      stop("The seminr model must be bootstrapped to obtain CIs and p-values.")
    }

    # ========================================
    # STEP 1: Extract PLS-SEM Results (cIPMA-identical)
    # ========================================
    smry <- summary(input)
    boot_total <- smry$bootstrapped_total_paths

    if (is.null(boot_total)) {
      stop("bootstrapped_total_paths not found in summary(input). Ensure 'input' is a bootstrapped seminr model.")
    }

    # Find rows targeting our outcome
    pattern <- paste0("\\s*->\\s*", target_construct, "$")
    target_rows_idx <- grep(pattern, rownames(boot_total))

    if (length(target_rows_idx) == 0) {
      stop(paste0("Target construct '", target_construct, "' not found as an outcome in the model."))
    }

    target_matrix <- boot_total[target_rows_idx, , drop = FALSE]
    predictors <- sub(pattern, "", rownames(target_matrix))

    # ========================================
    # STEP 2: Rescale LV Scores (cIPMA-identical)
    # ========================================
    # This is the EXACT rescaling procedure from cIPMA

    original_indicators <- as.data.frame(data)
    std_weights_matrix <- input$outer_weights

    # Convert to data frame for manipulation
    std_weights_df <- as.data.frame(as.table(std_weights_matrix))
    colnames(std_weights_df) <- c("Indicator", "Latent_variable", "Standardized_weight")
    std_weights_df <- std_weights_df[std_weights_df$Standardized_weight != 0, ]

    # Check for negative weights (warning only)
    negative_weights <- std_weights_df$Standardized_weight < 0
    if (any(negative_weights)) {
      neg_indicators <- std_weights_df$Indicator[negative_weights]
      warning(paste0(
        "Negative outer weights detected for: ",
        paste(neg_indicators, collapse = ", "),
        ". This may invert performance scores. Consider checking for reverse-coded items."
      ))
    }

    # Get standard deviations (use stored if available, else compute)
    if (!is.null(input$sdData)) {
      std_weights_df$SD <- input$sdData[as.character(std_weights_df$Indicator)]
    } else {
      std_weights_df$SD <- apply(
        data[, as.character(std_weights_df$Indicator), drop = FALSE],
        2, stats::sd, na.rm = TRUE
      )
    }

    # Unstandardize weights
    std_weights_df$Unstandardized_weight <- std_weights_df$Standardized_weight / std_weights_df$SD

    # Normalize within each LV
    sum_weights <- tapply(std_weights_df$Unstandardized_weight, std_weights_df$Latent_variable, sum)
    std_weights_df$Sum_weight <- sum_weights[as.character(std_weights_df$Latent_variable)]
    std_weights_df$Norm_Unstd_weight <- std_weights_df$Unstandardized_weight / std_weights_df$Sum_weight

    # Compute unscaled LV scores
    lv_unscaled <- data.frame(row.names = seq_len(nrow(data)))
    unique_lvs <- unique(as.character(std_weights_df$Latent_variable))

    for (lv in unique_lvs) {
      lv_dat <- std_weights_df[std_weights_df$Latent_variable == lv, ]
      inds <- as.character(lv_dat$Indicator)
      wts <- lv_dat$Norm_Unstd_weight
      ind_scores <- as.matrix(original_indicators[, inds, drop = FALSE])
      lv_unscaled[[as.character(lv)]] <- as.vector(ind_scores %*% wts)
    }

    # Rescale to 0-100 using THEORETICAL min/max from scales
    lv_df <- lv_unscaled
    for (lv in names(lv_unscaled)) {
      # Find first indicator for this LV to get scale
      related_ind <- as.character(
        std_weights_df$Indicator[std_weights_df$Latent_variable == lv][1]
      )

      if (related_ind %in% names(scales)) {
        min_val <- scales[[related_ind]][1]
        max_val <- scales[[related_ind]][2]
        norm_val <- ((lv_unscaled[[lv]] - min_val) / (max_val - min_val)) * 100
        norm_val <- pmax(0, pmin(100, norm_val))
        lv_df[[lv]] <- norm_val
      }
    }

    raw_data <- data

    # ========================================
    # STEP 3: Build PLS Results Table (cIPMA-identical)
    # ========================================
    pls_res <- data.frame(
      Construct = predictors,
      Importance = NA_real_,
      Performance = NA_real_,
      PLS_p_value = NA_real_,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(pls_res))) {
      pred_name <- pls_res$Construct[i]

      # Find matching row in target_matrix
      row_key <- rownames(target_matrix)[
        grep(paste0("^", pred_name, "\\s*->"), rownames(target_matrix))
      ][1]

      if (!is.na(row_key)) {
        # Extract importance (total effect)
        pls_res$Importance[i] <- target_matrix[row_key, "Original Est."]

        # Extract p-value
        cols <- colnames(target_matrix)
        if ("p" %in% cols) {
          pls_res$PLS_p_value[i] <- target_matrix[row_key, "p"]
        } else if ("T Stat." %in% cols) {
          pls_res$PLS_p_value[i] <- 2 * (1 - stats::pnorm(abs(target_matrix[row_key, "T Stat."])))
        }
      }

      # Calculate performance as mean of rescaled scores
      if (pred_name %in% colnames(lv_df)) {
        pls_res$Performance[i] <- mean(lv_df[[pred_name]], na.rm = TRUE)
      }
    }

  } else {
    stop("Input must be an 'fssm_object' (from fssm_extract) or a 'seminr_model'.")
  }

  # ========================================
  # STEP 4: NCA Analysis (cIPMA-identical)
  # ========================================
  nca_data <- lv_df[, c(predictors, target_construct), drop = FALSE]

  # Run NCA with permutation test for significance
  nca_res_sig <- NCA::nca_analysis(
    data = nca_data,
    x = predictors,
    y = target_construct,
    ceilings = "ce_fdh",
    corner = 1,
    test.rep = nca_rep,
    steps = 10
  )

  # Extract NCA statistics
  nca_stats <- data.frame(
    Construct = predictors,
    NCA_d = NA_real_,
    NCA_p_value = NA_real_,
    Is_Necessary = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(predictors)) {
    pred <- predictors[i]

    if (!is.null(nca_res_sig$tests[[pred]])) {
      nca_stats$NCA_d[i] <- nca_res_sig$tests[[pred]]$ce_fdh$observed
      nca_stats$NCA_p_value[i] <- nca_res_sig$tests[[pred]]$ce_fdh$p_value
    } else if (!is.null(nca_res_sig$summaries[[pred]])) {
      nca_stats$NCA_d[i] <- nca_res_sig$summaries[[pred]]$params[2]
    }

    d_val <- nca_stats$NCA_d[i]
    p_val <- nca_stats$NCA_p_value[i]
    if (!is.na(d_val) && !is.na(p_val)) {
      nca_stats$Is_Necessary[i] <- (d_val >= 0.10 & p_val < 0.05)
    }
  }

  # ========================================
  # STEP 5: Bottleneck Analysis (cIPMA-identical)
  # ========================================
  # IMPORTANT: Calculate fail count IMMEDIATELY after extracting bottleneck
  # (exactly like cIPMA) to ensure same threshold value is used

  bottleneck_levels <- list()
  bottleneck_fail_counts <- list()
  bottleneck_fail_pcts <- list()
  ceiling_params <- list()
  eps <- .Machine$double.eps^0.5
  n_cases <- nrow(lv_df)

  for (pred in predictors) {
    raw_threshold <- NA_real_

    # Method 1: Try to get ceiling line parameters from NCA summaries
    # NCA stores ceiling line as: Required_X = intercept + slope * Y
    if (!is.null(nca_res_sig$summaries[[pred]])) {
      smry <- nca_res_sig$summaries[[pred]]

      # Try different locations where NCA might store ceiling line params
      line_params <- NULL

      if (!is.null(smry$lines$ce_fdh)) {
        line_params <- smry$lines$ce_fdh
      } else if (!is.null(smry$ce_fdh$line)) {
        line_params <- smry$ce_fdh$line
      }

      if (!is.null(line_params) && length(line_params) >= 2) {
        intercept <- line_params[1]
        slope <- line_params[2]

        if (!is.na(intercept) && !is.na(slope)) {
          raw_threshold <- intercept + slope * target_level
          ceiling_params[[pred]] <- list(intercept = intercept, slope = slope)
        }
      }
    }

    # Method 2: Run dedicated NCA analysis for this predictor and extract from table
    # Use high-precision extraction
    if (is.na(raw_threshold)) {
      model_b <- NCA::nca_analysis(
        data = nca_data,
        x = pred,
        y = target_construct,
        ceilings = "ce_fdh",
        corner = 1,
        bottleneck.y = "actual",
        bottleneck.x = "actual",
        steps = c(target_level, NA)
      )

      # Try to get ceiling line from this dedicated analysis
      if (!is.null(model_b$summaries[[pred]]$lines$ce_fdh)) {
        line_params <- model_b$summaries[[pred]]$lines$ce_fdh
        if (length(line_params) >= 2) {
          intercept <- line_params[1]
          slope <- line_params[2]
          if (!is.na(intercept) && !is.na(slope)) {
            raw_threshold <- intercept + slope * target_level
            ceiling_params[[pred]] <- list(intercept = intercept, slope = slope)
          }
        }
      }

      # Fallback: extract from bottleneck table
      if (is.na(raw_threshold)) {
        full_bn_table <- model_b$bottlenecks$ce_fdh

        if (target_construct %in% colnames(full_bn_table)) {
          y_vals <- as.numeric(as.character(full_bn_table[[target_construct]]))
        } else {
          y_vals <- as.numeric(as.character(full_bn_table[, 1]))
        }

        target_idx <- which(abs(y_vals - target_level) < 0.001)[1]

        if (!is.na(target_idx)) {
          # Try to extract numeric value directly first (preserves precision)
          # Fall back to character conversion only if needed
          bn_val <- full_bn_table[target_idx, pred]
          if (is.numeric(bn_val)) {
            raw_threshold <- bn_val
          } else if (is.factor(bn_val)) {
            raw_threshold <- suppressWarnings(as.numeric(levels(bn_val)[bn_val]))
          } else {
            raw_threshold <- suppressWarnings(as.numeric(as.character(bn_val)))
          }
        }
      }
    }

    bottleneck_levels[[pred]] <- raw_threshold

    # Calculate fail count IMMEDIATELY (cIPMA-identical)
    # This ensures we use the exact same threshold value
    if (is.na(raw_threshold)) {
      count_below <- 0L
    } else {
      count_below <- sum(lv_df[[pred]] < (raw_threshold - eps), na.rm = TRUE)
    }
    bottleneck_fail_counts[[pred]] <- as.integer(count_below)
    bottleneck_fail_pcts[[pred]] <- (count_below / n_cases) * 100
  }

  # Full bottleneck table for output
  nca_full <- NCA::nca_analysis(
    data = nca_data,
    x = predictors,
    y = target_construct,
    ceilings = "ce_fdh",
    corner = 1,
    bottleneck.y = "actual",
    bottleneck.x = "actual",
    steps = 20
  )

  # Package NCA results
  nca_results <- list(
    nca_stats = nca_stats,
    bottleneck_levels = bottleneck_levels,
    bottleneck_fail_counts = bottleneck_fail_counts,
    bottleneck_fail_pcts = bottleneck_fail_pcts,
    bottleneck_full = nca_full$bottlenecks$ce_fdh,
    nca_full = nca_res_sig
  )

  # ========================================
  # STEP 6: FSSM-Specific Calculations

  # (Fuzzy membership, weighting, causal claims)
  # ========================================
  # NOTE: From here on, we apply FSSM-specific methodology
  # The baseline (pls_res, nca_stats, bottleneck_levels) is now cIPMA-identical

  n <- nrow(lv_df)

  # Determine confounders
  # If confounders = NULL, use uniform weights (no GPS adjustment)
  # If confounders = "auto", auto-detect non-indicator columns
  # If confounders = c("Age", "Gender"), use those specific columns
  if (is.null(confounders)) {
    # NULL means no GPS - uniform weights
    confounders <- character(0)
  } else if (length(confounders) == 1 && confounders == "auto") {
    # Auto-detect: use columns that are NOT indicators and NOT constructs
    indicator_cols <- unique(unlist(lapply(names(lv_df), function(lv) {
      if (exists("std_weights_df")) {
        as.character(std_weights_df$Indicator[std_weights_df$Latent_variable == lv])
      } else {
        character(0)
      }
    })))
    confounders <- setdiff(names(raw_data), c(names(lv_df), indicator_cols))
  }
  confounders <- intersect(confounders, names(raw_data))

  confounder_df <- if (length(confounders) > 0) {
    raw_data[, confounders, drop = FALSE]
  } else {
    data.frame(row.names = seq_len(n))
  }

  # Initialize storage
  claims_list <- list()
  membership_list <- list()
  weights_list <- list()
  weight_diagnostics <- list()

  for (pred in predictors) {
    # Get bottleneck threshold (T_X)
    T_X <- bottleneck_levels[[pred]]

    if (is.na(T_X)) {
      warning("NCA bottleneck not found for ", pred, ". Defaulting to target_level.")
      T_X <- target_level
    }

    # Calculate fuzzy membership
    membership <- calculate_membership(
      X_star = lv_df[[pred]],
      Y_star = lv_df[[target_construct]],
      T_X = T_X,
      T_Y = target_level
    )
    membership_list[[pred]] <- membership

    # Calculate weights (GPS or Binary)
    if (length(confounders) > 0 && ncol(confounder_df) > 0) {
      if (weighting == "gps") {
        wt_result <- calculate_gps_weights(membership$mu_S, confounder_df)
      } else {
        wt_result <- calculate_binary_weights(lv_df[[pred]], T_X, confounder_df)
      }
      weights_list[[pred]] <- wt_result$weights
      weight_diagnostics[[pred]] <- wt_result$diagnostics
    } else {
      # Uniform weights if no confounders
      weights_list[[pred]] <- rep(1 / n, n)
      weight_diagnostics[[pred]] <- list(note = "No confounders, uniform weights")
    }

    # Calculate all four causal claims
    claims_list[[pred]] <- calculate_all_claims(
      membership = membership,
      X_star = lv_df[[pred]],
      Y_star = lv_df[[target_construct]],
      T_X = T_X,
      T_Y = target_level,
      weights = weights_list[[pred]]
    )
  }

  # ========================================
  # STEP 7: Compile Results Table
  # ========================================
  results <- compile_fssm_results(
    predictors = predictors,
    pls_res = pls_res,
    nca_results = nca_results,
    claims_list = claims_list,
    lv_df = lv_df,
    target_level = target_level,
    thresholds = list(
      alpha = alpha_threshold,
      epsilon = epsilon_threshold,
      delta = delta_threshold,
      beta = beta_threshold
    )
  )

  # ========================================
  # STEP 8: Build Output Object
  # ========================================
  output <- list(
    results = results,
    claims = claims_list,
    nca_summary = nca_results$nca_stats,
    bottleneck_levels = nca_results$bottleneck_levels,
    bottleneck_full = nca_results$bottleneck_full,
    membership = membership_list,
    weights = weights_list,
    weight_diagnostics = weight_diagnostics,
    weighting_method = weighting,
    thresholds = list(
      alpha = alpha_threshold,
      epsilon = epsilon_threshold,
      delta = delta_threshold,
      beta = beta_threshold,
      target_level = target_level
    ),
    data_rescaled = lv_df,
    confounders = confounders,
    predictors = predictors,
    target_construct = target_construct,
    boot_model = if (inherits(input, "boot_seminr_model")) input else NULL,
    pls_model = pls_model
  )

  class(output) <- "fssm"
  output
}


#' Compile FSSM Results Table
#'
#' @keywords internal
compile_fssm_results <- function(predictors, pls_res, nca_results, claims_list,
                                 lv_df, target_level, thresholds) {

  # Start with PLS results
  results <- pls_res

  # Merge NCA statistics
  results <- merge(results, nca_results$nca_stats, by = "Construct", all.x = TRUE)

  # Add bottleneck levels
  results$Bottleneck_Level <- vapply(results$Construct, function(x) {
    bl <- nca_results$bottleneck_levels[[x]]
    if (is.null(bl)) NA_real_ else bl
  }, numeric(1))

  # Add crisp claims
  results$Crisp_Alpha <- vapply(results$Construct, function(x) {
    claims_list[[x]]$crisp$alpha
  }, numeric(1))

  results$Crisp_Epsilon <- vapply(results$Construct, function(x) {
    claims_list[[x]]$crisp$epsilon
  }, numeric(1))

  results$Crisp_Delta <- vapply(results$Construct, function(x) {
    claims_list[[x]]$crisp$delta
  }, numeric(1))

  results$Crisp_Beta <- vapply(results$Construct, function(x) {
    claims_list[[x]]$crisp$beta
  }, numeric(1))

  # Add fuzzy claims
  results$Fuzzy_Alpha <- vapply(results$Construct, function(x) {
    claims_list[[x]]$fuzzy$alpha
  }, numeric(1))

  results$Fuzzy_Epsilon <- vapply(results$Construct, function(x) {
    claims_list[[x]]$fuzzy$epsilon
  }, numeric(1))

  results$Fuzzy_Delta <- vapply(results$Construct, function(x) {
    claims_list[[x]]$fuzzy$delta
  }, numeric(1))

  results$Fuzzy_Beta <- vapply(results$Construct, function(x) {
    claims_list[[x]]$fuzzy$beta
  }, numeric(1))

  # Additional fuzzy metrics
  results$Violation_Mass <- vapply(results$Construct, function(x) {
    claims_list[[x]]$fuzzy$violation_mass
  }, numeric(1))

  results$Discriminative_Power <- vapply(results$Construct, function(x) {
    claims_list[[x]]$fuzzy$discriminative_power
  }, numeric(1))

  # Determine support for claims
  results$Typical_Suff_Supported <- results$Fuzzy_Alpha >= thresholds$alpha

  is_nec <- results$NCA_d >= 0.10 & results$NCA_p_value < 0.05
  results$Typical_Nec_Supported <- (results$Fuzzy_Epsilon <= thresholds$epsilon) &
    (!is.na(is_nec)) & is_nec

  results$Prob_Suff_Supported <- results$Fuzzy_Delta >= thresholds$delta
  results$Prob_Nec_Supported <- results$Fuzzy_Beta <= thresholds$beta

  # Calculate fail percentage (use pre-calculated values from main loop for precision)
  results$Fail_Percentage <- vapply(results$Construct, function(x) {
    pct <- nca_results$bottleneck_fail_pcts[[x]]
    if (is.null(pct)) NA_real_ else pct
  }, numeric(1))

  results$Fail_Count <- vapply(results$Construct, function(x) {
    cnt <- nca_results$bottleneck_fail_counts[[x]]
    if (is.null(cnt)) NA_integer_ else as.integer(cnt)
  }, integer(1))

  results
}
