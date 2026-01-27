# ==============================================================================
# FSSM Package - Utility Functions
# ==============================================================================
# These utilities support both PLS-SEM and raw data (entropy) workflows.
# For PLS-SEM route, rescaling follows cIPMA methodology exactly.
# ==============================================================================

# ------------------------------------------------------------------------------
# SYNTAX PARSING (for raw data / entropy route)
# ------------------------------------------------------------------------------

#' Internal Syntax Parser for MCDM/Entropy Mode
#'
#' @description Parses lavaan-style syntax to extract construct-indicator mappings.
#'
#' @param syntax Character string with syntax like "IM =~ item1 + item2; PU =~ item3 + item4"
#'
#' @return A named list mapping construct names to indicator vectors.
#' @keywords internal
.parse_mcdm_syntax <- function(syntax) {
  clean_syntax <- gsub("\n", "", syntax)
  lines <- strsplit(clean_syntax, ";")[[1]]

  mapping <- list()
  for (line in lines) {
    if (trimws(line) == "") next

    parts <- strsplit(line, "=~")[[1]]
    if (length(parts) == 2) {
      crit_name <- trimws(parts[1])
      items <- trimws(strsplit(parts[2], "\\+")[[1]])
      mapping[[crit_name]] <- items
    }
  }

  mapping
}


# ------------------------------------------------------------------------------
# RESCALING FUNCTIONS
# ------------------------------------------------------------------------------

#' Rescale Vector to 0-100 Range
#'
#' @description Simple linear rescaling using provided min/max values.
#'
#' @param x Numeric vector to rescale.
#' @param min_val Theoretical minimum value.
#' @param max_val Theoretical maximum value.
#'
#' @return Numeric vector rescaled to 0-100.
#' @keywords internal
rescale_0_100 <- function(x, min_val, max_val) {
  if (is.na(min_val) || is.na(max_val)) return(x)
  if (max_val == min_val) return(rep(50, length(x)))

  ((x - min_val) / (max_val - min_val)) * 100
}


#' Unwrap seminr Bootstrap Objects
#'
#' @description Extracts the base model from a bootstrapped seminr object.
#'
#' @param model A seminr model (possibly bootstrapped).
#'
#' @return The unwrapped base model.
#' @keywords internal
.unwrap_seminr_model <- function(model) {

  if (!is.null(model$original_model)) return(model$original_model)
  if (!is.null(model$model)) return(model$model)
  model
}


#' Get Rescaled Latent Variable Scores (cIPMA-Identical)
#'
#' @description Computes latent variable scores and rescales them to 0-100 using
#' theoretical min/max from the scales parameter. This function replicates the
#' exact methodology used in cIPMA for consistency.
#'
#' @param model A seminr model (base or bootstrapped).
#' @param data The raw indicator data frame.
#' @param scales A named list of theoretical min/max for indicators.
#'   Example: list(PU_01 = c(1, 7), PU_02 = c(1, 7), ...)
#'
#' @return A data frame with rescaled latent variable scores (0-100).
#'
#' @details
#' The rescaling process:
#' 1. Extract standardized outer weights from the model
#' 2. Unstandardize weights by dividing by indicator SD
#' 3. Normalize weights within each construct (sum to 1)
#' 4. Compute weighted sum of original indicators
#' 5. Rescale to 0-100 using THEORETICAL min/max (not observed!)
#'
#' This matches cIPMA exactly and differs from observed min-max rescaling.
#'
#' @export
get_rescaled_lv_scores <- function(model, data, scales) {

  # Unwrap if bootstrapped
  base_model <- .unwrap_seminr_model(model)

  if (is.null(base_model$outer_weights)) {
    stop("Outer weights not found in model.")
  }

  # Prepare data
  original_indicators <- as.data.frame(data)

  # Get standardized outer weights

  std_weights_matrix <- base_model$outer_weights
  std_weights_df <- as.data.frame(as.table(std_weights_matrix))
  colnames(std_weights_df) <- c("Indicator", "Latent_variable", "Standardized_weight")

  # Keep only non-zero weights
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

  # Get standard deviations
  # Use stored SDs if available (from seminr), else compute from data
  if (!is.null(base_model$sdData)) {
    std_weights_df$SD <- base_model$sdData[as.character(std_weights_df$Indicator)]
  } else {
    std_weights_df$SD <- apply(
      data[, as.character(std_weights_df$Indicator), drop = FALSE],
      2, stats::sd, na.rm = TRUE
    )
  }

  # Unstandardize weights: w_unstd = w_std / SD
  std_weights_df$Unstandardized_weight <- std_weights_df$Standardized_weight / std_weights_df$SD

  # Normalize within each LV (sum of unstandardized weights = 1)
  sum_weights <- tapply(
    std_weights_df$Unstandardized_weight,
    std_weights_df$Latent_variable,
    sum
  )
  std_weights_df$Sum_weight <- sum_weights[as.character(std_weights_df$Latent_variable)]
  std_weights_df$Norm_Unstd_weight <- std_weights_df$Unstandardized_weight / std_weights_df$Sum_weight

  # Compute unscaled LV scores as weighted sum of original indicators
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
    # Find first indicator for this LV to get the scale
    related_ind <- as.character(
      std_weights_df$Indicator[std_weights_df$Latent_variable == lv][1]
    )

    if (!is.null(scales) && related_ind %in% names(scales)) {
      min_val <- scales[[related_ind]][1]
      max_val <- scales[[related_ind]][2]

      # Linear rescaling to 0-100
      norm_val <- ((lv_unscaled[[lv]] - min_val) / (max_val - min_val)) * 100

      # Clip to [0, 100] to handle any extrapolation
      norm_val <- pmax(0, pmin(100, norm_val))

      lv_df[[lv]] <- norm_val
    } else {
      # If no scale provided, keep original (with warning)
      warning(paste0(
        "No scale found for construct '", lv,
        "' (first indicator: '", related_ind, "'). ",
        "Scores will not be rescaled to 0-100."
      ))
    }
  }

  lv_df
}


# ------------------------------------------------------------------------------
# ENTROPY WEIGHTING (for raw data route)
# ------------------------------------------------------------------------------

#' Calculate Entropy-Weighted Composite Score
#'
#' @description Computes Shannon entropy weights for indicators and returns
#' a weighted composite score rescaled to 0-100.
#'
#' @param df Data frame with indicator columns.
#'
#' @return Numeric vector of entropy-weighted scores (0-100).
#'
#' @details
#' The entropy weighting procedure:
#' 1. Normalize each column to sum to 1 (proportion matrix P_ij)
# NOTE: .calculate_entropy_score is defined in fssm_extract.R
# (removed duplicate from here to avoid conflicts)


# ------------------------------------------------------------------------------
# SAATY SCALING (for AHP-style analysis if needed)
# ------------------------------------------------------------------------------

#' Scale Vector to Saaty Scale (1-9)
#'
#' @description Rescales a numeric vector to the 1-9 Saaty scale used in AHP.
#'
#' @param vec Numeric vector to scale.
#'
#' @return Numeric vector on 1-9 scale.
#' @keywords internal
.scale_to_saaty <- function(vec) {

  if (any(vec < 0, na.rm = TRUE)) {
    stop("Negative values detected in raw data.")
  }

  # Treat NA and 99 as missing/zero
  vec[is.na(vec) | vec == 99] <- 0

  valid_mask <- vec > 0
  vals <- vec[valid_mask]

  if (length(vals) == 0) return(vec)

  min_v <- min(vals)
  max_v <- max(vals)

  if (min_v == max_v) {
    vec[valid_mask] <- 1
  } else {
    # Linear rescaling to 1-9
    vec[valid_mask] <- 1 + (vals - min_v) * (8 / (max_v - min_v))
  }

  vec
}


# ------------------------------------------------------------------------------
# DEPRECATED FUNCTIONS (kept for backward compatibility)
# ------------------------------------------------------------------------------

#' @title Get Rescaled Scores (Internal Fallback) - DEPRECATED
#'
#' @description This function uses OBSERVED min-max rescaling which differs from
#' cIPMA's theoretical rescaling. It is kept for backward compatibility but
#' should not be used for cIPMA-aligned analysis.
#'
#' @param model A seminr model.
#' @param scales Ignored in this function (uses observed range).
#'
#' @return Data frame with rescaled scores (0-100 based on observed range).
#'
#' @keywords internal
.get_rescaled_scores_internal_DEPRECATED <- function(model, scales) {

  .Deprecated("get_rescaled_lv_scores",
              msg = "This function uses observed min-max rescaling. Use get_rescaled_lv_scores() for cIPMA-aligned theoretical rescaling.")

  model <- .unwrap_seminr_model(model)

  if (is.null(model$construct_scores)) {
    stop("construct_scores not found in seminr model.")
  }

  cs <- model$construct_scores
  rescaled <- as.data.frame(cs)

  for (constr in colnames(cs)) {
    vals <- cs[, constr]
    rng <- range(vals, na.rm = TRUE)

    if (!is.finite(rng[1]) || !is.finite(rng[2])) {
      rescaled[[constr]] <- NA_real_
    } else if (rng[1] == rng[2]) {
      rescaled[[constr]] <- 50
    } else {
      # OBSERVED min-max rescaling (NOT theoretical!)
      rescaled[[constr]] <- (vals - rng[1]) / (rng[2] - rng[1]) * 100
    }
  }

  rescaled
}
