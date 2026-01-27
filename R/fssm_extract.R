#' @title Master FSSM Data Extractor
#'
#' @description Prepares a unified FSSM object from either a PLS-SEM model or raw data.
#' Implements Step 1 of the FSSM workflow:
#' - PLS-SEM route: Uses outer weights + theoretical min/max rescaling (cIPMA-identical)
#' - Raw data route: Uses entropy/equal weighting + observed or theoretical rescaling
#'
#' @param input Either a `boot_seminr_model` OR a `data.frame`.
#' @param syntax (Required for Raw Data) A string defining constructs (e.g., "IM =~ item1 + item2").
#' @param scales Optional. A named list of min/max scales (e.g., list(PU_01=c(1,5))).
#'   For SEM: Required for 0-100 normalization.
#'   For Raw: If provided, uses theoretical rescaling; if NULL, uses observed min/max.
#' @param data (Optional for SEM) Raw indicator data. If NULL, extracted from model.
#' @param target_construct String. Name of the outcome variable (Y).
#' @param treatment_construct String. Name of the treatment variable (X).
#' @param indicator_weighting String. How to weight indicators for raw data route:
#'   - "entropy" (default): Shannon entropy weights (more discriminating items get higher weight)
#'   - "equal": Equal weights (simple mean of items)
#'
#' @return An object of class `fssm_object`.
#' @export
fssm_extract <- function(input,
                         syntax = NULL,
                         scales = NULL,
                         data = NULL,
                         target_construct = NULL,
                         treatment_construct = NULL,
                         indicator_weighting = c("entropy", "equal")) {

  indicator_weighting <- match.arg(indicator_weighting)

  # ============================================================================
  # ROUTE A: PLS-SEM MODEL (Use PLS weights + Theoretical Rescaling)
  # ============================================================================
  if (inherits(input, "seminr_model") || inherits(input, "boot_seminr_model")) {

    message(">>> FSSM: Processing PLS-SEM Model...")

    if (is.null(scales)) {
      stop("Error: PLS-SEM models require 'scales' for 0-100 normalization.")
    }

    # Extract raw data from model if not provided
    if (is.null(data)) {
      raw_data <- input$data
    } else {
      raw_data <- data
    }

    if (is.null(raw_data)) {
      stop("Error: Could not extract raw data from model. Please provide 'data' argument.")
    }

    # Rescale to 0-100 using cIPMA-identical method (theoretical min/max)
    scores_df <- .get_rescaled_scores_cipma(input, raw_data, scales)

    meta <- list(type = "PLS-SEM", weighting = "Outer Weights (PLS)")
  }

  # ============================================================================
  # ROUTE B: RAW DATA (Use Entropy/Equal Weights + Rescaling)
  # ============================================================================
  else if (is.data.frame(input)) {

    weighting_label <- if (indicator_weighting == "entropy") {
      "Entropy (Shannon)"
    } else {
      "Equal Weights (Mean)"
    }
    message(paste0(">>> FSSM: Processing Raw Data (", weighting_label, ")..."))

    if (is.null(syntax)) {
      stop("Error: Raw data requires 'syntax' to define constructs.")
    }

    # Parse syntax to get construct-indicator mapping
    mapping <- .parse_mcdm_syntax(syntax)

    # Calculate weighted scores for each construct
    scores_list <- list()

    for (constr in names(mapping)) {
      items <- mapping[[constr]]

      # Check if items exist in data
      missing_items <- items[!items %in% names(input)]
      if (length(missing_items) > 0) {
        stop(paste("Items missing from data:", paste(missing_items, collapse = ", ")))
      }

      item_data <- input[, items, drop = FALSE]

      # Simple mean imputation for missing values
      if (any(is.na(item_data))) {
        for (i in seq_len(ncol(item_data))) {
          na_idx <- is.na(item_data[, i])
          if (any(na_idx)) {
            item_data[na_idx, i] <- mean(item_data[, i], na.rm = TRUE)
          }
        }
      }

      # Calculate composite score based on weighting method
      if (indicator_weighting == "entropy") {
        scores_list[[constr]] <- .calculate_entropy_score(item_data, scales)
      } else {
        # Equal weights (simple mean)
        scores_list[[constr]] <- .calculate_equal_score(item_data, scales)
      }
    }

    scores_df <- as.data.frame(scores_list)
    raw_data <- input

    meta <- list(type = "Raw-Composite", weighting = weighting_label)

  } else {
    stop("Input must be a seminr model or a data frame.")
  }

  # ============================================================================
  # VALIDATION
  # ============================================================================
  if (!is.null(target_construct)) {
    if (!target_construct %in% names(scores_df)) {
      stop(paste("Target construct", target_construct, "not found in scores."))
    }
  }

  if (!is.null(treatment_construct)) {
    if (!treatment_construct %in% names(scores_df)) {
      stop(paste("Treatment construct", treatment_construct, "not found in scores."))
    }
  }

  # ============================================================================
  # BUILD OUTPUT
  # ============================================================================
  structure(
    list(
      scores = scores_df,
      raw_data = raw_data,
      target = target_construct,
      treatment = treatment_construct,
      metadata = meta
    ),
    class = "fssm_object"
  )
}


# ==============================================================================
# INTERNAL HELPERS (used only by fssm_extract)
# ==============================================================================

#' Rescale PLS Scores Using cIPMA Method
#'
#' @description Uses the exact cIPMA rescaling procedure:
#' 1. Get standardized outer weights
#' 2. Unstandardize by dividing by indicator SD
#' 3. Normalize within each LV
#' 4. Compute weighted sum of original indicators
#' 5. Rescale using THEORETICAL min/max from scales
#'
#' @param model A seminr model (base or bootstrapped).
#' @param data Raw indicator data frame.
#' @param scales Named list of theoretical min/max for indicators.
#'
#' @return Data frame with rescaled LV scores (0-100).
#' @keywords internal
.get_rescaled_scores_cipma <- function(model, data, scales) {


  # Unwrap bootstrapped model if needed
  base_model <- model
  if (inherits(model, "boot_seminr_model")) {
    if (!is.null(model$original_model)) {
      base_model <- model$original_model
    } else if (!is.null(model$model)) {
      base_model <- model$model
    }
  }

  if (is.null(base_model$outer_weights)) {
    stop("Outer weights not found in model.")
  }

  # Prepare indicator data
  original_indicators <- as.data.frame(data)

  # Get standardized outer weights as data frame
  std_weights_matrix <- base_model$outer_weights
  std_weights_df <- as.data.frame(as.table(std_weights_matrix))
  colnames(std_weights_df) <- c("Indicator", "Latent_variable", "Standardized_weight")

  # Keep only non-zero weights
  std_weights_df <- std_weights_df[std_weights_df$Standardized_weight != 0, ]

  # Check for negative weights
  negative_weights <- std_weights_df$Standardized_weight < 0
  if (any(negative_weights)) {
    neg_indicators <- std_weights_df$Indicator[negative_weights]
    warning(paste0(
      "Negative outer weights detected for: ",
      paste(neg_indicators, collapse = ", "),
      ". This may invert performance scores."
    ))
  }

  # Get standard deviations (use stored if available)
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

  # Normalize within each LV
  sum_weights <- tapply(
    std_weights_df$Unstandardized_weight,
    std_weights_df$Latent_variable,
    sum
  )
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

      # Theoretical min/max rescaling (cIPMA-identical)
      norm_val <- ((lv_unscaled[[lv]] - min_val) / (max_val - min_val)) * 100

      # Clip to [0, 100]
      norm_val <- pmax(0, pmin(100, norm_val))

      lv_df[[lv]] <- norm_val
    } else {
      warning(paste0(
        "No scale found for construct '", lv, "'. Scores not rescaled."
      ))
    }
  }

  lv_df
}


#' Internal Syntax Parser for Entropy Mode
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


#' Calculate Entropy-Weighted Composite Score (0-100)
#'
#' @description For raw data route. Uses Shannon entropy to weight indicators,
#' then rescales to 0-100.
#'
#' @param df Data frame of indicator values for one construct.
#' @param scales Optional. Named list of min/max for theoretical rescaling.
#'   If NULL, uses observed min/max.
#'
#' @keywords internal
.calculate_entropy_score <- function(df, scales = NULL) {

  mtx <- as.matrix(df)
  original_mtx <- mtx  # Keep original for rescaling

  # Ensure positive values for log
  if (min(mtx, na.rm = TRUE) <= 0) {
    mtx <- mtx + abs(min(mtx, na.rm = TRUE)) + 0.001
  }

  # Normalize to proportions
  col_sums <- colSums(mtx, na.rm = TRUE)
  P_ij <- sweep(mtx, 2, col_sums, "/")

  # Entropy calculation
  m <- nrow(mtx)
  k <- 1 / log(m)

  P_log_P <- P_ij * log(P_ij)
  P_log_P[is.na(P_log_P) | is.infinite(P_log_P)] <- 0

  E_j <- -k * colSums(P_log_P, na.rm = TRUE)

  # Diversity and weights
  d_j <- 1 - E_j
  w_j <- d_j / sum(d_j)

  # Weighted composite (use original values for rescaling)
  weighted_sum <- as.vector(original_mtx %*% w_j)

  # Rescale to 0-100
  .rescale_composite(weighted_sum, colnames(df), scales)
}


#' Calculate Equal-Weighted Composite Score (0-100)
#'
#' @description For raw data route. Uses equal weights (simple mean) for indicators,
#' then rescales to 0-100.
#'
#' @param df Data frame of indicator values for one construct.
#' @param scales Optional. Named list of min/max for theoretical rescaling.
#'   If NULL, uses observed min/max.
#'
#' @keywords internal
.calculate_equal_score <- function(df, scales = NULL) {

  mtx <- as.matrix(df)
  n_items <- ncol(mtx)

  # Equal weights
  w_j <- rep(1 / n_items, n_items)

  # Weighted composite (= simple mean)
  weighted_sum <- as.vector(mtx %*% w_j)

  # Rescale to 0-100
  .rescale_composite(weighted_sum, colnames(df), scales)
}


#' Rescale Composite Score to 0-100
#'
#' @description Helper function to rescale composite scores.
#' Uses theoretical min/max if scales provided, otherwise observed.
#'
#' @keywords internal
.rescale_composite <- function(scores, item_names, scales = NULL) {

  # Try to get theoretical min/max from scales
  if (!is.null(scales)) {
    # Find first matching item in scales
    matched_item <- item_names[item_names %in% names(scales)][1]

    if (!is.na(matched_item)) {
      min_val <- scales[[matched_item]][1]
      max_val <- scales[[matched_item]][2]

      # Rescale using theoretical range
      rescaled <- ((scores - min_val) / (max_val - min_val)) * 100
      rescaled <- pmax(0, pmin(100, rescaled))
      return(rescaled)
    }
  }

  # Fallback: observed min/max
  min_s <- min(scores, na.rm = TRUE)
  max_s <- max(scores, na.rm = TRUE)

  if (max_s == min_s) {
    return(rep(50, length(scores)))
  }

  ((scores - min_s) / (max_s - min_s)) * 100
}


#' Print method for fssm_object
#' @export
print.fssm_object <- function(x, ...) {

  cat("\n=== FSSM Data Object ===\n\n")
  cat("Type:", x$metadata$type, "\n")
  cat("Weighting:", x$metadata$weighting, "\n")
  cat("Constructs:", ncol(x$scores), "\n")
  cat("Sample Size:", nrow(x$scores), "\n")

  if (!is.null(x$target)) {
    cat("Target:", x$target, "\n")
  }
  if (!is.null(x$treatment)) {
    cat("Treatment:", x$treatment, "\n")
  }

  cat("\nScore Summary (0-100 scale):\n")

  # Create numeric summary
  score_summary <- data.frame(
    Construct = names(x$scores),
    Min = sapply(x$scores, min, na.rm = TRUE),
    Mean = sapply(x$scores, mean, na.rm = TRUE),
    Median = sapply(x$scores, median, na.rm = TRUE),
    Max = sapply(x$scores, max, na.rm = TRUE),
    SD = sapply(x$scores, sd, na.rm = TRUE)
  )
  rownames(score_summary) <- NULL

  print(score_summary, digits = 2, row.names = FALSE)

  invisible(x)
}
