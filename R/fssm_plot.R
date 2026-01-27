# ==============================================================================
# FSSM Package - Visualization Functions
# ==============================================================================
# Plotting functions for FSSM results including:
# - IPMA-style plots (importance-performance maps)
# - Causal claims comparison plots
# - Dose-response curves
# - Membership distribution plots
# - NCA bottleneck plots
# ==============================================================================


#' Plot FSSM Results (IPMA-Style)
#'
#' @description Visualizes FSSM results using an importance-performance map
#' with optional causal claims overlays. Supports both cIPMA baseline metrics
#' and FSSM-specific causal claims as axes.
#'
#' @param x An fssm object.
#' @param x_axis Character. Variable for x-axis. Options:
#'   \itemize{
#'     \item "importance" (default) - PLS-SEM total effect (cIPMA baseline)
#'     \item "alpha" - Typical sufficiency
#'     \item "delta" - Probabilistic sufficiency (causal uplift)
#'     \item "epsilon" - Displayed as (1-epsilon) so higher = more necessary
#'     \item "beta" - Displayed as (1-beta) so higher = stronger constraint
#'     \item "nca_d" - NCA effect size (cIPMA baseline)
#'   }
#' @param y_axis Character. Variable for y-axis. Same options as x_axis,
#'   plus "performance" (default, cIPMA baseline).
#' @param bubble_size Character. Variable for bubble size. Options:
#'   "performance", "fail_percentage", "violation_mass", "discriminative_power", "none".
#' @param color_by Character. Variable for color grouping. Options:
#'   "necessity", "sufficiency", "claim_count", "none".
#' @param show_quadrants Logical. Show quadrant dividing lines at means. Default TRUE.
#' @param show_labels Logical. Show construct name labels. Default TRUE.
#' @param crisp Logical. Use crisp (TRUE) or fuzzy (FALSE, default) claims.
#' @param theme Character. Plot theme: "minimal" (default), "classic", "bw".
#' @param ... Additional arguments passed to ggplot.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Classic IPMA plot (cIPMA baseline)
#' plot(fssm_result)
#'
#' # Causal claims plot
#' plot(fssm_result, x_axis = "delta", y_axis = "alpha", color_by = "necessity")
#'
#' # NCA-focused plot
#' plot(fssm_result, x_axis = "nca_d", y_axis = "performance",
#'      bubble_size = "fail_percentage")
#' }
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
plot.fssm <- function(x,
                      x_axis = c("importance", "alpha", "delta", "epsilon", "beta", "nca_d"),
                      y_axis = c("performance", "alpha", "delta", "epsilon", "beta", "importance", "nca_d"),
                      bubble_size = c("performance", "fail_percentage", "violation_mass",
                                      "discriminative_power", "none"),
                      color_by = c("necessity", "sufficiency", "claim_count", "none"),
                      show_quadrants = TRUE,
                      show_labels = TRUE,
                      crisp = FALSE,
                      theme = c("minimal", "classic", "bw"),
                      ...) {

  x_axis <- match.arg(x_axis)
  y_axis <- match.arg(y_axis)
  bubble_size <- match.arg(bubble_size)
  color_by <- match.arg(color_by)
  theme <- match.arg(theme)

  df <- x$results
  target <- x$thresholds$target_level

  # ============================================================================
  # Helper: Get claim values
  # ============================================================================

  get_claim <- function(claim_name, use_crisp) {
    prefix <- if (use_crisp) "Crisp_" else "Fuzzy_"
    col_name <- paste0(prefix, claim_name)
    if (col_name %in% names(df)) {
      return(df[[col_name]])
    } else {
      return(rep(NA_real_, nrow(df)))
    }
  }

  # ============================================================================
  # X-axis variable
  # ============================================================================

  x_data <- switch(x_axis,
                   "importance" = {
                     list(var = df$Importance, label = "Importance (Total Effect)")
                   },
                   "alpha" = {
                     list(var = get_claim("Alpha", crisp),
                          label = paste0(if(crisp) "Crisp " else "Fuzzy ", "Typical Sufficiency (\u03b1)"))
                   },
                   "delta" = {
                     list(var = get_claim("Delta", crisp),
                          label = paste0(if(crisp) "Crisp " else "Fuzzy ", "Prob. Sufficiency (\u0394)"))
                   },
                   "epsilon" = {
                     list(var = 1 - get_claim("Epsilon", crisp),
                          label = paste0(if(crisp) "Crisp " else "Fuzzy ", "Necessity Rate (1-\u03b5)"))
                   },
                   "beta" = {
                     list(var = 1 - get_claim("Beta", crisp),
                          label = paste0(if(crisp) "Crisp " else "Fuzzy ", "Constraint Strength (1-\u03b2)"))
                   },
                   "nca_d" = {
                     list(var = df$NCA_d, label = "NCA Effect Size (d)")
                   }
  )

  df$X_var <- x_data$var
  x_label <- x_data$label

  # ============================================================================
  # Y-axis variable
  # ============================================================================

  y_data <- switch(y_axis,
                   "performance" = {
                     list(var = df$Performance, label = "Performance (0-100)")
                   },
                   "importance" = {
                     list(var = df$Importance, label = "Importance (Total Effect)")
                   },
                   "alpha" = {
                     list(var = get_claim("Alpha", crisp),
                          label = paste0(if(crisp) "Crisp " else "Fuzzy ", "Typical Sufficiency (\u03b1)"))
                   },
                   "delta" = {
                     list(var = get_claim("Delta", crisp),
                          label = paste0(if(crisp) "Crisp " else "Fuzzy ", "Prob. Sufficiency (\u0394)"))
                   },
                   "epsilon" = {
                     list(var = 1 - get_claim("Epsilon", crisp),
                          label = paste0(if(crisp) "Crisp " else "Fuzzy ", "Necessity Rate (1-\u03b5)"))
                   },
                   "beta" = {
                     list(var = 1 - get_claim("Beta", crisp),
                          label = paste0(if(crisp) "Crisp " else "Fuzzy ", "Constraint Strength (1-\u03b2)"))
                   },
                   "nca_d" = {
                     list(var = df$NCA_d, label = "NCA Effect Size (d)")
                   }
  )

  df$Y_var <- y_data$var
  y_label <- y_data$label

  # ============================================================================
  # Bubble size
  # ============================================================================

  if (bubble_size == "performance") {
    df$Size_var <- df$Performance
    size_label <- "Performance"
  } else if (bubble_size == "fail_percentage") {
    df$Size_var <- pmax(df$Fail_Percentage, 1)
    size_label <- "% Below Bottleneck"
  } else if (bubble_size == "violation_mass") {
    df$Size_var <- pmax(df$Violation_Mass * 1000, 1)
    size_label <- "Violation Mass"
  } else if (bubble_size == "discriminative_power") {
    df$Size_var <- pmin(df$Discriminative_Power, 100)
    size_label <- "Discriminative Power"
  } else {
    df$Size_var <- 10
    size_label <- NULL
  }

  # ============================================================================
  # Color grouping
  # ============================================================================

  if (color_by == "necessity") {
    df$Color_var <- ifelse(df$Typical_Nec_Supported, "Necessary", "Not Necessary")
    color_label <- "Necessity"
    color_values <- c("Necessary" = "#2166ac", "Not Necessary" = "#b2182b")
  } else if (color_by == "sufficiency") {
    df$Color_var <- ifelse(df$Typical_Suff_Supported, "Sufficient", "Not Sufficient")
    color_label <- "Sufficiency"
    color_values <- c("Sufficient" = "#2166ac", "Not Sufficient" = "#b2182b")
  } else if (color_by == "claim_count") {
    n_claims <- rowSums(cbind(
      df$Typical_Suff_Supported,
      df$Typical_Nec_Supported,
      df$Prob_Suff_Supported,
      df$Prob_Nec_Supported
    ), na.rm = TRUE)
    df$Color_var <- factor(n_claims, levels = 0:4)
    color_label <- "Claims Supported"
    color_values <- c("0" = "#d73027", "1" = "#fc8d59", "2" = "#fee090",
                      "3" = "#91bfdb", "4" = "#4575b4")
  } else {
    df$Color_var <- "All"
    color_label <- NULL
    color_values <- c("All" = "#4575b4")
  }

  # ============================================================================
  # Axis limits and quadrant lines
  # ============================================================================

  x_range <- range(df$X_var, na.rm = TRUE)
  y_range <- range(df$Y_var, na.rm = TRUE)

  x_spread <- max(diff(x_range), 0.1)
  y_spread <- max(diff(y_range), 10)

  x_pad <- x_spread * 0.15
  y_pad <- y_spread * 0.15

  x_limits <- c(x_range[1] - x_pad, x_range[2] + x_pad)

  # For performance, constrain to 0-100
  if (y_axis == "performance") {
    y_limits <- c(max(0, y_range[1] - y_pad), min(100, y_range[2] + y_pad))
  } else {
    y_limits <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  }

  x_mean <- mean(df$X_var, na.rm = TRUE)
  y_mean <- mean(df$Y_var, na.rm = TRUE)

  # ============================================================================
  # Build plot
  # ============================================================================

  p <- ggplot2::ggplot(df, ggplot2::aes(x = X_var, y = Y_var))

  # Quadrant lines
  if (show_quadrants) {
    p <- p +
      ggplot2::geom_vline(xintercept = x_mean, linetype = "dashed",
                          color = "grey50", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = y_mean, linetype = "dashed",
                          color = "grey50", linewidth = 0.5)
  }

  # Points
  if (bubble_size != "none" && !is.null(size_label)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(size = Size_var, fill = Color_var),
      shape = 21, color = "black", stroke = 0.8, alpha = 0.85
    )
  } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(fill = Color_var),
      shape = 21, size = 10, color = "black", stroke = 0.8, alpha = 0.85
    )
  }

  # Labels
  if (show_labels) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = Construct),
      size = 3.5, box.padding = 0.6, point.padding = 0.5,
      max.overlaps = 20, seed = 42
    )
  }

  # Scales
  p <- p + ggplot2::scale_fill_manual(values = color_values, name = color_label)

  if (bubble_size != "none" && !is.null(size_label)) {
    p <- p + ggplot2::scale_size_continuous(range = c(5, 15), name = size_label)
  }

  # Coordinates and labels
  p <- p +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits) +
    ggplot2::labs(
      title = paste0("FSSM Analysis (Target Level = ", target, ")"),
      subtitle = create_plot_subtitle(x_axis, y_axis, color_by, crisp),
      x = x_label,
      y = y_label
    )

  # Theme
  p <- p + switch(theme,
                  "minimal" = ggplot2::theme_minimal(),
                  "classic" = ggplot2::theme_classic(),
                  "bw" = ggplot2::theme_bw()
  )

  p <- p + ggplot2::theme(
    legend.position = "right",
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10, color = "grey40")
  )

  p
}


#' Create Plot Subtitle
#' @keywords internal
create_plot_subtitle <- function(x_axis, y_axis, color_by, crisp) {

  parts <- character(0)

  # Formulation type
  if (x_axis %in% c("alpha", "delta", "epsilon", "beta") ||
      y_axis %in% c("alpha", "delta", "epsilon", "beta")) {
    parts <- c(parts, if (crisp) "Crisp formulation" else "Fuzzy formulation")
  }

  # Color legend hint
  if (color_by == "necessity") {
    parts <- c(parts, "Blue = Necessary")
  } else if (color_by == "sufficiency") {
    parts <- c(parts, "Blue = Sufficient")
  } else if (color_by == "claim_count") {
    parts <- c(parts, "Color = # claims supported (0-4)")
  }

  if (length(parts) == 0) return(NULL)
  paste(parts, collapse = " | ")
}


#' Plot Dose-Response Curves
#'
#' @description Visualizes the dose-response relationship between sufficiency
#' membership and expected outcome for each predictor.
#'
#' @param x An fssm object.
#' @param constructs Character vector. Which constructs to plot. Default is all.
#' @param show_points Logical. Show endpoint markers. Default TRUE.
#' @param show_ci Logical. Show confidence bands (if available). Default FALSE.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plot_dose_response <- function(x, constructs = NULL, show_points = TRUE, show_ci = FALSE) {

  if (is.null(constructs)) constructs <- x$predictors

  # Validate constructs
  invalid <- constructs[!constructs %in% x$predictors]
  if (length(invalid) > 0) {
    warning("Constructs not found: ", paste(invalid, collapse = ", "))
    constructs <- constructs[constructs %in% x$predictors]
  }

  if (length(constructs) == 0) {
    stop("No valid constructs to plot.")
  }

  # Compile dose-response data
  dr_data <- do.call(rbind, lapply(constructs, function(pred) {
    dr <- x$claims[[pred]]$dose_response$predictions
    dr$Construct <- pred
    dr$Uplift <- x$claims[[pred]]$dose_response$uplift
    dr
  }))

  # Key endpoints
  key_points <- do.call(rbind, lapply(constructs, function(pred) {
    dr <- x$claims[[pred]]$dose_response
    data.frame(
      Construct = pred,
      mu_S = c(0, 1),
      E_mu_Y = c(dr$E_mu_Y_at_0, dr$E_mu_Y_at_1),
      Label = c("Baseline", "Full Effect"),
      stringsAsFactors = FALSE
    )
  }))

  # Build plot
  p <- ggplot2::ggplot(dr_data, ggplot2::aes(x = mu_S, y = E_mu_Y, color = Construct)) +
    ggplot2::geom_line(linewidth = 1.2)

  if (show_points) {
    p <- p + ggplot2::geom_point(data = key_points, size = 4)
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, 0.25),
      labels = c("0%", "25%", "50%", "75%", "100%")
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = "Dose-Response: Sufficiency vs Expected Outcome",
      subtitle = "Shows causal effect of increasing sufficiency membership",
      x = expression("Sufficiency Membership " * (mu[S])),
      y = expression("Expected Outcome " * E * "[" * mu[Y] * "]"),
      color = "Construct"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}


#' Plot Claims Comparison (Heatmap or Bar)
#'
#' @description Creates a comparison chart of all four causal claims across constructs.
#'
#' @param x An fssm object.
#' @param type Character. Plot type: "heatmap" (default) or "bar".
#' @param crisp Logical. Show crisp (TRUE) or fuzzy (FALSE, default) claims.
#' @param show_thresholds Logical. Show threshold lines (bar plot only). Default TRUE.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plot_claims_comparison <- function(x,
                                   type = c("heatmap", "bar"),
                                   crisp = FALSE,
                                   show_thresholds = TRUE) {

  type <- match.arg(type)
  df <- x$results
  thresholds <- x$thresholds

  # Prepare long format data
  prefix <- if (crisp) "Crisp_" else "Fuzzy_"

  # Transform epsilon and beta so higher = better for visualization
  claims_long <- data.frame(
    Construct = rep(df$Construct, 4),
    Claim = factor(
      rep(c("Alpha\n(Typ. Suff.)", "1-Epsilon\n(Typ. Nec.)",
            "Delta\n(Prob. Suff.)", "1-Beta\n(Prob. Nec.)"), each = nrow(df)),
      levels = c("Alpha\n(Typ. Suff.)", "1-Epsilon\n(Typ. Nec.)",
                 "Delta\n(Prob. Suff.)", "1-Beta\n(Prob. Nec.)")
    ),
    Value = c(
      df[[paste0(prefix, "Alpha")]],
      1 - df[[paste0(prefix, "Epsilon")]],
      df[[paste0(prefix, "Delta")]],
      1 - df[[paste0(prefix, "Beta")]]
    ),
    stringsAsFactors = FALSE
  )

  # Thresholds (transformed)
  threshold_values <- c(
    thresholds$alpha,
    1 - thresholds$epsilon,
    thresholds$delta,
    1 - thresholds$beta
  )

  claims_long$Threshold <- rep(threshold_values, each = nrow(df))
  claims_long$Supported <- claims_long$Value >= claims_long$Threshold

  # ============================================================================
  # Heatmap
  # ============================================================================

  if (type == "heatmap") {
    p <- ggplot2::ggplot(claims_long,
                         ggplot2::aes(x = Claim, y = Construct, fill = Value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", Value)),
        size = 3.5, color = "black"
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#d73027", mid = "#ffffbf", high = "#1a9850",
        midpoint = 0.5, limits = c(0, 1), name = "Value"
      ) +
      ggplot2::labs(
        title = paste("Causal Claims Comparison", if(crisp) "(Crisp)" else "(Fuzzy)"),
        subtitle = "Values transformed so higher = better for all claims",
        x = "", y = ""
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 9),
        panel.grid = ggplot2::element_blank()
      )

    # ============================================================================
    # Bar plot
    # ============================================================================

  } else {

    # Add threshold data for reference lines
    threshold_df <- data.frame(
      Claim = factor(
        c("Alpha\n(Typ. Suff.)", "1-Epsilon\n(Typ. Nec.)",
          "Delta\n(Prob. Suff.)", "1-Beta\n(Prob. Nec.)"),
        levels = levels(claims_long$Claim)
      ),
      Threshold = threshold_values
    )

    p <- ggplot2::ggplot(claims_long,
                         ggplot2::aes(x = Construct, y = Value, fill = Supported)) +
      ggplot2::geom_col(color = "black", width = 0.7, linewidth = 0.3) +
      ggplot2::facet_wrap(~ Claim, scales = "free_y", ncol = 2)

    if (show_thresholds) {
      p <- p + ggplot2::geom_hline(
        data = threshold_df,
        ggplot2::aes(yintercept = Threshold),
        linetype = "dashed", color = "blue", linewidth = 0.8
      )
    }

    p <- p +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#1a9850", "FALSE" = "#d73027"),
        labels = c("TRUE" = "Supported", "FALSE" = "Not Supported"),
        name = "Claim Status"
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        title = paste("Causal Claims Comparison", if(crisp) "(Crisp)" else "(Fuzzy)"),
        subtitle = "Blue dashed line = threshold. Values transformed so higher = better.",
        x = "", y = "Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )
  }

  p
}


#' Plot Membership Distribution
#'
#' @description Visualizes the fuzzy membership distribution for a single construct,
#' showing the NCA ceiling line and case positions.
#'
#' @param x An fssm object.
#' @param construct Character. Which construct to visualize.
#' @param color_by Character. What determines point color:
#'   "mu_S" (sufficiency, default), "mu_N" (necessity), "mu_Y" (outcome), "weight".
#' @param show_ceiling Logical. Show NCA ceiling line. Default TRUE.
#' @param show_thresholds Logical. Show threshold lines. Default TRUE.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plot_membership <- function(x,
                            construct,
                            color_by = c("mu_S", "mu_N", "mu_Y", "weight"),
                            show_ceiling = TRUE,
                            show_thresholds = TRUE) {

  color_by <- match.arg(color_by)

  if (!construct %in% x$predictors) {
    stop(paste("Construct", construct, "not found. Available:",
               paste(x$predictors, collapse = ", ")))
  }

  membership <- x$membership[[construct]]
  weights <- x$weights[[construct]]

  df <- data.frame(
    X_star = x$data_rescaled[[construct]],
    Y_star = x$data_rescaled[[x$target_construct]],
    mu_S = membership$mu_S,
    mu_N = membership$mu_N,
    mu_Y = membership$mu_Y,
    weight = weights
  )

  # Color variable
  color_var <- df[[color_by]]
  color_label <- switch(color_by,
                        "mu_S" = expression(mu[S] ~ "(Sufficiency)"),
                        "mu_N" = expression(mu[N] ~ "(Necessity)"),
                        "mu_Y" = expression(mu[Y] ~ "(Outcome)"),
                        "weight" = "Case Weight"
  )

  # Build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = X_star, y = Y_star))

  # NCA ceiling line
  if (show_ceiling) {
    p <- p + ggplot2::geom_abline(
      intercept = 0, slope = membership$bottleneck_slope,
      linetype = "dashed", color = "red", linewidth = 1
    )
  }

  # Threshold lines
  if (show_thresholds) {
    p <- p +
      ggplot2::geom_vline(xintercept = membership$T_X, linetype = "dotted",
                          color = "blue", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = membership$T_Y, linetype = "dotted",
                          color = "blue", linewidth = 0.8)
  }

  # Points
  p <- p + ggplot2::geom_point(
    ggplot2::aes_string(color = color_by),
    size = 3, alpha = 0.7
  )

  # Color scale
  p <- p + ggplot2::scale_color_gradient2(
    low = "#d73027", mid = "#ffffbf", high = "#1a9850",
    midpoint = 0.5, limits = c(0, 1), name = color_label
  )

  # Labels and theme
  p <- p +
    ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::labs(
      title = paste("Membership Distribution:", construct),
      subtitle = paste0(
        "Red dashed = NCA ceiling (slope = ", round(membership$bottleneck_slope, 3), "). ",
        "Blue dotted = thresholds (T_X = ", round(membership$T_X, 1),
        ", T_Y = ", membership$T_Y, ")."
      ),
      x = paste0(construct, " (0-100)"),
      y = paste0(x$target_construct, " (0-100)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}


#' Plot NCA Bottleneck Table
#'
#' @description Visualizes the NCA bottleneck levels for all predictors,
#' showing required X levels at different Y targets.
#'
#' @param x An fssm object.
#' @param highlight_target Logical. Highlight the target level row. Default TRUE.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plot_bottleneck <- function(x, highlight_target = TRUE) {

  bn_full <- x$bottleneck_full

  if (is.null(bn_full)) {
    stop("Bottleneck table not available in fssm object.")
  }

  # Convert to long format
  bn_df <- as.data.frame(bn_full)

  # Get Y column (first column is usually the target values)
  y_col <- names(bn_df)[1]
  y_values <- as.numeric(as.character(bn_df[[y_col]]))

  bn_long <- do.call(rbind, lapply(x$predictors, function(pred) {
    if (!pred %in% names(bn_df)) return(NULL)

    data.frame(
      Y_level = y_values,
      Construct = pred,
      Required_X = as.numeric(as.character(bn_df[[pred]])),
      stringsAsFactors = FALSE
    )
  }))

  bn_long <- bn_long[!is.na(bn_long$Required_X), ]

  # Build plot
  p <- ggplot2::ggplot(bn_long, ggplot2::aes(x = Y_level, y = Required_X,
                                             color = Construct, group = Construct)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2)

  # Highlight target level
  if (highlight_target) {
    target <- x$thresholds$target_level
    p <- p + ggplot2::geom_vline(
      xintercept = target, linetype = "dashed", color = "grey30", linewidth = 0.8
    ) +
      ggplot2::annotate("text", x = target, y = max(bn_long$Required_X, na.rm = TRUE),
                        label = paste0("Target = ", target), hjust = -0.1, vjust = 1)
  }

  p <- p +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 10)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 10)) +
    ggplot2::labs(
      title = "NCA Bottleneck Analysis",
      subtitle = "Required X level to achieve each Y level",
      x = paste0(x$target_construct, " Level (%)"),
      y = "Required Predictor Level (%)",
      color = "Construct"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  p
}


#' Plot FSSM Summary Dashboard
#'
#' @description Creates a multi-panel summary of FSSM results combining
#' IPMA, claims, and dose-response plots.
#'
#' @param x An fssm object.
#'
#' @return A combined ggplot object (requires patchwork package).
#'
#' @export
plot_dashboard <- function(x) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required for dashboard. Install with: install.packages('patchwork')")
  }

  # Create individual plots
  p1 <- plot.fssm(x, x_axis = "importance", y_axis = "performance",
                  color_by = "necessity", show_labels = TRUE) +
    ggplot2::ggtitle("IPMA: Importance vs Performance")

  p2 <- plot_claims_comparison(x, type = "heatmap") +
    ggplot2::ggtitle("Causal Claims Heatmap")

  p3 <- plot_dose_response(x) +
    ggplot2::ggtitle("Dose-Response Curves")

  p4 <- plot_bottleneck(x) +
    ggplot2::ggtitle("NCA Bottlenecks")

  # Combine using patchwork
  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)
}


# ==============================================================================
# Global Variables Declaration (for R CMD check)
# ==============================================================================

utils::globalVariables(c(
  "X_var", "Y_var", "Size_var", "Color_var", "Construct",
  "mu_S", "E_mu_Y", "Claim", "Value", "Supported", "Threshold",
  "X_star", "Y_star", "weight", "mu_N", "mu_Y",
  "Y_level", "Required_X", "Label", "Uplift"
))
