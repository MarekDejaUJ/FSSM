#' Run FSSM Proof of Concept (GPS vs Binary)
#'
#' @description This function runs the 20-case simulation demonstrating the
#' difference between GPS weighting and Binary IPW in the FSSM context.
#' It generates the exact tables found in the documentation/paper.
#'
#' @param verbose Logical. If TRUE, prints detailed step-by-step tables.
#' @return A list containing the proof of concept results.
#' @export
run_fssm_proof_of_concept <- function(verbose = TRUE) {

  # 1. Create Data
  data <- data.frame(
    Case = 1:20,
    X_star = c(82, 75, 91, 68, 55, 88, 45, 79, 62, 94, 38, 71, 85, 52, 77, 43, 89, 66, 81, 59),
    Y_star = c(88, 79, 94, 71, 62, 91, 48, 85, 58, 96, 42, 74, 72, 88, 82, 51, 92, 69, 86, 55),
    Role = factor(c("Faculty", "Faculty", "Faculty", "Staff", "Staff", "Faculty", "Staff", "Faculty",
                    "Staff", "Faculty", "Staff", "Staff", "Faculty", "Staff", "Faculty", "Staff",
                    "Faculty", "Staff", "Faculty", "Staff")),
    IPW_weight = c(1.12, 1.08, 0.95, 0.89, 0.94, 1.05, 1.02, 1.10, 0.91, 0.98,
                   1.15, 0.93, 1.04, 0.87, 1.06, 1.08, 0.96, 0.92, 1.03, 0.97)
  )

  # 2. Parameters
  T_Y <- 85; T_X <- 70

  # 3. Membership
  membership <- calculate_membership(data$X_star, data$Y_star, T_X, T_Y)

  # 4. GPS Weights
  gps_res <- calculate_gps_weights(membership$mu_S, data[, "Role", drop=FALSE])

  # 5. Claims (GPS)
  claims_gps <- calculate_all_claims(membership, data$X_star, data$Y_star, T_X, T_Y, gps_res$weights)

  # 6. Claims (Original Binary IPW) - Normalized
  w_orig <- data$IPW_weight / sum(data$IPW_weight)
  claims_orig <- calculate_all_claims(membership, data$X_star, data$Y_star, T_X, T_Y, w_orig)

  if(verbose) {
    cat("=== FSSM Proof of Concept ===\n")
    cat("Comparing GPS vs Binary IPW...\n\n")
    cat(sprintf("Alpha (GPS): %.3f vs (Binary): %.3f\n", claims_gps$fuzzy$alpha, claims_orig$fuzzy$alpha))
    cat(sprintf("Delta (GPS): %.3f vs (Binary): %.3f\n", claims_gps$fuzzy$delta, claims_orig$fuzzy$delta))
  }

  return(list(gps = claims_gps, original = claims_orig, data = data))
}
