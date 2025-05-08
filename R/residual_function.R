#' @title Compute Residuals from TMB Model
#' @description Computes Pearson residuals for both occupancy and abundance components.
#' @param model Output from run_full_TMB
#' @param data_array_filtered Original input array (Species x Sites x Replicates)
#' @return A list with residual matrices: occupancy_residuals, abundance_residuals
#' @export
#' @param type Character; residual type: "pearson", "response", "deviance"
#' @return List of residual matrices: named by type
compute_residuals_TMB <- function(model, data_array_filtered, type = c("pearson", "response", "deviance")) {
  type <- match.arg(type)

  Y <- to2D(data_array_filtered)
  observed_counts <- as.matrix(Y[, -(1:2)])
  site_ids <- as.numeric(Y$Site)

  rep <- model$TMBobj$report()
  fitted_abundance <- exp(rep$etaa)
  fitted_occupancy_site <- 1 - plogis(rep$etao)
  fitted_occupancy <- fitted_occupancy_site[site_ids, ]

  mu <- fitted_occupancy * fitted_abundance
  y <- observed_counts

  # Default: Pearson
  if (type == "pearson") {
    abund_resid <- (y - mu) / sqrt(pmax(mu, 1e-6))
    binary_pa <- (y > 0) * 1
    occ_resid <- (binary_pa - fitted_occupancy) /
      sqrt(pmax(fitted_occupancy * (1 - fitted_occupancy), 1e-6))
  }

  # Response residuals: raw difference
  else if (type == "response") {
    abund_resid <- y - mu
    occ_resid <- (y > 0) - fitted_occupancy
  }

  # Deviance residuals for ZINB approx (not exact)
  else if (type == "deviance") {
    eps <- 1e-6
    sign_val <- sign(y - mu)
    abund_resid <- sign_val * sqrt(2 * (
      y * log(pmax(y, eps) / pmax(mu, eps)) - (y - mu)
    ))
    binary_pa <- (y > 0) * 1
    occ_resid <- sign(binary_pa - fitted_occupancy) *
      sqrt(2 * (binary_pa * log(pmax(binary_pa, eps) / pmax(fitted_occupancy, eps)) +
                  (1 - binary_pa) * log(pmax(1 - binary_pa, eps) / pmax(1 - fitted_occupancy, eps))))
  }

  return(list(
    occupancy_residuals = occ_resid,
    abundance_residuals = abund_resid
  ))
}
