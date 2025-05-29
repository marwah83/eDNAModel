#' @title Compute Residuals from eDNAModel Fit
#' @description Computes residuals for both abundance and occupancy components of a hierarchical multispecies model fit using TMB.
#' @param model An object of class `eDNAModel` returned by `run_full_TMB()`.
#' @param y A numeric matrix of observed counts (replicates × species). Must match rownames in X.
#' @param X A data.frame of covariates corresponding to rows in `y`. Must include `Site` as a factor.
#' @param type Type of residuals: "pearson" (default), "response", or "deviance".
#' @return A list with two components:
#' \describe{
#'   \item{occupancy_residuals}{A matrix of residuals for the occupancy component.}
#'   \item{abundance_residuals}{A matrix of residuals for the abundance component.}
#' }
#' @export
#' @examples
#' residuals <- compute_residuals_TMB(model = model, y = y, X = X, type = "pearson")

compute_residuals_TMB <- function(model, y, X, type = c("pearson", "response", "deviance")) {
  type <- match.arg(type)

  # Input validation
  if (!is.matrix(y)) stop("'y' must be a matrix.")
  if (!all(rownames(y) %in% rownames(X))) stop("Row names of 'y' must match those in 'X'.")

  # Align X to y
  X <- X[rownames(y), , drop = FALSE]

  # Site index (0-based)
  if (!"Site" %in% colnames(X)) stop("Column 'Site' must exist in X.")
  site_ids <- as.numeric(factor(X$Site)) - 1

  # Extract model predictions
  rep <- model$TMBobj$report()
  if (is.null(rep$etaa) || is.null(rep$etao)) {
    stop("Model report must include 'etaa' (abundance) and 'etao' (occupancy).")
  }

  # Predicted values
  fitted_abundance <- exp(rep$etaa)
  fitted_occupancy_site <- 1 - plogis(rep$etao)
  fitted_occupancy <- fitted_occupancy_site[site_ids + 1, , drop = FALSE]

  mu <- fitted_occupancy * fitted_abundance
  obs <- as.matrix(y)

  # Residual types
  if (type == "pearson") {
    abundance_resid <- (obs - mu) / sqrt(pmax(mu, 1e-6))
    pa_obs <- (obs > 0) * 1
    occ_resid <- (pa_obs - fitted_occupancy) /
      sqrt(pmax(fitted_occupancy * (1 - fitted_occupancy), 1e-6))

  } else if (type == "response") {
    abundance_resid <- obs - mu
    occ_resid <- (obs > 0) - fitted_occupancy

  } else if (type == "deviance") {
    eps <- 1e-6
    sign_val <- sign(obs - mu)
    abundance_resid <- sign_val * sqrt(2 * (
      obs * log(pmax(obs, eps) / pmax(mu, eps)) - (obs - mu)
    ))

    pa_obs <- (obs > 0) * 1
    occ_resid <- sign(pa_obs - fitted_occupancy) * sqrt(
      2 * (
        pa_obs * log(pmax(pa_obs, eps) / pmax(fitted_occupancy, eps)) +
        (1 - pa_obs) * log(pmax(1 - pa_obs, eps) / pmax(1 - fitted_occupancy, eps))
      )
    )
  }

  return(list(
    occupancy_residuals = occ_resid,
    abundance_residuals = abundance_resid
  ))
}
