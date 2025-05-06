#' @title Compute Residuals from TMB Model
#' @description Computes Pearson residuals for both occupancy and abundance components.
#' @param model Output from run_full_TMB
#' @param data_array_filtered Original input array (Species x Sites x Replicates)
#' @return A list with residual matrices: occupancy_residuals, abundance_residuals
#' @export
compute_residuals_TMB <- function(model, data_array_filtered) {
  # Flatten to 2D: get Y and observed counts
  Y <- to2D(data_array_filtered)
  observed_counts <- as.matrix(Y[, -(1:2)])
  
  # Get number of observations and species
  n_obs <- nrow(observed_counts)
  n_species <- ncol(observed_counts)
  
  # Get model outputs
  rep <- model$TMBobj$report()
  fitted_abundance <- exp(rep$etaa)  # (n_obs x species)
  fitted_occupancy_site <- 1 - plogis(rep$etao)  # (n_sites x species)
  
  # Recover site IDs per row in Y
  site_ids <- as.numeric(Y$Site)
  
  # Expand site-level occupancy to match each replicate (per row in Y)
  fitted_occupancy <- fitted_occupancy_site[site_ids, ]
  
  # Now safe to multiply: both are (n_obs x species)
  expected_counts <- fitted_occupancy * fitted_abundance
  
  # Compute residuals
  abundance_residuals <- (observed_counts - expected_counts) / sqrt(pmax(expected_counts, 1e-6))
  
  # Binary PA matrix
  binary_pa <- (observed_counts > 0) * 1
  occupancy_residuals <- (binary_pa - fitted_occupancy) /
    sqrt(pmax(fitted_occupancy * (1 - fitted_occupancy), 1e-6))
  
  return(list(
    occupancy_residuals = occupancy_residuals,
    abundance_residuals = abundance_residuals
  ))
}
