#' @title Extract Fitted Values from TMB Model
#' @description Returns fitted values (expected values) for occupancy and abundance components.
#' @param model Output from run_full_TMB
#' @return A list with matrices: fitted_occupancy and fitted_abundance
#' @export
fitted_TMB <- function(model) {
  rep <- model$TMBobj$report()

  # Occupancy component: etao = linear predictor on logit scale
  fitted_occupancy <- 1 - plogis(rep$etao)  # apply inverse link

  # Abundance component: etaa = log(lambda) => apply exp()
  fitted_abundance <- exp(rep$etaa)

  # Assign proper species names if available
  species_names <- colnames(model$TMBobj$env$data$Y)
  if (!is.null(species_names)) {
    colnames(fitted_abundance) <- species_names
  }

  species_names_occ <- colnames(model$TMBobj$env$data$Ysites)
  if (!is.null(species_names_occ)) {
    colnames(fitted_occupancy) <- species_names_occ
  }

  return(list(
    fitted_occupancy = as.matrix(fitted_occupancy),
    fitted_abundance = as.matrix(fitted_abundance)
  ))
}


