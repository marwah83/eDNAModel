#' @title Convert 3D array to 2D data frame (safe version)
#'
#' @param array_3d A 3D array with dimensions (Species, Sites, Replicates) and appropriate dimnames.
#' @return A data frame with columns: Site (numeric), Replicate, Species1, Species2, ..., SpeciesN
#' @export
to2D <- function(array_3d) {
  species <- dim(array_3d)[1]
  sites <- dim(array_3d)[2]
  replicates <- dim(array_3d)[3]

  names(dimnames(array_3d)) <- c("Species", "Sites", "Replicates")

  # Use character site names directly without parsing to numeric
  site_column <- rep(dimnames(array_3d)$Sites, each = replicates)
  replicate_column <- rep(dimnames(array_3d)$Replicates, times = sites)

  # Flatten 3D array to matrix
  data_values <- matrix(array_3d, nrow = sites * replicates, ncol = species, byrow = FALSE)

  # Combine into data frame
  final_data <- data.frame(Site = site_column, Replicate = replicate_column, data_values)

  # Assign species column names
  colnames(final_data)[3:(2 + species)] <- dimnames(array_3d)$Species

  return(final_data)
}
