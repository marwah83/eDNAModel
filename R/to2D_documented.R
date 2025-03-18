#' Convert 3D array to 2D data frame (safe version)
#'
#' @param array_3d A 3D array with dimensions (Species, Sites, Replicates) and appropriate dimnames.
#' @return A data frame with columns: Site (numeric), Replicate, Species1, Species2, ..., SpeciesN
to2D <- function(array_3d) {
  # Extract dimensions
  species <- dim(array_3d)[1]
  sites <- dim(array_3d)[2]
  replicates <- dim(array_3d)[3]

  names(dimnames(array_3d)) <- c("Species",  "Sites", "Replicates")

  # Check dimnames exist
  if (is.null(dimnames(array_3d)) ||
      is.null(dimnames(array_3d)$Species) ||
      is.null(dimnames(array_3d)$Sites) ||
      is.null(dimnames(array_3d)$Replicates)) {
    stop("❌ data_array_filtered must have dimnames: 'Species', 'Sites', 'Replicates'.")
  }

  # ✅ Extract and parse site names to numeric
  site_column <- rep(dimnames(array_3d)$Sites, each = replicates)
  site_column_numeric <- as.numeric(gsub("\\D", "", site_column))  # Extract digits safely

  # Replicate column as factor or character
  replicate_column <- rep(dimnames(array_3d)$Replicates, times = sites)

  # Flatten the 3D array to 2D matrix
  data_values <- matrix(array_3d, nrow = sites * replicates, ncol = species, byrow = FALSE)

  # Combine into final data frame
  final_data <- data.frame(Site = site_column_numeric, Replicate = replicate_column, data_values)

  # Name species columns properly
  colnames(final_data)[3:(2 + species)] <- dimnames(array_3d)$Species

  # ✅ Optional: check for NAs in site parsing and warn
  if (any(is.na(final_data$Site))) {
    warning("⚠️ Some site names could not be parsed to numeric. Check the 'Sites' dimnames formatting.")
  }

  return(final_data)
}
