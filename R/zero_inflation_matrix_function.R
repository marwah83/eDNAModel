# ✅ Function to Create Zero-Inflation Matrix
create_zero_inflation_matrix <- function(species, sites, ZIP) {
  return(matrix(ZIP, nrow = species, ncol = sites))  # Fixed ZIP value
}
