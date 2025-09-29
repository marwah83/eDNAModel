# ✅ Function to Simulate ZIP Data
simulate_data <- function(species, sites, replicates, lambda, zero_inflation_prob_matrix) {
  tryCatch(
    simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix),
    error = function(e) return(NULL)
  )
}
