test_that("simulate_zip_data correctly applies zero-inflation", {
  species <- 3
  sites <- 3
  replicates <- 3
  lambda <- 5
  zero_inflation_prob_matrix <- matrix(runif(species * sites, 0.3, 0.8), nrow = species, ncol = sites)  # Ensure some > 0.5

  simulated_data <- simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)

  zero_indices <- which(zero_inflation_prob_matrix > 0.5, arr.ind = TRUE)

  # ✅ Only proceed if zero-inflated sites exist
  if (length(zero_indices) > 0) {
    for (i in 1:replicates) {
      zero_values <- simulated_data[zero_indices[, 1], zero_indices[, 2], i]
      print(zero_values)  # Debugging: Print values
      expect_true(any(zero_values == 0),
                  info = "Zero-inflated sites should contain zero values.")
    }
  } else {
    skip("No zero-inflated sites were generated (zero_inflation_prob_matrix > 0.5 is empty).")
  }
})
