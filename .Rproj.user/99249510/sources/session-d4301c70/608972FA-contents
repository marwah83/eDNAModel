test_that("simulate_zip_data correctly applies zero-inflation", {
  species <- 3
  sites <- 3
  replicates <- 2
  lambda <- 5

  # Create a zero-inflation matrix with high probability
  zero_inflation_prob_matrix <- matrix(runif(species * sites, 0.6, 0.9),
                                       nrow = species, ncol = sites)

  # Simulate data
  result <- simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)

  # Identify indices where zero-inflation should have been applied
  zero_indices <- which(zero_inflation_prob_matrix > 0.5, arr.ind = TRUE)

  # Ensure that at least some values at these locations are zero
  for (i in 1:replicates) {
    zero_values <- result[zero_indices[, 1], zero_indices[, 2], i]
    print(zero_values)  # Debugging: Print values
    expect_true(any(zero_values == 0))
  }
})
