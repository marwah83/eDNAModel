# ✅ Test compute_probabilities_wrapper function
test_that("compute_probabilities_wrapper returns valid probabilities", {
  species <- 5
  sites <- 5
  replicates <- 3
  lambda <- 2
  ZIP <- 0.5
  zero_inflation_matrix <- create_zero_inflation_matrix(species, sites, ZIP)

  simulated_data <- simulate_data(species, sites, replicates, lambda, zero_inflation_matrix)
  fit_results <- process_and_fit_data(simulated_data)
  probabilities <- compute_probabilities_wrapper(fit_results)

  expect_true(!is.null(probabilities))
  expect_true(all(probabilities$occ.prob >= 0 & probabilities$occ.prob <= 1))
  expect_true(all(probabilities$prob.detect >= 0 & probabilities$prob.detect <= 1))
})

