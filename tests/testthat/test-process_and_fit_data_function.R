# ✅ Test process_and_fit_data function
test_that("process_and_fit_data returns non-null fit results", {
  species <- 5
  sites <- 5
  replicates <- 3
  lambda <- 2
  ZIP <- 0.5
  zero_inflation_matrix <- create_zero_inflation_matrix(species, sites, ZIP)

  simulated_data <- simulate_data(species, sites, replicates, lambda, zero_inflation_matrix)
  fit_results <- process_and_fit_data(simulated_data)

  expect_true(!is.null(fit_results))
})

