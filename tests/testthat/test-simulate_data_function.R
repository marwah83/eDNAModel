# ✅ Test simulate_data function
test_that("simulate_data generates non-empty data", {
  species <- 5
  sites <- 5
  replicates <- 3
  lambda <- 2
  ZIP <- 0.5
  zero_inflation_matrix <- create_zero_inflation_matrix(species, sites, ZIP)

  simulated_data <- simulate_data(species, sites, replicates, lambda, zero_inflation_matrix)
  expect_true(!is.null(simulated_data))
  expect_equal(dim(simulated_data), c(species, sites, replicates))
})

