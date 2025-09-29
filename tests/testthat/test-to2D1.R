test_that("to2D correctly transforms 3D array into a data frame", {
  species <- 3
  sites <- 3
  replicates <- 3
  lambda <- 5
  zero_inflation_prob_matrix <- matrix(runif(species * sites, 0.1, 0.5), nrow = species, ncol = sites)

  simulated_data <- simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)
  Y <- to2D(simulated_data)

  # Check structure
  expect_true(is.data.frame(Y), info = "Y should be a data frame")

  # Check number of rows
  expect_equal(nrow(Y), sites * replicates, info = "Y should have sites * replicates rows")

  # Check column names
  expect_equal(colnames(Y), c("Site", "Replicate", paste0("Species_", 1:species)),
               info = "Column names should match expected format")
})
