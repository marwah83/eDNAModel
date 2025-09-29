test_that("prepare_tmb_data_fit correctly formats data for TMB", {
  species <- 3
  sites <- 3
  replicates <- 3
  lambda <- 5
  zero_inflation_prob_matrix <- matrix(runif(species * sites, 0.1, 0.5), nrow = species, ncol = sites)
  simulated_data <- simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)
  Y <- to2D(simulated_data)
  fit_data <- prepare_tmb_data_fit(Y)

  # Ensure expected elements are present
  expect_true("fit" %in% names(fit_data), info = "Output should contain 'fit'")
  expect_true("opt" %in% names(fit_data), info = "Output should contain 'opt'")

  # Ensure Xa and Xo matrices are correct
  expect_true(is.matrix(fit_data$fit$env$data$Xa), info = "Xa should be a matrix")
  expect_true(is.matrix(fit_data$fit$env$data$Xo), info = "Xo should be a matrix")

  # Check site structure
  expect_equal(length(fit_data$fit$env$data$sites), nrow(Y),
               info = "Sites vector should match number of rows in Y")
})
