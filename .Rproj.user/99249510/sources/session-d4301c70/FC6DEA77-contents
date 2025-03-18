library(testthat)
library(eDNAmodel)  # Ensure the package is loaded

test_that("TMB pipeline runs correctly", {
  # Simulate input data
  species <- 5
  sites <- 4
  replicates <- 3
  lambda <- 5

  set.seed(123)  # For reproducibility
  zero_inflation_prob_matrix <- matrix(runif(species * sites, 0.1, 0.5), nrow = species, ncol = sites)

  # Step 1: Simulate ZIP data
  simulated_data <- simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)

  # Step 2: Convert to 2D
  Y <- to2D(simulated_data)

  # Step 3: Prepare data for TMB
  data_list <- prepare_tmb_data_fit(Y)

  # Step 4: Ensure all required elements exist
  expect_true(is.matrix(data_list$y))
  expect_true(is.matrix(data_list$ysites))
  expect_true(is.matrix(data_list$Xa))
  expect_true(is(data_list$Za, "TsparseMatrix"))

  # Step 5: Run the TMB model
  fit_results <- run_TMB_model(data_list)

  # Step 6: Check optimization output
  expect_true(!is.null(fit_results$fit))
  expect_true(!is.null(fit_results$opt))
  expect_true(length(fit_results$opt$par) > 0)

  # Step 7: Compute probabilities
  plot_results <- plot_probabilities(fit_results$fit, Y)

  # Step 8: Validate output probabilities
  expect_true(all(plot_results$occ.prob >= 0 & plot_results$occ.prob <= 1))  # Occupancy prob must be in [0,1]
  expect_true(all(plot_results$prob.detect >= 0 & plot_results$prob.detect <= 1))  # Detection prob must be in [0,1]

  # Step 9: Ensure no missing values in final results
  expect_true(!any(is.na(plot_results$occ.prob)))
  expect_true(!any(is.na(plot_results$prob.detect)))
})

