library(testthat)
library(Matrix)


test_that("compute_probabilities returns expected outputs", {
  # ✅ Define test parameters
  species <- 3
  sites <- 3
  replicates <- 3
  lambda <- 5

  # ✅ Generate Zero-Inflated Poisson Data
  zero_inflation_prob_matrix <- matrix(runif(species * sites, 0.1, 0.5), nrow = species, ncol = sites)
  simulated_data <- simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)
  Y <- to2D(simulated_data)

  # ✅ Run TMB model fitting
  fit_results <- prepare_tmb_data_fit(Y)

  # ✅ Compute probabilities
  probabilities <- compute_probabilities(fit_results)

  # ✅ Validate dimensions
  expect_equal(dim(probabilities$occ.prob), dim(fit_results$fit$report(fit_results$fit$env$last.par.best)$etao),
               info = "Mismatch: occ.prob dimensions do not match etao dimensions.")

  expect_equal(nrow(probabilities$prob.detect), length(unique(fit_results$fit$env$data$sites)),
               info = "Mismatch: prob.detect number of rows should match unique sites.")

  # ✅ Ensure probability values are within valid range
  expect_true(all(probabilities$occ.prob >= 0 & probabilities$occ.prob <= 1),
              info = "Error: occ.prob should be between 0 and 1.")

  expect_true(all(probabilities$prob.detect >= 0 & probabilities$prob.detect <= 1),
              info = "Error: prob.detect should be between 0 and 1.")
})
