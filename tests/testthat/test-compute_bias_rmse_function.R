# ✅ Test compute_bias_rmse function
test_that("compute_bias_rmse handles matching dimensions", {
  species <- 5
  sites <- 5
  lambda <- 2
  ZIP <- 0.5
  true_occ_prob <- create_zero_inflation_matrix(species, sites, ZIP)
  est_occ_prob <- true_occ_prob + matrix(runif(species * sites, -0.05, 0.05), nrow = species)

  true_prob_detect <- 1 - true_occ_prob * exp(-lambda)
  est_prob_detect <- true_prob_detect + matrix(runif(species * sites, -0.05, 0.05), nrow = species)

  results <- compute_bias_rmse(true_occ_prob, est_occ_prob, true_prob_detect, est_prob_detect)

  expect_true(!is.null(results$bias_occ))
  expect_true(!is.null(results$rmse_occ))
  expect_true(!is.null(results$bias_detect))
  expect_true(!is.null(results$rmse_detect))
})

