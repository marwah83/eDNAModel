# ✅ Test run_full_simulation function
test_that("run_full_simulation completes successfully", {
  species <- 5
  sites <- 5
  replicates <- 3
  lambda <- 2
  ZIP <- 0.5

  results <- run_full_simulation(species, sites, replicates, lambda, ZIP)

  expect_true(!is.null(results))
  expect_equal(ncol(results), 12)
  expect_true(all(!is.na(results$Bias_Occupancy)))
})

