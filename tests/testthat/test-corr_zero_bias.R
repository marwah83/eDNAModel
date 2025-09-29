# ✅ Test correlation between zero proportion and bias
test_that("correlation between zero proportion and bias is computed correctly", {
  species <- 5
  sites <- 5
  replicates <- 3
  lambda <- 2
  ZIP <- 0.5
  results <- run_full_simulation(species, sites, replicates, lambda, ZIP)

  valid_rows <- complete.cases(results$zero,results$Bias_Occupancy)
  correlation <- cor(results$zero[valid_rows], results$Bias_Occupancy[valid_rows])
  expect_true(!is.na(correlation))
  expect_true(correlation >= -1 & correlation <= 1)
})



