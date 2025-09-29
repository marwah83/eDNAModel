test_that("compute_rmse calculates correctly", {
  actual <- c(1, 2, 3, 4, 5)
  estimated <- c(1.1, 2.2, 2.9, 4.1, 5.0)
  expected_rmse <- sqrt(mean((actual - estimated) ^ 2))
  expect_equal(compute_rmse(actual, estimated), expected_rmse, tolerance = 1e-6)
})
