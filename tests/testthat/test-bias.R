test_that("compute_bias calculates correctly", {
  actual <- c(1, 2, 3, 4, 5)
  estimated <- c(1.1, 2.2, 2.9, 4.1, 5.0)
  expected_bias <- mean(estimated - actual)
  expect_equal(compute_bias(actual, estimated), expected_bias, tolerance = 1e-6)
})
