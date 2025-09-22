test_that("extract_fixed extracts coefficients correctly", {
  skip_if_not_installed("glmmTMB")
  
  # Simulate minimal dataset
  set.seed(42)
  data <- data.frame(
    y = rpois(30, lambda = 3),
    treatment = factor(rep(c("A", "B"), each = 15)),
    site = factor(rep(c("S1", "S2", "S3"), 10))
  )
  
  # Fit multiple glmmTMB models
  model_list <- lapply(1:2, function(i) {
    glmmTMB::glmmTMB(y ~ treatment + (1 | site), family = poisson, data = data)
  })
  
  # Run extract_fixed
  fixed_df <- extract_fixed(model_list, model_name = "poisson")
  
  # Check structure
  expect_s3_class(fixed_df, "data.frame")
  expect_true(all(c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "term", "iter", "model") %in% colnames(fixed_df)))
  
  # Expect 2 models × (intercept + treatmentB)
  expect_equal(nrow(fixed_df), 2 * 2)
  
  # Check model name column
  expect_true(all(fixed_df$model == "poisson"))
  
  # Check that each iteration is labeled
  expect_true(all(fixed_df$iter %in% 1:2))
})
