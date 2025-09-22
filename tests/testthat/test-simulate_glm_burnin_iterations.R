test_that("simulate_glm_burnin_iterations returns expected structure", {
  skip_if_not_installed("glmmTMB")
  
  # Simulate a small toy dataset
  set.seed(123)
  mock_data <- data.frame(
    y = rpois(20, lambda = 2),
    i = factor(rep(c("OTU1", "OTU2"), each = 10)),
    treatment = factor(rep(c("A", "B"), 10)),
    Site = factor(rep(c("S1", "S2"), 10)),
    Sample = factor(rep(1:4, length.out = 20)),
    Replicate = factor(rep(1:2, length.out = 20))
  )
  
  # Use minimal formulas for speed
  poisson_formula <- y ~ treatment + (1 | Site)
  binomial_formula <- z_sim ~ treatment + (1 | Site)
  
  result <- simulate_glm_burnin_iterations(
    data_glm = mock_data,
    poisson_formula = poisson_formula,
    binomial_formula = binomial_formula,
    num_iterations = 6,
    burn_in = 3
  )
  
  # --- Checks ---
  expect_type(result, "list")
  expect_named(result, c("poisson_models", "binomial_models"))
  
  expect_length(result$poisson_models, 3)
  expect_length(result$binomial_models, 3)
  
  expect_s3_class(result$poisson_models[[1]], "glmmTMB")
  expect_s3_class(result$binomial_models[[1]], "glmmTMB")
})
