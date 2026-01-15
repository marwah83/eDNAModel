test_that("FitModel_gllvm runs and returns expected structure without noisy warnings", {
  skip_if_not_installed("gllvm")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("phyloseq")
  
  # Check test data exists
  expect_true(exists("physeq_one"), info = "physeq_one is not available")
  expect_true("Name" %in% colnames(sample_data(physeq_one)), info = "Column 'Name' missing from sample_data")

  # Reduce complexity: fewer levels in covariates (if needed)
  physeq_test <- physeq_one
  sample_data(physeq_test)$Samplingmonth <- forcats::fct_lump_n(
    sample_data(physeq_test)$Samplingmonth, n = 3
  )
  
  # Run the model, suppressing expected warnings
  out <- suppressWarnings({
    FitModel_gllvm(
      phyloseq = physeq_test,
      site_col = "Name",
      abundance_rhs = (1 | Samplingmonth / OTU) + (1 | OTU),
      occupancy_covars = "Samplingmonth",
      abundance_family = "poisson",
      n_iter = 3,
      burn_in = 1
    )
  })

  # Check result
  expect_type(out, "list")
  
  expected_components <- c(
    "summary", "psi_list", "lambda_list", "p_detect_list",
    "occupancy_models", "abundance_models", "reduced_data",
    "lv_sites", "lv_species", "mean_lv_sites", "mean_lv_species"
  )
  expect_true(all(expected_components %in% names(out)))
  
  # Additional content checks
  expect_gt(nrow(out$summary), 0)
  expect_s3_class(out$lv_sites, "data.frame")
  expect_s3_class(out$mean_lv_species, "data.frame")
})
