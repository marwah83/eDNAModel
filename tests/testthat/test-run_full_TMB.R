# Load required package
library(testthat)

# ---------------------------
# ✅ Final Unit Test for TMB Pipeline Function with Covariates + Formulas
# ---------------------------
test_that("run_full_TMB executes properly with covariates and returns expected structure", {

  # ✅ Mock a simplified 3D data array (Species x Sites x Replicates)
  mock_data_array <- array(
    sample(0:10, 3 * 3 * 4, replace = TRUE),
    dim = c(3, 3, 4),  # 3 species, 3 sites, 4 replicates
    dimnames = list(
      Species = paste0("Species_", 1:3),
      Sites = as.character(1:3),
      Replicates = paste0("r", 1:4)
    )
  )

  # ✅ Mock covariate data
  covariates <- data.frame(
    site = factor(rep(1:3, each = 4)),       # Sites repeated per replicate
    replicate = factor(rep(1:4, times = 3))  # 4 replicates per site
  )

  # ✅ Suppress warnings from compile/link if TMB is not compiled yet
  suppressWarnings({
    result <- run_full_TMB(
      data_array_filtered = mock_data_array,
      covariate_data = covariates,
      o.formula = ~ 1,
      a.formula = ~ 1
    )
  })

  # ---------------------------
  # ✅ Structural Checks
  # ---------------------------
  expect_true(is.list(result), label = "Result should be a list.")

  expect_true(all(c(
    "optimization", "occupancy_probability", "lambda",
    "detection_probability", "check_large_occupancy", "check_small_occupancy"
  ) %in% names(result)),
  label = "Result must contain all expected components.")

  # ---------------------------
  # ✅ Dimension Checks
  # ---------------------------
  expected_species <- length(dimnames(mock_data_array)$Species)       # 3
  expected_sites <- length(dimnames(mock_data_array)$Sites)           # 3
  expected_replicates <- length(dimnames(mock_data_array)$Replicates) # 4
  expected_lambda_rows <- expected_sites * expected_replicates

  expect_equal(nrow(result$lambda), expected_lambda_rows,
               label = "Lambda matrix should have rows = Sites x Replicates.")
  expect_equal(ncol(result$lambda), expected_species,
               label = "Lambda matrix should have columns = Species.")

  expect_equal(nrow(result$occupancy_probability), expected_sites,
               label = "Occupancy matrix should have rows = Sites.")
  expect_equal(ncol(result$occupancy_probability), expected_species,
               label = "Occupancy matrix should have columns = Species.")

  expect_equal(nrow(result$detection_probability), expected_sites,
               label = "Detection matrix should have rows = Sites.")
  expect_equal(ncol(result$detection_probability), expected_species,
               label = "Detection matrix should have columns = Species.")

  # ---------------------------
  # ✅ Logical Output Checks
  # ---------------------------
  expect_type(result$check_large_occupancy, "logical")
  expect_type(result$check_small_occupancy, "logical")

  # ✅ Final confirmation message
  message("✅ Unit test for 'run_full_TMB' with covariates passed all checks!")
})
