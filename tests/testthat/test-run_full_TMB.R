# Load required package
library(testthat)

# ---------------------------
# ✅ Final Unit Test for TMB Pipeline Function
# ---------------------------
test_that("run_full_TMB_pipeline_with_results executes properly and returns correct output structure", {

  # ✅ Mock a simplified 3D data array (Species x Sites x Replicates)
  mock_data_array <- array(
    sample(0:10, 3 * 3 * 4, replace = TRUE),  # Random counts for simplicity
    dim = c(3, 3, 4),  # 3 species, 3 sites, 4 replicates
    dimnames = list(
      Species = paste0("Species_", 1:3),
      Sites = as.character(1:3),  # Numeric-like characters as site IDs
      Replicates = paste0("r", 1:4)
    )
  )

  # ✅ Suppress compile/read warnings and run the function
  suppressWarnings({
    result <- run_full_TMB(mock_data_array)
  })

  # ---------------------------
  # ✅ Structural Checks
  # ---------------------------

  # ✅ Check output type
  expect_true(is.list(result), label = "Result should be a list.")

  # ✅ Check required components are present
  expect_true(all(c(
    "optimization", "occupancy_probability", "lambda",
    "detection_probability", "check_large_occupancy", "check_small_occupancy"
  ) %in% names(result)),
  label = "Result must contain all expected components."
  )

  # ---------------------------
  # ✅ Dimension Checks
  # ---------------------------
  expected_species <- length(dimnames(mock_data_array)$Species)       # 3 species
  expected_sites <- length(dimnames(mock_data_array)$Sites)           # 3 sites
  expected_replicates <- length(dimnames(mock_data_array)$Replicates) # 4 replicates
  expected_lambda_rows <- expected_sites * expected_replicates        # Sites x Replicates

  # ✅ Lambda matrix (Sites x Replicates) x Species
  expect_equal(nrow(result$lambda), expected_lambda_rows,
               label = "Lambda matrix should have rows equal to Sites x Replicates.")
  expect_equal(ncol(result$lambda), expected_species,
               label = "Lambda matrix should have columns equal to number of species.")

  # ✅ Occupancy probability matrix (Sites x Species)
  expect_equal(nrow(result$occupancy_probability), expected_sites,
               label = "Occupancy probability matrix should have rows equal to number of Sites.")
  expect_equal(ncol(result$occupancy_probability), expected_species,
               label = "Occupancy probability matrix should have columns equal to number of species.")

  # ✅ Detection probability matrix (Sites x Species) -- as per clarified final structure
  expect_equal(nrow(result$detection_probability), expected_sites,
               label = "Detection probability matrix should have rows equal to number of Sites.")
  expect_equal(ncol(result$detection_probability), expected_species,
               label = "Detection probability matrix should have columns equal to number of species.")

  # ---------------------------
  # ✅ Logical Output Checks
  # ---------------------------
  expect_type(result$check_large_occupancy, "logical")
  expect_type(result$check_small_occupancy, "logical")

  # ✅ Final confirmation message
  message("✅ Unit test for 'run_full_TMB' passed all structural checks!")
})
