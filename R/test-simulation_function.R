library(testthat)

# Example of smaller grid for unit testing
param_grid_test <- expand.grid(
  Lambda = c(0.5, 1),
  ZIP = 0.5,
  Species = c(5, 10),
  Sites = c(5, 10),
  Replicates = 2
)

test_that("run_full_simulation runs and returns correct structure", {

  # Start time
  start_time <- Sys.time()

  # Run the function on test grid
  simulation_results_test <- do.call(rbind, lapply(1:nrow(param_grid_test), function(i) {
    with(param_grid_test[i, ], run_full_simulation(Species, Sites, Replicates, Lambda, ZIP))
  }))

  end_time <- Sys.time()

  # ✅ Test 1: Check output is a data frame
  expect_s3_class(simulation_results_test, "data.frame")

  # ✅ Test 2: Check expected columns exist
  expected_columns <- c(
    "Species", "Sites", "Replicates", "Lambda", "ZIP",
    "Mean_Occ_Prob", "Mean_Prob_Detect",
    "Bias_Occupancy", "RMSE_Occupancy",
    "Bias_Detection", "RMSE_Detection",
    "zero"
  )
  expect_true(all(expected_columns %in% colnames(simulation_results_test)))

  # ✅ Test 3: Check number of rows equals parameter grid
  expect_equal(nrow(simulation_results_test), nrow(param_grid_test))

  # ✅ Test 4: Check that numeric columns are numeric
  numeric_columns <- c(
    "Mean_Occ_Prob", "Mean_Prob_Detect",
    "Bias_Occupancy", "RMSE_Occupancy",
    "Bias_Detection", "RMSE_Detection", "zero"
  )
  expect_true(all(sapply(simulation_results_test[numeric_columns], is.numeric)))

  # ✅ Test 5: Check no row is completely NA (at least Lambda should not be NA)
  expect_false(any(is.na(simulation_results_test$Lambda)))

  # ✅ Test 6: Optional - Check if the computation is within a reasonable time for small test
  expect_true(as.numeric(difftime(end_time, start_time, units = "secs")) < 60)  # e.g., less than 60 seconds for small grid

  # ✅ Print result for inspection (optional in testing framework)
  print("✅ Unit test passed for simulation grid execution!")
  print(head(simulation_results_test))

})
