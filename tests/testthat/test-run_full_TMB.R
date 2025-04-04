library(testthat)

# ---------------------------
# ✅ Unit Test for `run_full_TMB1`
# ---------------------------
test_that("run_full_TMB executes and returns expected structure", {

  # ✅ Mock 3D data array (Species x Sites x Replicates)
  mock_data_array <- array(
    sample(0:10, 2 * 3 * 4, replace = TRUE),
    dim = c(2, 3, 4),
    dimnames = list(
      Species = paste0("Species_", 1:2),
      Sites = as.character(1:3),
      Replicates = paste0("r", 1:4)
    )
  )

  # ✅ Create mock covariate data (site and replicate)
  covariate_data <- data.frame(
    site = factor(rep(1:3, each = 4)),
    replicate = factor(rep(1:4, times = 3))
  )

  # ✅ Suppress warnings from compilation or TMB internals
  result <- suppressWarnings(
    run_full_TMB(
      data_array_filtered = mock_data_array,
      covariate_data = covariate_data,
      a.formula = ~ site + (1 | replicate),
      o.formula = ~ site
    )
  )

  # ✅ Check result is a list
  expect_true(is.list(result), label = "Result should be a list")

  # ✅ Check result has expected names
  expect_true(all(c("optimization", "occupancy_probability",
                    "detection_probability") %in% names(result)),
              label = "List must contain required outputs")

  # ✅ Check expected dimensions
  expect_true(is.matrix(result$occupancy_probability), label = "Occupancy prob must be matrix")
  expect_true(is.data.frame(result$detection_probability), label = "Detection prob must be dataframe")
})

