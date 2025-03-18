library(testthat)

# --- Sample Mock Data Generator for to2D Testing ---
create_mock_array <- function(species = 3, sites = 4, replicates = 2) {
  array(
    rep(1:0, length.out = species * sites * replicates),  # Alternating 1/0 pattern
    dim = c(species, sites, replicates),
    dimnames = list(
      Species = paste0("Sp", 1:species),
      Sites = paste0("S", 1:sites),  # "S1", "S2", ...
      Replicates = paste0("R", 1:replicates)
    )
  )
}

# --- Unit Test for to2D ---
test_that("to2D correctly converts 3D array to 2D data frame", {

  # Create mock data
  mock_array <- create_mock_array()

  # Run to2D
  result <- to2D(mock_array)

  # ✅ Check if result is a data.frame
  expect_true(is.data.frame(result), info = "Output should be a data.frame")

  # ✅ Check number of rows (sites * replicates)
  expected_rows <- dim(mock_array)[2] * dim(mock_array)[3]
  expect_equal(nrow(result), expected_rows, info = "Number of rows should equal sites * replicates")

  # ✅ Check number of columns (Site + Replicate + species)
  expected_cols <- 2 + dim(mock_array)[1]
  expect_equal(ncol(result), expected_cols, info = "Number of columns should equal 2 + number of species")

  # ✅ Check Site column is numeric
  expect_true(is.numeric(result$Site), info = "Site column should be numeric")

  # ✅ Check Replicate column is character (or factor)
  expect_true(is.character(result$Replicate) || is.factor(result$Replicate),
              info = "Replicate column should be character or factor")

  # ✅ Check that species columns are named correctly
  expected_species_names <- dimnames(mock_array)$Species
  expect_equal(colnames(result)[3:ncol(result)], expected_species_names,
               info = "Species columns should be named as in dimnames")

  # ✅ Optional: Print result for visual inspection (during development)
  print(head(result))

})

# --- Unit Test for missing dimnames handling ---
test_that("to2D throws error if dimnames are missing", {
  # Create array without dimnames
  malformed_array <- array(
    rep(1:0, length.out = 3 * 4 * 2),
    dim = c(3, 4, 2)
  )

  # Expect error when running to2D on malformed array
  expect_error(
    to2D(malformed_array),
    regexp = "must have dimnames",
    info = "Should throw error when dimnames are missing"
  )
})
