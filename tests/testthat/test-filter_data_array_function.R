# Load required library
library(testthat)

# Unit test for filter_data_array function
test_that("filter_data_array correctly filters 3D array based on sum thresholds", {

  # ---- Mock Data Setup ----
  # Create a mock 3D array (Species x Sites x Replicates)
  set.seed(123)
  mock_array <- array(
    sample(0:10, 3 * 4 * 2, replace = TRUE),  # 3 species, 4 sites, 2 replicates
    dim = c(3, 4, 2),
    dimnames = list(
      Species = paste0("Sp", 1:3),
      Sites = paste0("Site", 1:4),
      Replicates = paste0("r", 1:2)
    )
  )

  # ---- Set Temporary File Path ----
  temp_file <- tempfile(fileext = ".Rdata")  # create temp file path

  # ---- Apply Function ----
  filtered_array <- filter_data_array(mock_array, min_species_sum = 10, save_path = temp_file)  # use temp file path

  # ---- Test Structure ----
  expect_true(is.array(filtered_array))                # Check if output is still an array
  expect_equal(length(dim(filtered_array)), 3)        # Check that output remains 3D
  expect_true(all(dim(filtered_array) <= dim(mock_array)))  # Check that dimensions did not increase

  # ---- Test Filtering Logic ----
  species_sums <- apply(filtered_array, 1, sum)  # Sum over sites and replicates for each species
  expect_true(all(species_sums >= 10))           # Ensure all species have sum >= 10 (min_species_sum)

  site_sums <- apply(filtered_array, 2, sum)     # Sum over species and replicates for each site
  expect_true(all(site_sums > 0))                # Ensure no site is all-zero

  # ---- Test Content ----
  expect_gt(sum(filtered_array), 0)              # Ensure the filtered array has some non-zero values

  # ---- Test Save Path ----
  expect_true(file.exists(temp_file))           # ✅ Check that the file was saved correctly

})
