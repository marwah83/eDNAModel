library(testthat)

test_that("prepare_tmb_data_fit correctly processes input data", {
  # Define test parameters
  species <- 3
  sites <- 2
  replicates <- 4

  # Create a sample dataset (3D array)
  set.seed(123)
  sample_array <- array(
    data = rpois(species * sites * replicates, lambda = 5),
    dim = c(species, sites, replicates),
    dimnames = list(
      Species = paste("Species", 1:species, sep = "_"),
      Sites = paste("Site", 1:sites, sep = "_"),
      Replicates = paste("Replicate", 1:replicates, sep = "_")
    )
  )

  # Convert 3D array to 2D
  Y <- to2D(sample_array)

  # Run function
  prepared_data <- prepare_tmb_data_fit(Y)

  # ✅ Check if result is a list
  expect_true(is.list(prepared_data), "Output should be a list.")

  # ✅ Check if required keys exist
  expected_keys <- c("y", "ysites", "Xa", "Xo", "Za", "Zo", "sites", "csa", "cso")
  expect_true(all(expected_keys %in% names(prepared_data)), "Missing expected keys in the output list.")

  # ✅ Check dimensions of returned matrices
  expect_equal(dim(prepared_data$y), c(sites * replicates, species), "Incorrect dimensions for y matrix.")
  expect_equal(dim(prepared_data$ysites), c(sites, species), "Incorrect dimensions for ysites matrix.")
  expect_equal(dim(prepared_data$Xa)[1], sites * replicates, "Incorrect row count for Xa.")
  expect_equal(dim(prepared_data$Xo)[1], sites, "Incorrect row count for Xo.")

  # ✅ Check that `sites` matches the expected structure
  expect_length(prepared_data$sites, sites * replicates, "Mismatch in length of sites vector.")

  # ✅ Check if `csa`, `cso`, `Za`, and `Zo` are valid matrices
  expect_true(is.matrix(prepared_data$csa), "csa should be a matrix.")
  expect_true(is.matrix(prepared_data$cso), "cso should be a matrix.")
  expect_s4_class(prepared_data$Za, "dgCMatrix", "Za should be a sparse matrix.")
  expect_s4_class(prepared_data$Zo, "dgCMatrix", "Zo should be a sparse matrix.")

  # ✅ Print debugging information
  print(str(prepared_data))  # Debugging print

})
