# File: tests/testthat/test-data_array_phyloseq.R

test_that("data_array_phyloseq correctly creates a 3D array from a phyloseq object", {

  # Mock a small phyloseq object for testing
  library(phyloseq)

  # Create mock OTU table (Species x Samples)
  otu_mat <- matrix(1:12, nrow = 4, dimnames = list(paste0("OTU", 1:4), paste0("Sample", 1:3)))
  otu_table_test <- otu_table(otu_mat, taxa_are_rows = TRUE)

  # Create sample metadata
  sample_metadata <- data.frame(
    sample = paste0("Sample", 1:3),
    samplerep = c("r1", "r2", "r3"),
    sampletype = rep("biologicalsample", 3)
  )
  rownames(sample_metadata) <- sample_metadata$sample
  sample_data_test <- sample_data(sample_metadata)

  # Assemble phyloseq object
  physeq_test <- phyloseq(otu_table_test, sample_data_test)

  # Save as a temporary RDS file
  temp_rds_path <- tempfile(fileext = ".RDS")
  saveRDS(physeq_test, temp_rds_path)

  # Run the function
  result_array <- data_array_phyloseq(temp_rds_path, verbose = FALSE)

  # ✅ Test 1: Check if the output is an array
  expect_true(is.array(result_array), label = "Output should be an array.")

  # ✅ Test 2: Check dimensions (Species x Sites x Replicates)
  expect_equal(length(dim(result_array)), 3, label = "Array should be 3-dimensional.")

  # ✅ Test 3: Check dimension names
  expect_equal(length(dimnames(result_array)$Species), 4, label = "Should have 4 species (OTUs).")
  expect_equal(length(dimnames(result_array)$Sites), 3, label = "Should have 3 sites (Samples).")
  expect_equal(length(dimnames(result_array)$Replicates), 3, label = "Should have 3 replicates.")

  # ✅ Test 4: Check if data is numeric
  expect_true(is.numeric(result_array), label = "Array content should be numeric.")

  # ✅ Test 5: Check if array contains non-zero values (since input was 1:12)
  expect_gt(sum(result_array), 0, label = "Array should contain non-zero values.")

  # ✅ Test 6: Check structure correctness (optional detailed test)
  expect_equal(
    dim(result_array),
    c(4, 3, 3),
    label = "Array dimensions should be (Species x Sites x Replicates) -> (4 x 3 x 3)."
  )

  # ✅ Test 7: Check that OTU names are preserved
  expect_equal(
    dimnames(result_array)$Species,
    paste0("OTU", 1:4),
    label = "Species names should match OTU names."
  )

  # ✅ Clean up temporary RDS file
  unlink(temp_rds_path)

  # ✅ Print confirmation message
  message("✅ Unit test for `data_array_phyloseq` passed all checks.")
})
