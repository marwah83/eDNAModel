test_that("prepare_long_data processes a phyloseq object correctly", {
  skip_if_not_installed("phyloseq")
  
  # Mock OTU matrix with sufficient abundance
  otu <- matrix(
    c(10, 5, 5, 10,
      3, 4, 6, 7),  # each OTU appears in all 4 samples
    nrow = 2,
    byrow = TRUE
  )
  rownames(otu) <- c("OTU1", "OTU2")
  colnames(otu) <- c("S1_r1", "S2_r1", "S3_r1", "S4_r1")
  otu_table_obj <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  
  # Sample metadata
  sample_data_df <- data.frame(
    sampletype = rep("biologicalsample", 4),
    location = c("Site1", "Site1", "Site2", "Site2"),
    treatment = c("A", "A", "B", "B"),
    row.names = colnames(otu)
  )
  sample_data_obj <- phyloseq::sample_data(sample_data_df)
  
  ps <- phyloseq::phyloseq(otu_table_obj, sample_data_obj)
  
  # Run the function
  result <- prepare_long_data(
    physeq_obj = ps,
    min_species_sum = 10,  # Lower this to keep some OTUs
    sampletype_var = "sampletype",
    sampletype_keep = "biologicalsample",
    location_var = "location",
    treatment_var = "treatment"
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("physeq_filtered", "long_df"))
  
  # Check long_df structure
  expect_s3_class(result$long_df, "data.frame")
  expect_true(all(c("i", "Site", "Sample", "Replicate", "treatment", "y") %in% names(result$long_df)))
  expect_true(is.factor(result$long_df$i))
  expect_true(is.integer(result$long_df$y))
})
