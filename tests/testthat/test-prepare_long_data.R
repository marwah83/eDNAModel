test_that("prepare_long_data processes a phyloseq object correctly", {
  skip_if_not_installed("phyloseq")
  
  # Create mock OTU table with sufficient total abundance
  otu <- matrix(
    c(10, 2, 1, 5,   # OTU1 total = 18
      3, 4, 0, 1),   # OTU2 total = 8
    nrow = 2,
    byrow = TRUE
  )
  rownames(otu) <- c("OTU1", "OTU2")
  colnames(otu) <- c("S1_r1", "S2_r1", "S3_r1", "S4_r1")
  otu_table_obj <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  
  # Create matching sample metadata
  sample_data_df <- data.frame(
    sampletype = rep("biologicalsample", 4),
    location = c("Site1", "Site1", "Site2", "Site2"),
    treatment = c("A", "A", "B", "B"),
    row.names = colnames(otu)
  )
  sample_data_obj <- phyloseq::sample_data(sample_data_df)
  
  # Build phyloseq object
  ps <- phyloseq::phyloseq(otu_table_obj, sample_data_obj)
  
  # Run function
  result <- prepare_long_data(
    physeq_obj = ps,
    min_species_sum = 1,
    sampletype_var = "sampletype",
    sampletype_keep = "biologicalsample",
    location_var = "location",
    treatment_var = "treatment"
  )
  
  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("physeq_filtered", "long_df"))
  
  # Check long_df structure
  expect_s3_class(result$long_df, "data.frame")
  expect_true(all(c("i", "Site", "Sample", "Replicate", "treatment", "y") %in% names(result$long_df)))
  
  # Check column types
  expect_true(is.factor(result$long_df$i))
  expect_true(is.integer(result$long_df$y))
  expect_true(is.factor(result$long_df$Sample))
  expect_true(is.factor(result$long_df$Replicate))
  expect_true(is.factor(result$long_df$Site))
})
