test_that("prepare_long_data errors when all OTUs are filtered out", {
  skip_if_not_installed("phyloseq")
  
  # OTU matrix with all zeros
  otu <- matrix(c(0, 0, 0, 0), nrow = 1)
  rownames(otu) <- "OTU1"
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
  
  # Function should error due to all OTUs being filtered out
  expect_error(
    prepare_long_data(
      physeq_obj = ps,
      min_species_sum = 5,
      sampletype_var = "sampletype",
      sampletype_keep = "biologicalsample",
      location_var = "location",
      treatment_var = "treatment"
    ),
    regexp = "No species left after removing all-zero species"  # ✔ Match exact message
  )
})
