test_that("filter_phyloseq_data errors when no species meet min threshold", {
  skip_if_not_installed("phyloseq")
  
  # Mock input: all species too rare
  otu <- matrix(c(1, 1, 1, 1), nrow = 1)
  rownames(otu) <- "OTU1"
  colnames(otu) <- c("S1", "S2", "S3", "S4")
  otu_table_obj <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  
  sample_data_df <- data.frame(
    dummy = 1:4,
    row.names = colnames(otu)
  )
  sample_data_obj <- phyloseq::sample_data(sample_data_df)
  
  ps <- phyloseq::phyloseq(otu_table_obj, sample_data_obj)
  
  expect_error(
    filter_phyloseq_data(ps, min_species_sum = 10, save_path = NULL),
    regexp = "No species meet the min_species_sum"
  )
})
