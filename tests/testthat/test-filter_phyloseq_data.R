test_that("filter_phyloseq_data filters low-quality samples and taxa correctly", {
  skip_if_not_installed("phyloseq")
  
  # Create mock OTU table with 2 taxa and 4 samples
  otu <- matrix(
    c(10, 0, 0, 0,
      5,  3, 1, 2),
    nrow = 2,
    byrow = TRUE
  )
  rownames(otu) <- c("OTU1", "OTU2")
  colnames(otu) <- c("S1", "S2", "S3", "S4")
  
  otu_table_obj <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  
  # Sample metadata
  sample_data_df <- data.frame(
    dummy = 1:4,
    row.names = colnames(otu)
  )
  sample_data_obj <- phyloseq::sample_data(sample_data_df)
  
  # Assemble phyloseq object
  physeq <- phyloseq::phyloseq(otu_table_obj, sample_data_obj)
  
  # Run filtering
  filtered <- filter_phyloseq_data(physeq, min_species_sum = 10, save_path = NULL)
  
  # Basic checks
  expect_s4_class(filtered, "phyloseq")
  expect_gt(phyloseq::nsamples(filtered), 0)
  expect_gt(phyloseq::ntaxa(filtered), 0)
})

test_that("filter_phyloseq_data errors on invalid input", {
  expect_error(
    filter_phyloseq_data("not a phyloseq object"),
    regexp = "must be a valid `phyloseq` object"
  )
})

test_that("filter_phyloseq_data warns when no species meet min threshold", {
  skip_if_not_installed("phyloseq")
  
  # OTU matrix with low counts
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
  
  expect_warning(
    filter_phyloseq_data(ps, min_species_sum = 10, save_path = NULL),
    regexp = "No species meet the min_species_sum"
  )
})
