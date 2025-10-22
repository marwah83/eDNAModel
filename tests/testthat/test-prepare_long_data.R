test_that("prepare_long_data returns correctly formatted output", {
  set.seed(123)

  # Simulate metadata
  sites <- paste0("Site", 1:2)
  otus <- paste0("OTU", 1:5)
  treatments <- c("Control", "Rats")
  replicates <- 1:2

  meta_data <- expand.grid(
    Site = sites,
    treatment = treatments,
    Replicate = replicates
  )

  meta_data$SampleID <- paste0(meta_data$Site, "_", meta_data$treatment, "_r", meta_data$Replicate)
  meta_data$sampletype <- "biologicalsample"
  meta_data$location <- meta_data$Site

  rownames(meta_data) <- meta_data$SampleID

  # Simulate OTU count matrix: 5 OTUs × n samples
  otu_mat <- matrix(
    data = rpois(n = 5 * nrow(meta_data), lambda = 10),
    nrow = 5,
    ncol = nrow(meta_data),
    dimnames = list(otus, meta_data$SampleID)
  )

  # Construct phyloseq object
  otu <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  sample_df <- phyloseq::sample_data(meta_data)
  physeq <- phyloseq::phyloseq(otu, sample_df)

  # Run the function
  result <- prepare_long_data(physeq_obj = physeq, min_species_sum = 10, sampletype_keep = "biologicalsample")

  # Validate structure
  expect_type(result, "list")
  expect_named(result, c("physeq_filtered", "long_df"))

  # Check filtered phyloseq object
  expect_s4_class(result$physeq_filtered, "phyloseq")

  # Check long_df
  long_df <- result$long_df
  expect_s3_class(long_df, "data.frame")
  expect_true(all(c("Site", "OTU", "SampleRep", "Sample", "Replicate", "y") %in% names(long_df)))
  expect_type(long_df$y, "integer")
  expect_s3_class(long_df$OTU, "factor")
  expect_s3_class(long_df$Site, "factor")
})
