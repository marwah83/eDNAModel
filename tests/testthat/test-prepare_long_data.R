test_that("prepare_long_data returns correctly formatted output", {
  set.seed(123)

  # -----------------------------
  # Simulate metadata
  # -----------------------------
  sites <- paste0("Site", 1:2)
  otus <- paste0("OTU", 1:5)
  treatments <- c("Control", "Rats")
  reps <- 1:2

  meta_data <- expand.grid(
    site_field = sites,                  # simulate site
    treat_field = treatments,            # simulate treatment
    rep_field = reps                     # simulate replicate
  )

  meta_data$unique_id <- paste0(meta_data$site_field, "_", meta_data$treat_field, "_r", meta_data$rep_field)
  meta_data$sampletype <- "biologicalsample"
  meta_data$location <- meta_data$site_field

  rownames(meta_data) <- meta_data$unique_id

  # -----------------------------
  # Simulate OTU table (5 OTUs × N samples)
  # -----------------------------
  otu_mat <- matrix(
    rpois(n = 5 * nrow(meta_data), lambda = 10),
    nrow = 5,
    ncol = nrow(meta_data),
    dimnames = list(otus, meta_data$unique_id)
  )

  otu <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  sample_tab <- phyloseq::sample_data(meta_data)
  physeq <- phyloseq::phyloseq(otu, sample_tab)

  # -----------------------------
  # Run the function
  # -----------------------------
  result <- prepare_long_data(
    physeq_obj = physeq,
    min_species_sum = 10,
    site_col = "site_field"
  )

  # -----------------------------
  # Structural checks
  # -----------------------------
  expect_type(result, "list")
  expect_named(result, c("physeq_filtered", "long_df"))

  # Check physeq_filtered is a phyloseq object
  expect_s4_class(result$physeq_filtered, "phyloseq")

  # Check long_df structure
  long_df <- result$long_df
  expect_s3_class(long_df, "data.frame")

  expect_true(all(c("SampleRep", "OTU", "y", "Site") %in% names(long_df)))

  expect_type(long_df$y, "integer")
  expect_s3_class(long_df$OTU, "factor")
  expect_s3_class(long_df$Site, "factor")
})
