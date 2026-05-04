test_that("prepare_long_data returns correctly formatted output", {

  set.seed(123)

  library(phyloseq)
  library(dplyr)

  # -----------------------------
  # Simulate metadata
  # -----------------------------
  sites <- paste0("Site", 1:2)
  otus <- paste0("OTU", 1:5)
  treatments <- c("Control", "Rats")
  reps <- 1:2

  meta_data <- expand.grid(
    site_field  = sites,
    treat_field = treatments,
    rep_field   = reps
  )

  meta_data$unique_id <- paste0(
    meta_data$site_field, "_",
    meta_data$treat_field, "_r",
    meta_data$rep_field
  )

  rownames(meta_data) <- meta_data$unique_id

  # -----------------------------
  # Simulate OTU table
  # -----------------------------
  otu_mat <- matrix(
    rpois(n = length(otus) * nrow(meta_data), lambda = 10),
    nrow = length(otus),
    ncol = nrow(meta_data),
    dimnames = list(otus, meta_data$unique_id)
  )

  physeq <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    sample_data(meta_data)
  )

  # -----------------------------
  # Run function
  # -----------------------------
  result <- prepare_long_data(
    physeq_obj = physeq,
    site_col = "site_field",
    nested_cols = c("site_field", "treat_field", "rep_field")
  )

  # -----------------------------
  # Structural checks
  # -----------------------------
  expect_type(result, "list")
  expect_named(result, c("physeq", "long_df"), ignore.order = TRUE)
  expect_s4_class(result$physeq, "phyloseq")

  long_df <- result$long_df
  expect_s3_class(long_df, "data.frame")

  # -----------------------------
  # Column checks
  # -----------------------------
  expect_true(all(c("SampleRep", "OTU", "y", "Site") %in% names(long_df)))

  # -----------------------------
  # Type checks
  # -----------------------------
  expect_s3_class(long_df$OTU, "factor")
  expect_true(is.numeric(long_df$y))

  # -----------------------------
  # SampleRep correctness (MATCH FUNCTION)
  # -----------------------------
  # IMPORTANT: SampleRep = rownames(meta_data)
  expected_SampleRep <- rownames(meta_data)

  expect_setequal(
    unique(long_df$SampleRep),
    expected_SampleRep
  )

  # -----------------------------
  # Site column checks
  # -----------------------------
  expect_true(all(!is.na(long_df$Site)))
  expect_true(all(long_df$Site == long_df$site_field))

  # -----------------------------
  # Dimension checks
  # -----------------------------
  expect_equal(
    nrow(long_df),
    length(otus) * nrow(meta_data)
  )

  # -----------------------------
  # No filtering check
  # -----------------------------
  expect_equal(
    length(unique(long_df$OTU)),
    length(otus)
  )

  # -----------------------------
  # Metadata preservation
  # -----------------------------
  expect_true(all(c("treat_field", "rep_field") %in% names(long_df)))

  # -----------------------------
  # Sanity checks
  # -----------------------------
  expect_true(all(long_df$y >= 0))
  expect_true(all(is.finite(long_df$y)))

})
