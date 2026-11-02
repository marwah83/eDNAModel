
test_that("prepare_long_data returns correctly formatted output", {
  
  set.seed(123)
  
  # --------------------------------------------------
  # 1. Simulate metadata
  # --------------------------------------------------
  
  sites <- paste0("Site", 1:2)
  otus <- paste0("OTU", 1:5)
  treatments <- c("Control", "Rats")
  reps <- 1:2
  
  meta_data <- expand.grid(
    site_field  = sites,
    treat_field = treatments,
    rep_field   = reps,
    stringsAsFactors = FALSE
  )
  
  meta_data$unique_id <- paste0(
    meta_data$site_field, "_",
    meta_data$treat_field, "_r",
    meta_data$rep_field
  )
  
  rownames(meta_data) <- meta_data$unique_id
  
  # --------------------------------------------------
  # 2. Simulate OTU table
  # --------------------------------------------------
  
  otu_mat <- matrix(
    rpois(
      n = length(otus) * nrow(meta_data),
      lambda = 10
    ),
    nrow = length(otus),
    ncol = nrow(meta_data),
    dimnames = list(
      otus,
      meta_data$unique_id
    )
  )
  
  physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(
      otu_mat,
      taxa_are_rows = TRUE
    ),
    phyloseq::sample_data(meta_data)
  )
  
  # --------------------------------------------------
  # 3. Run prepare_long_data()
  # --------------------------------------------------
  
  result <- prepare_long_data(
    physeq_obj = physeq,
    site_col = "site_field",
    nested_cols = c(
      "site_field",
      "treat_field",
      "rep_field"
    )
  )
  
  # --------------------------------------------------
  # 4. Structural checks
  # --------------------------------------------------
  
  expect_type(result, "list")
  
  expect_named(
    result,
    c("physeq", "long_df"),
    ignore.order = TRUE
  )
  
  expect_s4_class(result$physeq, "phyloseq")
  
  long_df <- result$long_df
  
  expect_s3_class(long_df, "data.frame")
  
  # --------------------------------------------------
  # 5. Required columns
  # --------------------------------------------------
  
  required_cols <- c(
    "SampleID",
    "SampleRep",
    "OTU",
    "y",
    "Site",
    "site_field",
    "treat_field",
    "rep_field"
  )
  
  expect_true(
    all(required_cols %in% names(long_df))
  )
  
  # --------------------------------------------------
  # 6. Column types
  # --------------------------------------------------
  
  expect_s3_class(long_df$OTU, "factor")
  
  expect_true(is.numeric(long_df$y))
  
  expect_true(
    is.character(long_df$SampleID) ||
      is.factor(long_df$SampleID)
  )
  
  expect_true(
    is.character(long_df$SampleRep) ||
      is.factor(long_df$SampleRep)
  )
  
  # --------------------------------------------------
  # 7. SampleID: original sample identifiers
  # --------------------------------------------------
  
  expected_sample_ids <- rownames(meta_data)
  
  expect_setequal(
    unique(as.character(long_df$SampleID)),
    expected_sample_ids
  )
  
  expect_equal(
    length(unique(long_df$SampleID)),
    nrow(meta_data)
  )
  
  # Each sample must contain one row per OTU
  sample_counts <- table(long_df$SampleID)
  
  expect_true(
    all(sample_counts == length(otus))
  )
  
  # --------------------------------------------------
  # 8. SampleRep: nested sampling identifier
  # --------------------------------------------------
  
  # Use sep = "_" to match prepare_long_data()
  
  expected_nested <- interaction(
    meta_data[
      c("site_field", "treat_field", "rep_field")
    ],
    drop = TRUE,
    sep = "_"
  )
  
  expected_map <- data.frame(
    SampleID = rownames(meta_data),
    expected_SampleRep = as.character(expected_nested),
    stringsAsFactors = FALSE
  )
  
  observed_map <- unique(
    long_df[, c("SampleID", "SampleRep")]
  )
  
  observed_map$SampleID <- as.character(
    observed_map$SampleID
  )
  
  observed_map$SampleRep <- as.character(
    observed_map$SampleRep
  )
  
  # Match identifiers independently of row order
  idx <- match(
    expected_map$SampleID,
    observed_map$SampleID
  )
  
  expect_false(anyNA(idx))
  
  expect_equal(
    observed_map$SampleRep[idx],
    expected_map$expected_SampleRep
  )
  
  # Check that the number of nested groups is correct
  expect_equal(
    length(unique(long_df$SampleRep)),
    length(unique(expected_nested))
  )
  
  # --------------------------------------------------
  # 9. Site column correctness
  # --------------------------------------------------
  
  expect_false(anyNA(long_df$Site))
  
  expect_equal(
    as.character(long_df$Site),
    as.character(long_df$site_field)
  )
  
  expect_setequal(
    unique(as.character(long_df$Site)),
    sites
  )
  
  # --------------------------------------------------
  # 10. Dimensions
  # --------------------------------------------------
  
  expect_equal(
    nrow(long_df),
    length(otus) * nrow(meta_data)
  )
  
  expect_equal(
    length(unique(long_df$OTU)),
    length(otus)
  )
  
  # --------------------------------------------------
  # 11. Metadata preservation
  # --------------------------------------------------
  
  expect_true(
    all(
      c("treat_field", "rep_field", "unique_id") %in%
        names(long_df)
    )
  )
  
  # --------------------------------------------------
  # 12. Check exact count preservation
  # --------------------------------------------------
  
  expected_counts <- as.data.frame(
    as.table(otu_mat),
    stringsAsFactors = FALSE
  )
  
  names(expected_counts) <- c(
    "OTU",
    "SampleID",
    "expected_y"
  )
  
  expected_counts$OTU <- as.character(
    expected_counts$OTU
  )
  
  expected_counts$SampleID <- as.character(
    expected_counts$SampleID
  )
  
  observed_counts <- data.frame(
    OTU = as.character(long_df$OTU),
    SampleID = as.character(long_df$SampleID),
    y = as.numeric(long_df$y)
  )
  
  merged_counts <- merge(
    observed_counts,
    expected_counts,
    by = c("OTU", "SampleID"),
    all = TRUE
  )
  
  expect_equal(
    nrow(merged_counts),
    length(otus) * nrow(meta_data)
  )
  
  expect_false(anyNA(merged_counts$y))
  
  expect_false(anyNA(merged_counts$expected_y))
  
  expect_equal(
    merged_counts$y,
    merged_counts$expected_y
  )
  
  # --------------------------------------------------
  # 13. Sanity checks
  # --------------------------------------------------
  
  expect_true(all(long_df$y >= 0))
  
  expect_true(all(is.finite(long_df$y)))
  
  expect_false(anyNA(long_df$OTU))
  
  expect_false(anyNA(long_df$SampleID))
  
  expect_false(anyNA(long_df$SampleRep))
  
})


# ====================================================
# Additional test: repeated nested sampling groups
# ====================================================

test_that("prepare_long_data preserves repeated nested groups", {
  
  set.seed(123)
  
  # --------------------------------------------------
  # 1. Metadata with repeated nested groups
  # --------------------------------------------------
  
  meta <- data.frame(
    site_field = c(
      "Site1", "Site1", "Site1", "Site2"
    ),
    treat_field = c(
      "Control", "Control", "Rats", "Control"
    ),
    rep_field = c(1, 1, 1, 1),
    row.names = paste0("Sample", 1:4)
  )
  
  # --------------------------------------------------
  # 2. OTU table
  # --------------------------------------------------
  
  otu_mat <- matrix(
    rpois(3 * 4, lambda = 5),
    nrow = 3,
    ncol = 4,
    dimnames = list(
      paste0("OTU", 1:3),
      rownames(meta)
    )
  )
  
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(
      otu_mat,
      taxa_are_rows = TRUE
    ),
    phyloseq::sample_data(meta)
  )
  
  # --------------------------------------------------
  # 3. Run function
  # --------------------------------------------------
  
  result <- prepare_long_data(
    physeq_obj = ps,
    site_col = "site_field",
    nested_cols = c(
      "site_field",
      "treat_field",
      "rep_field"
    )
  )
  
  df <- result$long_df
  
  # --------------------------------------------------
  # 4. Check original sample IDs
  # --------------------------------------------------
  
  expect_equal(
    length(unique(df$SampleID)),
    4L
  )
  
  # --------------------------------------------------
  # 5. Check nested groups
  # --------------------------------------------------
  
  # Four samples but only three nested groups
  expect_equal(
    length(unique(df$SampleRep)),
    3L
  )
  
  # Sample1 and Sample2 must share SampleRep
  group1 <- unique(
    as.character(
      df$SampleRep[df$SampleID == "Sample1"]
    )
  )
  
  group2 <- unique(
    as.character(
      df$SampleRep[df$SampleID == "Sample2"]
    )
  )
  
  expect_length(group1, 1L)
  expect_length(group2, 1L)
  
  expect_equal(group1, group2)
  
  # Check the exact nested identifier
  expect_equal(
    group1,
    "Site1_Control_1"
  )
  
  # --------------------------------------------------
  # 6. Check dimensions and counts
  # --------------------------------------------------
  
  # --------------------------------------------------
  # 6. Check dimensions and counts
  # --------------------------------------------------
  
  expect_equal(nrow(df), 12L)
  
  expect_equal(
    length(unique(df$OTU)),
    3L
  )
  
  expect_true(all(df$y >= 0))
  expect_true(all(is.finite(df$y)))
  
})
  
