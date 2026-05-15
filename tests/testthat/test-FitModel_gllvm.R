test_that("FitModel_gllvm runs and returns expected structure with synthetic data (3 Sites)", {

  skip_if_not_installed("gllvm")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("phyloseq")

  # ------------------------------------------------------------
  # Synthetic OTU table: 3 OTUs × 6 Samples
  # ------------------------------------------------------------
  otu_mat <- matrix(
    c(5, 2, 3, 4, 6, 1,
      1, 4, 4, 3, 3, 2,
      3, 1, 2, 2, 5, 4),
    nrow = 3, byrow = TRUE
  )
  rownames(otu_mat) <- paste0("OTU", 1:3)
  colnames(otu_mat) <- paste0("S", 1:6)

  otu_tab <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

  # ------------------------------------------------------------
  # Sample metadata: 3 Sites × 2 replicates
  # ------------------------------------------------------------
  sample_df <- data.frame(
    Site          = rep(c("Loc1", "Loc2", "Loc3"), each = 2),
    Name          = paste0("Sample", 1:6),
    Samplingmonth = factor(c("Jan", "Feb", "Jan", "Feb", "Mar", "Mar")),
    Replicate     = rep(1:2, 3),
    row.names     = paste0("S", 1:6)
  )

  sample_tab <- phyloseq::sample_data(sample_df)

  physeq <- phyloseq::phyloseq(otu_tab, sample_tab)

  # ------------------------------------------------------------
  # Run model (simplified + stable specification)
  # ------------------------------------------------------------
  out <- suppressWarnings(
    FitModel_gllvm(
      phyloseq = physeq,
      site_col = "Site",
      abundance_rhs = (1 | OTU) + (1 | Name / OTU),
      occupancy_covars = NULL,   # <-- FIXED (must be site-level)
      abundance_family = "poisson",
      min_species_sum = 1,
      min_detection_replicates = 1,
      n_iter = 3,
      burn_in = 1,
      num_lv_c = 1
    )
  )

  # ------------------------------------------------------------
  # Structure tests
  # ------------------------------------------------------------
  expect_type(out, "list")

  expected_components <- c(
    "summary", "psi_list", "lambda_list", "p_detect_list",
    "occupancy_models", "abundance_models", "reduced_data",
    "lv_sites", "lv_species", "mean_lv_sites", "mean_lv_species",
    "filter_summary"   # <-- FIXED
  )

  expect_true(all(expected_components %in% names(out)))

  # ------------------------------------------------------------
  # Basic output checks
  # ------------------------------------------------------------
  expect_gt(nrow(out$summary), 0)

  expect_s3_class(out$lv_sites, "data.frame")
  expect_s3_class(out$mean_lv_species, "data.frame")

  # ------------------------------------------------------------
  # Iteration-related checks
  # ------------------------------------------------------------
  expect_length(out$psi_list, 2)       # n_iter - burn_in = 2
  expect_length(out$lambda_list, 2)
  expect_length(out$p_detect_list, 2)

  # ------------------------------------------------------------
  # Summary columns exist
  # ------------------------------------------------------------
  expect_true(all(c("psi_mean", "lambda_mean", "p_detect_mean") %in% names(out$summary)))

})
