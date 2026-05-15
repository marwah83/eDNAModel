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
  # Run 3-level model
  # ------------------------------------------------------------
  out <- suppressWarnings(
    FitModel_gllvm(
      phyloseq = physeq,
      site_col = "Site",
      sample_col = "Name",
      replicate_col = "Replicate",

      abundance_rhs = (1 | OTU) + (1 | Name / OTU),

      capture_formula = a_sim ~ 1 + (1 | OTU),

      occupancy_covars = NULL,
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
    "summary",
    "psi_list", "capture_list", "lambda_list", "p_detect_list",
    "occupancy_models", "capture_models", "abundance_models",
    "site_data", "sample_data", "long_df",
    "lv_sites", "lv_species", "mean_lv_sites", "mean_lv_species",
    "filter_summary"
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
  expected_length <- 3 - 1  # n_iter - burn_in

  expect_length(out$psi_list, expected_length)
  expect_length(out$capture_list, expected_length)
  expect_length(out$lambda_list, expected_length)
  expect_length(out$p_detect_list, expected_length)

  # ------------------------------------------------------------
  # Summary columns exist
  # ------------------------------------------------------------
  expect_true(all(c(
    "psi_mean",
    "capture_mean",
    "lambda_mean",
    "p_detect_mean"
  ) %in% names(out$summary)))

  # ------------------------------------------------------------
  # Probability sanity checks
  # ------------------------------------------------------------
  expect_true(all(out$summary$psi_mean >= 0 & out$summary$psi_mean <= 1, na.rm = TRUE))
  expect_true(all(out$summary$capture_mean >= 0 & out$summary$capture_mean <= 1, na.rm = TRUE))
  expect_true(all(out$summary$p_detect_mean >= 0 & out$summary$p_detect_mean <= 1, na.rm = TRUE))

  expect_true(all(out$summary$lambda_mean >= 0, na.rm = TRUE))

})
