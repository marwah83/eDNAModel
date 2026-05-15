test_that("FitModel_gllvm runs and returns expected structure with synthetic data (3 Sites)", {

  skip_if_not_installed("gllvm")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("phyloseq")

  # ----------------------------
  # Data
  # ----------------------------
  otu_mat <- matrix(
    c(5,2,3,4,6,1,
      1,4,4,3,3,2,
      3,1,2,2,5,4),
    nrow = 3, byrow = TRUE
  )

  rownames(otu_mat) <- paste0("OTU", 1:3)
  colnames(otu_mat) <- paste0("S", 1:6)

  otu_tab <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

  sample_df <- data.frame(
    Site      = rep(c("Loc1","Loc2","Loc3"), each = 2),
    Name      = paste0("Sample", 1:6),
    Replicate = rep(1:2, 3),
    row.names = paste0("S", 1:6)
  )

  physeq <- phyloseq::phyloseq(
    otu_tab,
    phyloseq::sample_data(sample_df)
  )

  # ----------------------------
  # Run model (FIXED: add sample_col + capture_formula)
  # ----------------------------
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
      num_lv_c = 1,
      verbose = FALSE
    )
  )

  # ----------------------------
  # Structure (UPDATED for 3-level model)
  # ----------------------------
  expected_components <- c(
    "summary",
    "capture",
    "capture_site",

    "psi_list",
    "capture_list",
    "lambda_list",
    "p_detect_list",

    "occupancy_models",
    "capture_models",
    "abundance_models",

    "reduced_data",
    "sample_data",
    "long_df",

    "lv_sites",
    "lv_species",
    "mean_lv_sites",
    "mean_lv_species",

    "filter_summary",
    "note"
  )

  expect_true(all(expected_components %in% names(out)))

  # ----------------------------
  # Basic checks
  # ----------------------------
  expect_gt(nrow(out$summary), 0)

  expect_s3_class(out$lv_sites, "data.frame")
  expect_s3_class(out$mean_lv_species, "data.frame")

  # ----------------------------
  # Iteration length
  # ----------------------------
  expected_length <- 3 - 1  # n_iter - burn_in

  expect_length(out$psi_list, expected_length)
  expect_length(out$capture_list, expected_length)
  expect_length(out$lambda_list, expected_length)
  expect_length(out$p_detect_list, expected_length)

  # ----------------------------
  # Summary columns (UPDATED)
  # ----------------------------
  expect_true(all(c(
    "psi_mean",
    "lambda_mean",
    "p_detect_mean"
  ) %in% names(out$summary)))

  # capture is NOT inside summary anymore
  expect_true(all(c(
    "capture_mean",
    "capture_median",
    "capture_lwr",
    "capture_upr"
  ) %in% names(out$capture)))

  # ----------------------------
  # Sanity checks
  # ----------------------------
  expect_true(all(out$summary$psi_mean >= 0 & out$summary$psi_mean <= 1, na.rm = TRUE))
  expect_true(all(out$summary$p_detect_mean >= 0 & out$summary$p_detect_mean <= 1, na.rm = TRUE))
  expect_true(all(out$summary$lambda_mean >= 0, na.rm = TRUE))

})
