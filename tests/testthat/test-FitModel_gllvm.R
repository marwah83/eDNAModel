test_that("FitModel_gllvm runs and returns expected structure with synthetic data", {
  skip_if_not_installed("gllvm")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("phyloseq")

  # Create synthetic OTU matrix (samples x OTUs)
  otu_mat <- matrix(
    sample(0:10, 30, replace = TRUE),
    nrow = 10, ncol = 3,
    dimnames = list(
      paste0("Sample", 1:10),
      paste0("OTU", 1:3)
    )
  )

  # Create synthetic sample metadata
  meta_df <- data.frame(
    Name = rownames(otu_mat),
    Samplingmonth = sample(c("Jan", "Feb", "Mar"), 10, replace = TRUE),
    stringsAsFactors = TRUE
  )

  # Build phyloseq object
  OTU <- phyloseq::otu_table(otu_mat, taxa_are_rows = FALSE)
  SAM <- phyloseq::sample_data(meta_df)
  physeq_obj <- phyloseq::phyloseq(OTU, SAM)

  # Run the model with suppressed warnings for test cleanliness
  out <- suppressWarnings({
    FitModel_gllvm(
      phyloseq = physeq_obj,
      site_col = "Name",
      abundance_rhs = (1 | Samplingmonth / OTU) + (1 | OTU),
      occupancy_covars = "Samplingmonth",
      abundance_family = "poisson",
      n_iter = 3,
      burn_in = 1
    )
  })

  # === Assertions ===
  expect_type(out, "list")

  expected_components <- c(
    "summary", "psi_list", "lambda_list", "p_detect_list",
    "occupancy_models", "abundance_models", "reduced_data",
    "lv_sites", "lv_species", "mean_lv_sites", "mean_lv_species"
  )
  expect_true(all(expected_components %in% names(out)))

  expect_gt(nrow(out$summary), 0)
  expect_s3_class(out$lv_sites, "data.frame")
  expect_s3_class(out$mean_lv_species, "data.frame")
})
