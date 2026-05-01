test_that("FitModel runs with minimal arguments and inferred metadata", {

  library(phyloseq)
  library(glmmTMB)
  library(dplyr)

  # -------------------------------
  # Simulate small eDNA-like dataset
  # -------------------------------
  otu_mat <- matrix(
    c(5, 2, 3, 1, 2, 1,
      1, 2, 4, 1, 3, 1,
      3, 1, 1, 2, 1, 1),
    nrow = 3, byrow = TRUE
  )

  rownames(otu_mat) <- paste0("OTU", 1:3)
  colnames(otu_mat) <- paste0("S", 1:6)

  otu_tab <- otu_table(otu_mat, taxa_are_rows = TRUE)

  sample_df <- data.frame(
    Site      = rep(c("Loc1", "Loc2"), each = 3),
    Name      = paste0("Sample", 1:6),
    treatment = c("control", "control", "control", "rats", "rats", "rats"),
    Replicate = rep(1:3, 2),
    row.names = paste0("S", 1:6)
  )

  physeq <- phyloseq(otu_tab, sample_data(sample_df))

  # -------------------------------
  # Run FitModel
  # -------------------------------
  n_iter_test <- 5
  burn_in_test <- 2

  result <- suppressWarnings(
    FitModel(
      phyloseq = physeq,
      site_col = "Site",
      sample_col = "Name",
      replicate_col = "Replicate",
      otu_col = "OTU",
      count_col = "y",

      occupancy_formula = z_sim ~ 1 + (1 | OTU),
      capture_formula   = a_sim ~ 1 + (1 | OTU),
      abundance_formula = y ~ (1 | OTU),

      abundance_family = "poisson",
      min_species_sum = 1,
      abundance_threshold = 1,
      n_iter = n_iter_test,
      burn_in = burn_in_test
    )
  )

  # -------------------------------
  # Output structure
  # -------------------------------
  expect_type(result, "list")

  expect_named(result, c(
    "psi", "capture", "lambda", "p_detect",
    "psi_list", "capture_list", "lambda_list", "p_detect_list",
    "occupancy_models", "capture_models", "abundance_models",
    "site_data", "sample_data", "long_df",
    "filter_summary", "diagnostic_AIC", "note"
  ), ignore.order = TRUE)

  # -------------------------------
  # Check main outputs exist
  # -------------------------------
  expect_s3_class(result$psi, "data.frame")
  expect_s3_class(result$capture, "data.frame")
  expect_s3_class(result$lambda, "data.frame")
  expect_s3_class(result$p_detect, "data.frame")

  expect_gt(nrow(result$psi), 0)
  expect_gt(nrow(result$lambda), 0)

  # -------------------------------
  # Column names (CRITICAL)
  # -------------------------------
  expect_true("psi_mean" %in% colnames(result$psi))
  expect_true("capture_mean" %in% colnames(result$capture))
  expect_true("lambda_mean" %in% colnames(result$lambda))
  expect_true("p_detect_mean" %in% colnames(result$p_detect))

  # -------------------------------
  # Models stored correctly
  # -------------------------------
  expect_length(result$occupancy_models, n_iter_test)
  expect_length(result$capture_models, n_iter_test)
  expect_length(result$abundance_models, n_iter_test)

  expect_true(all(vapply(result$occupancy_models, inherits, logical(1), "glmmTMB")))
  expect_true(all(vapply(result$capture_models, inherits, logical(1), "glmmTMB")))
  expect_true(all(vapply(result$abundance_models, inherits, logical(1), "glmmTMB")))

  # -------------------------------
  # Posterior draws after burn-in
  # -------------------------------
  expected_length <- n_iter_test - burn_in_test

  expect_length(result$psi_list, expected_length)
  expect_length(result$capture_list, expected_length)
  expect_length(result$lambda_list, expected_length)
  expect_length(result$p_detect_list, expected_length)

  # -------------------------------
  # Sanity checks (ROBUST)
  # -------------------------------
  expect_true(is.numeric(result$psi$psi_mean))
  expect_true(is.numeric(result$capture$capture_mean))
  expect_true(is.numeric(result$lambda$lambda_mean))
  expect_true(is.numeric(result$p_detect$p_detect_mean))

  # Allow NA (tiny data), but require at least some valid values
  expect_true(any(is.finite(result$psi$psi_mean)))
  expect_true(any(is.finite(result$lambda$lambda_mean)))

  # -------------------------------
  # Logical constraints (NEW — VERY IMPORTANT)
  # -------------------------------
  # Probabilities must be in [0,1]
  expect_true(all(result$psi$psi_mean >= 0 & result$psi$psi_mean <= 1, na.rm = TRUE))
  expect_true(all(result$capture$capture_mean >= 0 & result$capture$capture_mean <= 1, na.rm = TRUE))
  expect_true(all(result$p_detect$p_detect_mean >= 0 & result$p_detect$p_detect_mean <= 1, na.rm = TRUE))

  # Lambda must be non-negative
  expect_true(all(result$lambda$lambda_mean >= 0, na.rm = TRUE))

})
