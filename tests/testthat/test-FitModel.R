test_that("FitModel runs with deterministic dataset and stable outputs", {
  
  library(phyloseq)
  library(glmmTMB)
  
  # -------------------------------
  # Deterministic eDNA dataset
  # -------------------------------
  
  otu_mat <- matrix(
    c(
      # OTU1 (low)
      1, 1, 1, 1, 1, 1,
      
      # OTU2 (medium)
      5, 5, 5, 5, 5, 5,
      
      # OTU3 (high)
      15, 15, 15, 15, 15, 15
    ),
    nrow = 3,
    byrow = TRUE
  )
  
  rownames(otu_mat) <- paste0("OTU", 1:3)
  colnames(otu_mat) <- paste0("S", 1:6)
  
  otu_tab <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  
  sample_df <- data.frame(
    Site      = rep(c("Loc1", "Loc2"), each = 3),
    Name      = paste0("Sample", 1:6),
    Replicate = rep(1:3, 2),
    row.names = paste0("S", 1:6)
  )
  
  physeq <- phyloseq::phyloseq(
    otu_tab,
    phyloseq::sample_data(sample_df)
  )
  
  # -------------------------------
  # Run FitModel
  # -------------------------------
  
  n_iter_test <- 6
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
      burn_in = burn_in_test,
      verbose = FALSE
    )
  )
  
  # -------------------------------
  # Expected names
  # -------------------------------
  
  expected_names <- c(
    "psi", "capture", "lambda", "p_detect",
    "psi_list", "capture_list", "lambda_list", "p_detect_list",
    "occupancy_models", "capture_models", "abundance_models",
    "site_data", "sample_data", "long_df",
    "filter_summary", "diagnostic_AIC", "note"
  )
  
  expect_true(all(expected_names %in% names(result)))
  
  # -------------------------------
  # Core outputs exist
  # -------------------------------
  
  expect_s3_class(result$psi, "data.frame")
  expect_s3_class(result$capture, "data.frame")
  expect_s3_class(result$lambda, "data.frame")
  expect_s3_class(result$p_detect, "data.frame")
  
  expect_gt(nrow(result$psi), 0)
  expect_gt(nrow(result$lambda), 0)
  
  # -------------------------------
  # Column names
  # -------------------------------
  
  expect_true(all(c("psi_mean","psi_median","psi_lwr","psi_upr") %in% names(result$psi)))
  expect_true(all(c("capture_mean","capture_median","capture_lwr","capture_upr") %in% names(result$capture)))
  expect_true(all(c("lambda_mean","lambda_median","lambda_lwr","lambda_upr") %in% names(result$lambda)))
  expect_true(all(c("p_detect_mean","p_detect_median","p_detect_lwr","p_detect_upr") %in% names(result$p_detect)))
  
  # -------------------------------
  # Model storage
  # -------------------------------
  
  expected_length <- n_iter_test - burn_in_test
  
  expect_length(result$occupancy_models, expected_length)
  expect_length(result$capture_models, expected_length)
  expect_length(result$abundance_models, expected_length)
  
  expect_length(result$psi_list, expected_length)
  expect_length(result$capture_list, expected_length)
  expect_length(result$lambda_list, expected_length)
  expect_length(result$p_detect_list, expected_length)
  
  # -------------------------------
  # Numeric sanity
  # -------------------------------
  
  expect_true(any(is.finite(result$psi$psi_mean)))
  expect_true(any(is.finite(result$lambda$lambda_mean)))
  
  # -------------------------------
  # Probability constraints
  # -------------------------------
  
  expect_true(all(result$psi$psi_mean >= 0 & result$psi$psi_mean <= 1, na.rm = TRUE))
  expect_true(all(result$capture$capture_mean >= 0 & result$capture$capture_mean <= 1, na.rm = TRUE))
  expect_true(all(result$p_detect$p_detect_mean >= 0 & result$p_detect$p_detect_mean <= 1, na.rm = TRUE))
  
  expect_true(all(result$lambda$lambda_mean >= 0, na.rm = TRUE))
  
  # -------------------------------
  # STRONG RELATIONSHIP TEST
  # λ → p_detect (monotonic)
  # -------------------------------
  
  cor_val <- suppressWarnings(
    cor(result$lambda$lambda_mean, result$p_detect$p_detect_mean, use = "complete.obs")
  )
  
  expect_true(is.finite(cor_val))
  expect_gt(cor_val, 0.9)   # deterministic dataset ensures this
  
  # -------------------------------
  # Functional relationship check
  # p_detect ≈ 1 - exp(-lambda)
  # -------------------------------
  
  approx_pd <- 1 - exp(-result$lambda$lambda_mean)
  
  expect_true(
    mean(abs(result$p_detect$p_detect_mean - approx_pd), na.rm = TRUE) < 1e-2
  )
  
})
