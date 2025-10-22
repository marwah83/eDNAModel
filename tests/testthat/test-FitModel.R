test_that("FitModel", {
  library(phyloseq)
  library(glmmTMB)
  library(dplyr)
  library(tidyr)

  # -------------------------------
  # Simulate small eDNA-like dataset
  # -------------------------------

  # OTU table: 3 OTUs × 6 samples
  otu_mat <- matrix(
    c(5, 2, 3, 1, 2, 1,
      1, 2, 4, 1, 3, 1,
      3, 1, 1, 2, 1, 1),
    nrow = 3, byrow = TRUE
  )
  rownames(otu_mat) <- paste0("OTU", 1:3)
  colnames(otu_mat) <- paste0("S", 1:6)
  otu_tab <- otu_table(otu_mat, taxa_are_rows = TRUE)

  # Sample metadata
  sample_df <- data.frame(
    sampletype = rep("biologicalsample", 6),
    location   = rep(c("Loc1", "Loc2"), each = 3),
    treatment  = c("control", "control", "control", "rats", "rats", "rats"),
    Replicate  = rep(1:3, 2),
    row.names  = paste0("S", 1:6)
  )
  sample_tab <- sample_data(sample_df)

  physeq <- phyloseq(otu_tab, sample_tab)

  # -------------------------------
  # Run FitModel
  # -------------------------------
  n_iter_test <- 5
  burn_in_test <- 2

  result <- suppressWarnings(
    FitModel(
      phyloseq = physeq,
      poisson_rhs = quote((1 | OTU) + (1 | Site/OTU) + (1 | Sample/OTU) + (1 | Replicate/OTU) + treatment * OTU),
      binomial_rhs = quote((1 | OTU) + (1 | Site/OTU)),
      sampletype_keep = "biologicalsample",
      min_species_sum = 1,
      abundance_threshold = 1,
      treatment_exclude = NULL,
      n_iter = n_iter_test,
      burn_in = burn_in_test
    )
  )

  # -------------------------------
  # Structure checks
  # -------------------------------
  expect_type(result, "list")

  expect_named(result, c(
    "summary", "psi_list", "lambda_list", "p_detect_list",
    "binomial_models", "poisson_models", "reduced_data"
  ))

  expect_s3_class(result$summary, "data.frame")
  expect_gt(nrow(result$summary), 0)

  expect_true(all(c("psi_mean", "lambda_mean", "p_detect_mean") %in% names(result$summary)))

  # -------------------------------
  # Model storage checks
  # -------------------------------
  expect_length(result$binomial_models, n_iter_test)
  expect_length(result$poisson_models, n_iter_test)

  expect_true(all(sapply(result$binomial_models, inherits, "glmmTMB")))
  expect_true(all(sapply(result$poisson_models, inherits, "glmmTMB")))

  # -------------------------------
  # Posterior list lengths (after burn-in)
  # -------------------------------
  expected_length <- n_iter_test - burn_in_test

  expect_length(result$psi_list, expected_length)
  expect_length(result$lambda_list, expected_length)
  expect_length(result$p_detect_list, expected_length)
})
