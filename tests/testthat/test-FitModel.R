test_that("FitModel runs with minimal arguments and inferred metadata", {
  library(phyloseq)
  library(glmmTMB)
  library(dplyr)
  library(tidyr)

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
    Site       = rep(c("Loc1", "Loc2"), each = 3),
    SampleName = paste0("Sample", 1:6),
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
      site_col = "Site",
      abundance_rhs = (1 | OTU) + (1 | Site / OTU) + (1 | SampleName / OTU) + (1 | Replicate / OTU) + treatment * OTU,
      occupancy_rhs = (1 | OTU) + (1 | Site / OTU),
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
    "summary", "psi_list", "lambda_list", "p_detect_list",
    "occupancy_models", "abundance_models", "reduced_data"
  ))

  expect_s3_class(result$summary, "data.frame")
  expect_gt(nrow(result$summary), 0)

  expect_true(any(grepl("psi_mean", names(result$summary))))
  expect_true(any(grepl("lambda_mean", names(result$summary))))
  expect_true(any(grepl("p_detect_mean", names(result$summary))))

  # -------------------------------
  # Models stored correctly
  # -------------------------------
  expect_length(result$occupancy_models, n_iter_test)
  expect_length(result$abundance_models, n_iter_test)
  expect_true(all(sapply(result$occupancy_models, inherits, "glmmTMB")))
  expect_true(all(sapply(result$abundance_models, inherits, "glmmTMB")))

  # -------------------------------
  # Posterior draws after burn-in
  # -------------------------------
  expected_length <- n_iter_test - burn_in_test
  expect_length(result$psi_list, expected_length)
  expect_length(result$lambda_list, expected_length)
  expect_length(result$p_detect_list, expected_length)
})
