
library(testthat)

test_that("FitModel_joint correctly prepares hierarchical TMB data", {

  skip_if_not_installed("TMB")
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("dplyr")

  # ============================================================
  # 1. Synthetic data: species, sites, biological samples,
  #    and sequencing replicates
  # ============================================================

  set.seed(123)

  n_otu <- 3
  n_sites <- 2
  samples_per_site <- 2
  replicates_per_sample <- 4

  n_bio_samples <- n_sites * samples_per_site
  n_libraries <- n_bio_samples * replicates_per_sample

  otu_names <- paste0("OTU", seq_len(n_otu))
  library_names <- paste0("Library", seq_len(n_libraries))

  sample_df <- data.frame(
    Site = rep(
      paste0("Site", seq_len(n_sites)),
      each = samples_per_site * replicates_per_sample
    ),
    Name = rep(
      paste0("Sample", seq_len(n_bio_samples)),
      each = replicates_per_sample
    ),
    Replicate = rep(
      paste0("R", seq_len(replicates_per_sample)),
      times = n_bio_samples
    ),
    row.names = library_names,
    stringsAsFactors = FALSE
  )

  counts <- matrix(
    rpois(n_otu * n_libraries, lambda = 5),
    nrow = n_otu,
    ncol = n_libraries,
    dimnames = list(otu_names, library_names)
  )

  # Ensure every OTU has detections
  counts[, 1] <- counts[, 1] + 2

  ps_test <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::sample_data(sample_df)
  )

  # ============================================================
  # 2. Fit the joint model
  # ============================================================

  fit <- FitModel_joint(
    phyloseq = ps_test,
    site_col = "Site",
    sample_col = "Name",
    replicate_col = "Replicate",

    occupancy_formula = ~ 1,
    capture_formula = ~ 1,
    abundance_formula = ~ 1,

    abundance_family = "poisson",

    min_species_sum = 1,
    min_detection_replicates = 1,

    random_occ_otu = TRUE,
    random_capture_otu = TRUE,
    random_abund_otu = TRUE,

    random_sample = FALSE,
    random_sample_otu = FALSE,
    random_zi_otu = FALSE,

    max_restarts = 0,
    iter_max = 1000,
    eval_max = 2000,

    get_report_covariance = FALSE,
    get_joint_precision = FALSE,
    otu_intervals = FALSE,

    verbose = FALSE
  )

  # ============================================================
  # 3. Check returned TMB data
  # ============================================================

  expect_s3_class(fit, "eDNAModel_joint")

  prepared_data <- fit$data_tmb

  expect_true(is.list(prepared_data))

  expected_keys <- c(
    "y",
    "X_occ",
    "X_cap",
    "X_abund",
    "offset_abund",
    "occ_otu",
    "sample_otu",
    "row_sample_group",
    "sample_site_group",
    "row_otu",
    "row_sample_id",
    "row_sample_otu_re",
    "sample_positive",
    "site_positive",
    "n_site_groups",
    "n_sample_groups",
    "family_code"
  )

  expect_true(
    all(expected_keys %in% names(prepared_data))
  )

  # ============================================================
  # 4. Check hierarchy dimensions
  # ============================================================

  expect_length(
    prepared_data$y,
    n_otu * n_libraries
  )

  expect_equal(
    nrow(prepared_data$X_occ),
    n_sites * n_otu
  )

  expect_equal(
    nrow(prepared_data$X_cap),
    n_bio_samples * n_otu
  )

  expect_equal(
    nrow(prepared_data$X_abund),
    n_libraries * n_otu
  )

  expect_equal(
    prepared_data$n_site_groups,
    n_sites * n_otu
  )

  expect_equal(
    prepared_data$n_sample_groups,
    n_bio_samples * n_otu
  )

  # ============================================================
  # 5. Check index mappings
  # ============================================================

  expect_length(
    prepared_data$row_sample_group,
    length(prepared_data$y)
  )

  expect_length(
    prepared_data$sample_site_group,
    prepared_data$n_sample_groups
  )

  expect_true(all(
    prepared_data$row_sample_group >= 0 &
    prepared_data$row_sample_group <
      prepared_data$n_sample_groups
  ))

  expect_true(all(
    prepared_data$sample_site_group >= 0 &
    prepared_data$sample_site_group <
      prepared_data$n_site_groups
  ))

  # ============================================================
  # 6. Check detection indicators
  # ============================================================

  expect_length(
    prepared_data$sample_positive,
    prepared_data$n_sample_groups
  )

  expect_length(
    prepared_data$site_positive,
    prepared_data$n_site_groups
  )

  expect_true(all(
    prepared_data$sample_positive %in% c(0L, 1L)
  ))

  expect_true(all(
    prepared_data$site_positive %in% c(0L, 1L)
  ))

  # ============================================================
  # 7. Check Poisson family and offset
  # ============================================================

  expect_equal(prepared_data$family_code, 0L)

  expect_true(all(
    prepared_data$offset_abund == 0
  ))

  # ============================================================
  # 8. Check TMB object
  # ============================================================

  expect_true(is.list(fit$tmb_object))

  expect_true(is.function(fit$tmb_object$fn))
  expect_true(is.function(fit$tmb_object$gr))

  expect_true(is.finite(fit$fit$objective))

  # ============================================================
  # 9. Check memory settings
  # ============================================================

  expect_false(
    fit$sdreport_settings$getReportCovariance
  )

  expect_false(
    fit$sdreport_settings$getJointPrecision
  )

})
