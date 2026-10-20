
library(testthat)

test_that("FitModel_joint correctly prepares TMB data", {

  skip_if_not_installed("phyloseq")
  skip_if_not_installed("TMB")
  skip_if_not_installed("dplyr")

  # ============================================================
  # 1. Create synthetic data
  # ============================================================

  set.seed(123)

  n_otu <- 3
  n_sites <- 2
  n_replicates <- 4

  n_libraries <- n_sites * n_replicates

  otu_names <- paste0("Species_", seq_len(n_otu))
  library_names <- paste0("Library_", seq_len(n_libraries))

  counts <- matrix(
    rpois(n_otu * n_libraries, lambda = 5),
    nrow = n_otu,
    ncol = n_libraries,
    dimnames = list(otu_names, library_names)
  )

  # Ensure each OTU has positive counts
  counts[, 1] <- counts[, 1] + 2

  # ============================================================
  # 2. Construct sample metadata
  # ============================================================

  sample_df <- data.frame(
    Site = rep(
      paste0("Site_", seq_len(n_sites)),
      each = n_replicates
    ),

    Name = rep(
      paste0("Sample_", seq_len(n_sites)),
      each = n_replicates
    ),

    Replicate = rep(
      paste0("R", seq_len(n_replicates)),
      times = n_sites
    ),

    row.names = library_names,
    stringsAsFactors = FALSE
  )

  # ============================================================
  # 3. Create phyloseq object
  # ============================================================

  ps_test <- phyloseq::phyloseq(
    phyloseq::otu_table(
      counts,
      taxa_are_rows = TRUE
    ),
    phyloseq::sample_data(sample_df)
  )

  expect_s4_class(ps_test, "phyloseq")

  # ============================================================
  # 4. Fit simple joint model
  # ============================================================

  # Capture numerical warnings separately because this test
  # is intended to validate TMB data structures.

  fit_warnings <- character()

  fit <- withCallingHandlers(

    FitModel_joint(
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

      random_occ_otu = FALSE,
      random_capture_otu = FALSE,
      random_abund_otu = FALSE,

      random_sample = FALSE,
      random_sample_otu = FALSE,
      random_zi_otu = FALSE,

      max_restarts = 0L,
      iter_max = 1000L,
      eval_max = 2000L,

      get_report_covariance = FALSE,
      get_joint_precision = FALSE,
      otu_intervals = FALSE,

      verbose = FALSE
    ),

    warning = function(w) {
      fit_warnings <<- c(
        fit_warnings,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )

  # ============================================================
  # 5. Check returned data structure
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
    "sample_positive",
    "site_positive",
    "n_site_groups",
    "n_sample_groups",
    "family_code"
  )

  expect_true(
    all(expected_keys %in% names(prepared_data)),
    info = paste(
      "Missing fields:",
      paste(
        setdiff(expected_keys, names(prepared_data)),
        collapse = ", "
      )
    )
  )

  # ============================================================
  # 6. Check dimensions
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
    n_sites * n_otu
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
    n_sites * n_otu
  )

  # ============================================================
  # 7. Check grouping indices
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
  # 8. Check detection indicators
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
  # 9. Check model settings
  # ============================================================

  expect_equal(prepared_data$family_code, 0L)

  expect_true(all(
    prepared_data$offset_abund == 0
  ))

  expect_equal(prepared_data$use_occ_otu, 0L)
  expect_equal(prepared_data$use_cap_otu, 0L)
  expect_equal(prepared_data$use_abund_otu, 0L)

  # ============================================================
  # 10. Check TMB objective
  # ============================================================

  expect_true(is.list(fit$tmb_object))

  expect_true(is.function(fit$tmb_object$fn))
  expect_true(is.function(fit$tmb_object$gr))

  expect_true(is.finite(fit$fit$objective))

  # ============================================================
  # 11. Check convergence diagnostics exist
  # ============================================================

  expect_true(is.list(fit$convergence))

  expect_true(
    fit$convergence$overall_status %in%
      c("PASS", "MARGINAL", "FAIL")
  )

  # ============================================================
  # 12. Check memory settings
  # ============================================================

  expect_false(
    fit$sdreport_settings$getReportCovariance
  )

  expect_false(
    fit$sdreport_settings$getJointPrecision
  )

  # ============================================================
  # 13. Report numerical warnings
  # ============================================================

  if (length(fit_warnings) > 0L) {
    message(
      "Numerical warnings captured: ",
      paste(unique(fit_warnings), collapse = " | ")
    )
  }

})
