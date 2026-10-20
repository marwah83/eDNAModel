
library(testthat)

test_that("FitModel_joint prepares valid TMB structures from real data", {

  skip_if_not_installed("phyloseq")
  skip_if_not_installed("TMB")
  skip_if_not_installed("dplyr")

  # ============================================================
  # 1. Load real data
  # ============================================================

  data_file <- "~/Desktop/Diversea/data.Rdata"

  if (!file.exists(data_file)) {
    skip("Real data.Rdata file is unavailable")
  }

  env <- new.env(parent = globalenv())
  load(data_file, envir = env)

  expect_true(
    exists("y", envir = env, inherits = FALSE)
  )

  y <- get("y", envir = env, inherits = FALSE)

  expect_equal(length(dim(y)), 3L)

  # Original dimensions:
  # Species x Sites x Replicates

  n_otu <- dim(y)[1]
  n_sites <- dim(y)[2]
  n_reps <- dim(y)[3]

  expect_gt(n_otu, 0)
  expect_gt(n_sites, 0)
  expect_gt(n_reps, 0)

  expect_true(all(is.finite(y)))
  expect_true(all(y >= 0))
  expect_true(all(y == floor(y)))

  # ============================================================
  # 2. Convert 3D array into phyloseq
  # ============================================================

  otu_names <- dimnames(y)[[1]]

  if (is.null(otu_names)) {
    otu_names <- paste0("OTU", seq_len(n_otu))
  }

  site_names <- dimnames(y)[[2]]

  if (is.null(site_names)) {
    site_names <- paste0("Site", seq_len(n_sites))
  }

  library_names <- paste0(
    "Library",
    seq_len(n_sites * n_reps)
  )

  counts <- matrix(
    0,
    nrow = n_otu,
    ncol = n_sites * n_reps,
    dimnames = list(otu_names, library_names)
  )

  metadata <- data.frame(
    Site = rep(site_names, each = n_reps),

    Name = rep(
      paste0("Sample", seq_len(n_sites)),
      each = n_reps
    ),

    Replicate = rep(
      paste0("R", seq_len(n_reps)),
      times = n_sites
    ),

    row.names = library_names,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_sites)) {
    for (r in seq_len(n_reps)) {

      column <- (i - 1L) * n_reps + r

      counts[, column] <- y[, i, r]
    }
  }

  expect_equal(
    dim(counts),
    c(n_otu, n_sites * n_reps)
  )

  expect_equal(
    sum(counts),
    sum(y)
  )

  expect_true(all(is.finite(counts)))
  expect_true(all(counts >= 0))

  ps_test <- phyloseq::phyloseq(
    phyloseq::otu_table(
      counts,
      taxa_are_rows = TRUE
    ),
    phyloseq::sample_data(metadata)
  )

  expect_s4_class(ps_test, "phyloseq")

  expect_equal(
    phyloseq::ntaxa(ps_test),
    n_otu
  )

  expect_equal(
    phyloseq::nsamples(ps_test),
    n_sites * n_reps
  )

  # ============================================================
  # 3. Fit joint model and capture numerical warnings
  # ============================================================

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

      # Simplified model for structural testing
      random_occ_otu = FALSE,
      random_capture_otu = FALSE,
      random_abund_otu = FALSE,

      random_sample = FALSE,
      random_sample_otu = FALSE,
      random_zi_otu = FALSE,

      gradient_tol = 1e-3,
      gradient_marginal_factor = 5,

      max_restarts = 2L,
      iter_max = 5000L,
      eval_max = 10000L,
      rel_tol = 1e-10,

      objective_rel_tol = 1e-8,
      gradient_improvement_tol = 1e-4,

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
  # 4. Check returned object
  # ============================================================

  expect_s3_class(fit, "eDNAModel_joint")

  expect_true(is.list(fit$data_tmb))
  expect_true(is.list(fit$convergence))

  expect_true(is.finite(fit$fit$objective))

  # ============================================================
  # 5. Check required TMB data fields
  # ============================================================

  dat <- fit$data_tmb

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
    "family_code",
    "use_occ_otu",
    "use_cap_otu",
    "use_abund_otu",
    "use_sample_re",
    "use_sample_otu_re",
    "use_zi_otu"
  )

  expect_true(
    all(expected_keys %in% names(dat)),
    info = paste(
      "Missing fields:",
      paste(
        setdiff(expected_keys, names(dat)),
        collapse = ", "
      )
    )
  )

  # ============================================================
  # 6. Check retained OTUs
  # ============================================================

  n_retained_otu <- length(fit$retained_otus)

  expect_gt(n_retained_otu, 0)

  expect_lte(n_retained_otu, n_otu)

  # ============================================================
  # 7. Check TMB hierarchy dimensions
  # ============================================================

  expect_length(
    dat$y,
    n_sites * n_reps * n_retained_otu
  )

  expect_equal(
    nrow(dat$X_occ),
    n_sites * n_retained_otu
  )

  expect_equal(
    nrow(dat$X_cap),
    n_sites * n_retained_otu
  )

  expect_equal(
    nrow(dat$X_abund),
    length(dat$y)
  )

  expect_equal(
    dat$n_site_groups,
    n_sites * n_retained_otu
  )

  expect_equal(
    dat$n_sample_groups,
    n_sites * n_retained_otu
  )

  # ============================================================
  # 8. Check fixed-effect design matrices
  # ============================================================

  expect_true(is.matrix(dat$X_occ))
  expect_true(is.matrix(dat$X_cap))
  expect_true(is.matrix(dat$X_abund))

  expect_equal(ncol(dat$X_occ), 1L)
  expect_equal(ncol(dat$X_cap), 1L)
  expect_equal(ncol(dat$X_abund), 1L)

  expect_true(all(dat$X_occ == 1))
  expect_true(all(dat$X_cap == 1))
  expect_true(all(dat$X_abund == 1))

  # ============================================================
  # 9. Check grouping indices
  # ============================================================

  expect_length(
    dat$row_sample_group,
    length(dat$y)
  )

  expect_length(
    dat$sample_site_group,
    dat$n_sample_groups
  )

  expect_length(
    dat$row_otu,
    length(dat$y)
  )

  expect_length(
    dat$row_sample_id,
    length(dat$y)
  )

  expect_length(
    dat$sample_otu,
    dat$n_sample_groups
  )

  expect_length(
    dat$occ_otu,
    dat$n_site_groups
  )

  expect_true(all(
    dat$row_sample_group >= 0 &
      dat$row_sample_group < dat$n_sample_groups
  ))

  expect_true(all(
    dat$sample_site_group >= 0 &
      dat$sample_site_group < dat$n_site_groups
  ))

  expect_true(all(
    dat$row_otu >= 0 &
      dat$row_otu < n_retained_otu
  ))

  expect_true(all(
    dat$occ_otu >= 0 &
      dat$occ_otu < n_retained_otu
  ))

  expect_true(all(
    dat$sample_otu >= 0 &
      dat$sample_otu < n_retained_otu
  ))

  # ============================================================
  # 10. Check observed detection indicators
  # ============================================================

  expect_length(
    dat$site_positive,
    dat$n_site_groups
  )

  expect_length(
    dat$sample_positive,
    dat$n_sample_groups
  )

  expect_true(all(
    dat$site_positive %in% c(0L, 1L)
  ))

  expect_true(all(
    dat$sample_positive %in% c(0L, 1L)
  ))

  # Observed sample detections must imply site detections.
  expect_true(all(
    dat$sample_positive <=
      dat$site_positive[dat$sample_site_group + 1L]
  ))

  # ============================================================
  # 11. Check abundance and model settings
  # ============================================================

  expect_equal(dat$family_code, 0L)

  expect_true(all(dat$offset_abund == 0))

  expect_equal(dat$use_occ_otu, 0L)
  expect_equal(dat$use_cap_otu, 0L)
  expect_equal(dat$use_abund_otu, 0L)

  expect_equal(dat$use_sample_re, 0L)
  expect_equal(dat$use_sample_otu_re, 0L)
  expect_equal(dat$use_zi_otu, 0L)

  # ============================================================
  # 12. Check TMB objective
  # ============================================================

  expect_true(is.list(fit$tmb_object))

  expect_true(is.function(fit$tmb_object$fn))
  expect_true(is.function(fit$tmb_object$gr))

  expect_true(is.finite(fit$fit$objective))

  # ============================================================
  # 13. Check numerical diagnostic structure
  # ============================================================

  conv <- fit$convergence

  expected_convergence_fields <- c(
    "optimizer_code",
    "optimizer_ok",
    "max_abs_gradient",
    "gradient_status",
    "gradient_finite",
    "pd_hessian",
    "finite_standard_errors",
    "overall_status",
    "optimization_history",
    "n_optimization_passes",
    "selected_pass"
  )

  expect_true(all(
    expected_convergence_fields %in% names(conv)
  ))

  expect_true(
    conv$overall_status %in%
      c("PASS", "MARGINAL", "FAIL")
  )

  expect_true(
    conv$gradient_status %in%
      c("PASS", "MARGINAL", "FAIL")
  )

  expect_true(is.logical(conv$optimizer_ok))
  expect_true(is.logical(conv$pd_hessian))
  expect_true(is.logical(conv$finite_standard_errors))

  expect_true(is.list(conv$optimization_history))

  expect_gte(conv$n_optimization_passes, 1L)
  expect_lte(conv$n_optimization_passes, 3L)

  expect_gte(conv$selected_pass, 1L)
  expect_lte(
    conv$selected_pass,
    conv$n_optimization_passes
  )

  # ============================================================
  # 14. Check memory settings
  # ============================================================

  expect_false(
    fit$sdreport_settings$getReportCovariance
  )

  expect_false(
    fit$sdreport_settings$getJointPrecision
  )

  # ============================================================
  # 15. Report diagnostics and captured warnings
  # ============================================================

  diagnostics <- data.frame(
    optimizer_code = conv$optimizer_code,
    max_abs_gradient = conv$max_abs_gradient,
    gradient_status = conv$gradient_status,
    pd_hessian = conv$pd_hessian,
    finite_standard_errors = conv$finite_standard_errors,
    overall_status = conv$overall_status
  )

  print(diagnostics)

  if (length(fit_warnings) > 0L) {

    message(
      "\nNumerical warnings captured: ",
      length(fit_warnings)
    )

    for (w in unique(fit_warnings)) {
      message(" - ", w)
    }
  }

})
