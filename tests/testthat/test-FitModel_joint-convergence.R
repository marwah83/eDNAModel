# Place in tests/testthat/test-fitmodel-joint.R
# Requires the package's compiled TMB DLL and internal prepare_long_data().

testthat::test_that("FitModel_joint preserves the three-level hierarchy and diagnostics", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("phyloseq")
  testthat::skip_if_not_installed("dplyr")

  set.seed(123)

  # Six sites, three biological samples per site, two PCR/sequencing
  # replicates per biological sample, and six OTUs.
  n_sites <- 6L
  bio_per_site <- 3L
  pcr_per_bio <- 2L
  n_otu <- 6L
  n_bio <- n_sites * bio_per_site
  n_libraries <- n_bio * pcr_per_bio

  biological_id <- paste0("Bio", seq_len(n_bio))
  library_id <- paste0("Lib", seq_len(n_libraries))
  otu_id <- paste0("OTU", seq_len(n_otu))

  sample_df <- data.frame(
    Site = rep(paste0("Site", seq_len(n_sites)),
               each = bio_per_site * pcr_per_bio),
    Name = rep(biological_id, each = pcr_per_bio),
    Replicate = rep(c("R1", "R2"), times = n_bio),
    row.names = library_id,
    stringsAsFactors = FALSE
  )

  # Include genuine zero observations while avoiding an all-positive
  # matrix. Keep all OTUs represented and count data integer-valued.
  counts <- matrix(
    stats::rpois(n_otu * n_libraries, lambda = 5),
    nrow = n_otu,
    ncol = n_libraries,
    dimnames = list(otu_id, library_id)
  )
  counts[matrix(stats::runif(length(counts)) < 0.25,
                nrow = n_otu)] <- 0L
  counts[, 1L] <- counts[, 1L] + 2L

  ps_test <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::sample_data(sample_df)
  )

  testthat::expect_equal(phyloseq::nsamples(ps_test), n_libraries)
  testthat::expect_equal(phyloseq::ntaxa(ps_test), n_otu)
  testthat::expect_equal(length(unique(sample_df$Site)), n_sites)
  testthat::expect_equal(length(unique(sample_df$Name)), n_bio)
  testthat::expect_true(all(table(sample_df$Name) == pcr_per_bio))

  # The implementation analytically marginalizes occupancy/capture
  # and fits one joint TMB likelihood (not an iterative gllvm sampler).
  # Keep covariance-intensive OTU intervals off for a fast unit test.
  fit <- suppressWarnings(FitModel_joint(
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
    random_zi_otu = TRUE,  # ignored for Poisson
    gradient_tol = 1e-3,
    gradient_marginal_factor = 5,
    max_restarts = 2L,
    iter_max = 2000L,
    eval_max = 4000L,
    objective_rel_tol = 1e-8,
    gradient_improvement_tol = 1e-4,
    get_report_covariance = FALSE,
    get_joint_precision = FALSE,
    otu_intervals = FALSE,
    verbose = FALSE,
    tmb_verbose = FALSE
  ))

  testthat::expect_s3_class(fit, "eDNAModel_joint")
  required <- c(
    "fit", "tmb_object", "sdreport", "report", "fixed_effects",
    "derived", "site_data", "sample_data", "long_df", "otu_stats",
    "retained_otus", "data_tmb", "parameters_start", "parameter_map",
    "random_effects", "formulas", "abundance_family", "convergence",
    "diagnostics", "sdreport_settings", "zero_inflation_structure",
    "otu_probabilities", "otu_probability_note", "note"
  )
  testthat::expect_true(all(required %in% names(fit)))

  # Verify grouping and replicate mapping, including counts of groups.
  testthat::expect_equal(nrow(fit$site_data), n_sites * n_otu)
  testthat::expect_equal(nrow(fit$sample_data), n_bio * n_otu)
  testthat::expect_equal(nrow(fit$long_df), n_libraries * n_otu)
  testthat::expect_equal(length(fit$retained_otus), n_otu)
  testthat::expect_equal(length(fit$data_tmb$row_sample_group),
                         nrow(fit$long_df))
  testthat::expect_equal(length(unique(fit$data_tmb$row_sample_group)),
                         nrow(fit$sample_data))
  testthat::expect_equal(length(fit$data_tmb$sample_site_group),
                         nrow(fit$sample_data))
  testthat::expect_true(all(fit$data_tmb$row_sample_group >= 0L))
  testthat::expect_true(all(fit$data_tmb$sample_site_group >= 0L))
  testthat::expect_equal(fit$abundance_family, "poisson")
  testthat::expect_false(fit$zero_inflation_structure$enabled)
  testthat::expect_false(fit$sdreport_settings$getReportCovariance)
  testthat::expect_false(fit$sdreport_settings$getJointPrecision)

  conv <- fit$convergence
  conv_fields <- c(
    "converged", "strict_convergence", "acceptable_convergence",
    "overall_status", "optimizer_ok", "optimizer_status",
    "optimizer_code", "optimizer_message", "objective", "gradient",
    "max_abs_gradient", "gradient_tolerance", "gradient_marginal_factor",
    "gradient_ok", "gradient_status", "gradient_finite",
    "gradient_table", "largest_gradients", "sdreport_ok",
    "sdreport_status", "pd_hessian", "hessian_status",
    "finite_standard_errors", "se_status", "optimization_history",
    "n_optimization_passes", "selected_pass", "max_restarts",
    "objective_rel_tol", "gradient_improvement_tol", "stalled"
  )
  testthat::expect_true(all(conv_fields %in% names(conv)))
  testthat::expect_identical(conv$optimizer_ok, conv$optimizer_code == 0L)

  # Check the reported gradient against an independent TMB evaluation.
  obj <- fit$tmb_object
  opt <- fit$fit
  invisible(obj$fn(opt$par))
  g_manual <- as.numeric(obj$gr(opt$par))
  testthat::expect_equal(as.numeric(conv$gradient), g_manual,
                         tolerance = 1e-6)
  testthat::expect_identical(conv$gradient_finite,
                             length(g_manual) > 0L && all(is.finite(g_manual)))
  testthat::expect_equal(conv$max_abs_gradient,
                         if (all(is.finite(g_manual))) max(abs(g_manual)) else Inf,
                         tolerance = 1e-6)

  gmax <- conv$max_abs_gradient
  tol <- conv$gradient_tolerance
  expected_gradient_status <- if (!is.finite(gmax)) {
    "FAIL"
  } else if (gmax <= tol) {
    "PASS"
  } else if (gmax <= conv$gradient_marginal_factor * tol) {
    "MARGINAL"
  } else {
    "FAIL"
  }
  testthat::expect_identical(conv$gradient_status, expected_gradient_status)
  testthat::expect_identical(conv$gradient_ok,
                             is.finite(gmax) && gmax <= tol)

  testthat::expect_identical(conv$sdreport_ok, !is.null(fit$sdreport))
  testthat::expect_identical(conv$pd_hessian,
                             isTRUE(fit$sdreport$pdHess))
  fixed_se_ok <- nrow(fit$fixed_effects) > 0L &&
    "Std. Error" %in% names(fit$fixed_effects) &&
    all(is.finite(fit$fixed_effects[["Std. Error"]]))
  testthat::expect_identical(conv$finite_standard_errors, fixed_se_ok)

  expected_strict <- conv$optimizer_ok &&
    conv$gradient_status == "PASS" && conv$sdreport_ok &&
    conv$pd_hessian && conv$finite_standard_errors
  expected_acceptable <- conv$optimizer_ok &&
    conv$gradient_status %in% c("PASS", "MARGINAL") &&
    conv$sdreport_ok && conv$pd_hessian && conv$finite_standard_errors
  testthat::expect_identical(conv$strict_convergence, expected_strict)
  testthat::expect_identical(conv$acceptable_convergence, expected_acceptable)
  testthat::expect_identical(conv$converged, expected_acceptable)
  testthat::expect_identical(conv$overall_status,
    if (expected_strict) "PASS" else if (expected_acceptable) "MARGINAL" else "FAIL")

  # The chosen optimization pass must be among the available passes,
  # and its parameters/objective must match the returned optimizer.
  history <- conv$optimization_history
  testthat::expect_true(is.list(history))
  testthat::expect_gte(length(history), 1L)
  testthat::expect_lte(length(history), 3L)
  testthat::expect_equal(conv$n_optimization_passes, length(history))
  testthat::expect_gte(conv$selected_pass, 1L)
  testthat::expect_lte(conv$selected_pass, length(history))
  testthat::expect_equal(opt$par, history[[conv$selected_pass]]$par,
                         tolerance = 1e-8)
  testthat::expect_equal(opt$objective,
                         history[[conv$selected_pass]]$objective,
                         tolerance = 1e-8)

  objectives <- vapply(history, function(x) x$objective, numeric(1))
  gradients <- vapply(history, function(x) x$max_abs_gradient, numeric(1))
  finite_obj <- is.finite(objectives)
  testthat::expect_true(any(finite_obj))
  min_obj <- min(objectives[finite_obj])
  equivalent <- vapply(seq_along(objectives), function(i) {
    is.finite(objectives[i]) &&
      abs(objectives[i] - min_obj) <=
      conv$objective_rel_tol * max(1, abs(objectives[i]), abs(min_obj))
  }, logical(1))
  candidates <- which(equivalent & is.finite(gradients))
  expected_pass <- if (length(candidates)) {
    candidates[which.min(gradients[candidates])]
  } else {
    which.min(ifelse(finite_obj, objectives, Inf))
  }
  testthat::expect_equal(conv$selected_pass, expected_pass)
})
