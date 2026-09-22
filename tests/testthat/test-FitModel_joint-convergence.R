test_that("FitModel_joint runs and returns valid convergence diagnostics", {

  skip_if_not_installed("TMB")
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("dplyr")

  # ============================================================
  # 1. Create synthetic count data
  # ============================================================

  set.seed(123)

  n_otu <- 6
  n_sites <- 6
  reps_per_site <- 3

  n_samples <- n_sites * reps_per_site

  otu_names <- paste0("OTU", seq_len(n_otu))
  sample_names <- paste0("S", seq_len(n_samples))

  # Moderately dense synthetic data so that this unit test is
  # primarily testing implementation/numerics rather than an
  # intentionally pathological sparse-data case.
  counts <- matrix(
    stats::rpois(
      n_otu * n_samples,
      lambda = 8
    ),
    nrow = n_otu,
    ncol = n_samples
  )

  # Ensure every OTU has positive observations.
  counts[, 1] <- counts[, 1] + 2

  rownames(counts) <- otu_names
  colnames(counts) <- sample_names


  # ============================================================
  # 2. Construct sample metadata
  # ============================================================

  sample_df <- data.frame(

    Site = rep(
      paste0("Site", seq_len(n_sites)),
      each = reps_per_site
    ),

    Name = sample_names,

    Replicate = rep(
      seq_len(reps_per_site),
      times = n_sites
    ),

    row.names = sample_names,

    stringsAsFactors = FALSE
  )


  # ============================================================
  # 3. Construct phyloseq object
  # ============================================================

  otu_tab <- phyloseq::otu_table(
    counts,
    taxa_are_rows = TRUE
  )

  samp_tab <- phyloseq::sample_data(
    sample_df
  )

  ps_test <- phyloseq::phyloseq(
    otu_tab,
    samp_tab
  )


  # ============================================================
  # 4. Fit a simple Poisson version
  #
  # Keep the model deliberately simple. More complex ZIP/ZINB
  # behavior should be tested separately.
  # ============================================================

  fit <- suppressWarnings(
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

      random_occ_otu = TRUE,
      random_capture_otu = TRUE,
      random_abund_otu = TRUE,

      random_sample = TRUE,
      random_sample_otu = FALSE,

      # Ignored for Poisson, but deliberately leave default-style
      # behavior switched on to check that this is handled.
      random_zi_otu = TRUE,

      gradient_tol = 1e-3,
      gradient_marginal_factor = 5,

      max_restarts = 2,

      iter_max = 2000,
      eval_max = 4000,

      objective_rel_tol = 1e-8,
      gradient_improvement_tol = 1e-4,

      get_report_covariance = FALSE,

      # Keep this FALSE in routine unit tests to reduce memory/time.
      get_joint_precision = FALSE,

      verbose = FALSE
    )
  )


  # ============================================================
  # 5. Basic returned-object structure
  # ============================================================

  expect_s3_class(
    fit,
    "eDNAModel_joint"
  )

  expected_components <- c(
    "fit",
    "tmb_object",
    "sdreport",
    "report",
    "fixed_effects",
    "derived",
    "site_data",
    "sample_data",
    "long_df",
    "otu_stats",
    "retained_otus",
    "data_tmb",
    "parameters_start",
    "parameter_map",
    "random_effects",
    "formulas",
    "abundance_family",
    "convergence",
    "diagnostics",
    "sdreport_settings",
    "zero_inflation_structure",
    "note"
  )

  expect_true(
    all(expected_components %in% names(fit))
  )


  # ============================================================
  # 6. Test convergence structure
  # ============================================================

  conv <- fit$convergence

  expected_convergence_fields <- c(
    "converged",
    "strict_convergence",
    "acceptable_convergence",
    "overall_status",

    "optimizer_ok",
    "optimizer_status",

    "optimizer_code",
    "optimizer_message",

    "objective",

    "gradient",
    "max_abs_gradient",

    "gradient_tolerance",
    "gradient_marginal_factor",

    "gradient_ok",
    "gradient_status",
    "gradient_finite",

    "gradient_table",
    "largest_gradients",

    "sdreport_ok",
    "sdreport_status",

    "pd_hessian",
    "hessian_status",

    "finite_standard_errors",
    "se_status",

    "optimization_history",
    "n_optimization_passes",
    "selected_pass",

    "max_restarts",
    "objective_rel_tol",
    "gradient_improvement_tol",

    "stalled"
  )

  expect_true(
    all(
      expected_convergence_fields %in%
        names(conv)
    )
  )


  # ============================================================
  # 7. Optimizer code should agree with optimizer_ok
  # ============================================================

  expect_identical(
    conv$optimizer_ok,
    conv$optimizer_code == 0L
  )


  # ============================================================
  # 8. Gradient must be finite
  # ============================================================

  expect_true(
    is.numeric(conv$gradient)
  )

  expect_true(
    length(conv$gradient) > 0
  )

  expect_true(
    all(is.finite(conv$gradient))
  )

  expect_true(
    conv$gradient_finite
  )


  # ============================================================
  # 9. Verify max absolute gradient manually
  # ============================================================

  manual_max_gradient <- max(
    abs(conv$gradient)
  )

  expect_equal(
    conv$max_abs_gradient,
    manual_max_gradient,
    tolerance = 1e-8
  )


  # ============================================================
  # 10. Independently recalculate gradient from TMB object
  # ============================================================

  obj <- fit$tmb_object
  opt <- fit$fit

  # Make sure the random-effect mode corresponds to opt$par.
  invisible(
    obj$fn(opt$par)
  )

  g_manual <- obj$gr(
    opt$par
  )

  expect_equal(
    as.numeric(g_manual),
    as.numeric(conv$gradient),
    tolerance = 1e-7
  )

  expect_equal(
    max(abs(g_manual)),
    conv$max_abs_gradient,
    tolerance = 1e-7
  )


  # ============================================================
  # 11. Test gradient classification
  # ============================================================

  gmax <- conv$max_abs_gradient
  tol <- conv$gradient_tolerance

  marginal_limit <-
    conv$gradient_marginal_factor * tol

  expected_gradient_status <- if (
    !is.finite(gmax)
  ) {

    "FAIL"

  } else if (
    gmax <= tol
  ) {

    "PASS"

  } else if (
    gmax <= marginal_limit
  ) {

    "MARGINAL"

  } else {

    "FAIL"
  }

  expect_identical(
    conv$gradient_status,
    expected_gradient_status
  )


  # ============================================================
  # 12. gradient_ok must represent STRICT gradient convergence
  # ============================================================

  expect_identical(
    conv$gradient_ok,
    is.finite(gmax) &&
      gmax <= tol
  )


  # ============================================================
  # 13. Check sdreport
  # ============================================================

  expect_true(
    conv$sdreport_ok
  )

  expect_false(
    is.null(fit$sdreport)
  )


  # ============================================================
  # 14. Hessian diagnostic must agree with TMB
  # ============================================================

  expect_identical(
    conv$pd_hessian,
    isTRUE(fit$sdreport$pdHess)
  )


  # ============================================================
  # 15. Verify finite standard errors manually
  # ============================================================

  fixed <- fit$fixed_effects

  expect_true(
    "Std. Error" %in% names(fixed)
  )

  manual_finite_se <- all(
    is.finite(
      fixed[["Std. Error"]]
    )
  )

  expect_identical(
    conv$finite_standard_errors,
    manual_finite_se
  )


  # ============================================================
  # 16. Verify strict convergence definition
  # ============================================================

  expected_strict <-
    conv$optimizer_ok &&
    conv$gradient_status == "PASS" &&
    conv$sdreport_ok &&
    conv$pd_hessian &&
    conv$finite_standard_errors

  expect_identical(
    conv$strict_convergence,
    expected_strict
  )


  # ============================================================
  # 17. Verify acceptable convergence definition
  # ============================================================

  expected_acceptable <-
    conv$optimizer_ok &&
    conv$gradient_status %in%
      c("PASS", "MARGINAL") &&
    conv$sdreport_ok &&
    conv$pd_hessian &&
    conv$finite_standard_errors

  expect_identical(
    conv$acceptable_convergence,
    expected_acceptable
  )


  # ============================================================
  # 18. Backward-compatible converged field
  # ============================================================

  expect_identical(
    conv$converged,
    conv$acceptable_convergence
  )


  # ============================================================
  # 19. Overall PASS / MARGINAL / FAIL status
  # ============================================================

  expected_overall <- if (
    conv$strict_convergence
  ) {

    "PASS"

  } else if (
    conv$acceptable_convergence
  ) {

    "MARGINAL"

  } else {

    "FAIL"
  }

  expect_identical(
    conv$overall_status,
    expected_overall
  )


  # ============================================================
  # 20. Test optimization history
  # ============================================================

  history <- conv$optimization_history

  expect_true(
    is.list(history)
  )

  expect_gte(
    length(history),
    1
  )

  expect_lte(
    length(history),
    3
  )

  expect_identical(
    conv$n_optimization_passes,
    length(history)
  )


  # ============================================================
  # 21. selected_pass must be valid
  # ============================================================

  expect_true(
    conv$selected_pass >= 1L
  )

  expect_true(
    conv$selected_pass <= length(history)
  )


  # ============================================================
  # 22. Verify selected pass is the best equivalent solution
  #
  # This is particularly important because the reviewer observed
  # gradients such as:
  #
  #      0.018 -> 0.003 -> 0.004
  #
  # with effectively identical likelihoods.
  # ============================================================

  objectives <- vapply(
    history,
    function(x) x$objective,
    numeric(1)
  )

  gradients <- vapply(
    history,
    function(x) x$max_abs_gradient,
    numeric(1)
  )

  min_obj <- min(
    objectives,
    na.rm = TRUE
  )

  equivalent <- vapply(
    seq_along(objectives),

    function(i) {

      obj_tol <-
        conv$objective_rel_tol *
        max(
          1,
          abs(objectives[i]),
          abs(min_obj)
        )

      abs(
        objectives[i] - min_obj
      ) <= obj_tol
    },

    logical(1)
  )

  candidate_passes <- which(
    equivalent &
      is.finite(gradients)
  )

  expected_best_pass <-
    candidate_passes[
      which.min(
        gradients[candidate_passes]
      )
    ]

  expect_identical(
    conv$selected_pass,
    expected_best_pass
  )


  # ============================================================
  # 23. Returned fit must actually correspond to selected pass
  # ============================================================

  expect_equal(
    fit$fit$objective,
    history[[conv$selected_pass]]$objective,
    tolerance = 1e-8
  )

  expect_equal(
    fit$fit$par,
    history[[conv$selected_pass]]$par,
    tolerance = 1e-8
  )


  # ============================================================
  # 24. Memory settings
  # ============================================================

  expect_false(
    fit$sdreport_settings$getReportCovariance
  )

  expect_false(
    fit$sdreport_settings$getJointPrecision
  )


  # ============================================================
  # 25. Poisson should not use zero inflation
  # ============================================================

  expect_false(
    fit$zero_inflation_structure$enabled
  )

})
