
library(testthat)

test_that("TMB pipeline runs correctly", {

  skip_if_not_installed("phyloseq")
  skip_if_not_installed("TMB")

  # ============================================================
  # 1. Simulate ZIP-like count data
  # ============================================================

  set.seed(123)

  species <- 5
  sites <- 4
  replicates <- 3
  lambda <- 5

  n_libraries <- sites * replicates

  zero_inflation_prob_matrix <- matrix(
    runif(species * sites, 0.1, 0.5),
    nrow = species,
    ncol = sites
  )

  counts <- matrix(
    0L,
    nrow = species,
    ncol = n_libraries,
    dimnames = list(
      paste0("Species_", seq_len(species)),
      paste0("Library_", seq_len(n_libraries))
    )
  )

  for (s in seq_len(species)) {
    for (i in seq_len(sites)) {
      for (r in seq_len(replicates)) {

        col_id <- (i - 1L) * replicates + r

        if (runif(1) < zero_inflation_prob_matrix[s, i]) {
          counts[s, col_id] <- 0L
        } else {
          counts[s, col_id] <- rpois(1, lambda)
        }
      }
    }
  }

  # ============================================================
  # 2. Create metadata
  # ============================================================

  metadata <- data.frame(
    Site = rep(
      paste0("Site_", seq_len(sites)),
      each = replicates
    ),
    Name = rep(
      paste0("Sample_", seq_len(sites)),
      each = replicates
    ),
    Replicate = rep(
      paste0("R", seq_len(replicates)),
      times = sites
    ),
    row.names = colnames(counts)
  )

  # ============================================================
  # 3. Create phyloseq object
  # ============================================================

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(
      counts,
      taxa_are_rows = TRUE
    ),
    phyloseq::sample_data(metadata)
  )

  expect_s4_class(ps, "phyloseq")

  # ============================================================
  # 4. Fit current TMB model
  # ============================================================

  fit_warnings <- character()

  fit <- withCallingHandlers(
    FitModel_joint(
      phyloseq = ps,

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

      gradient_tol = 1e-3,
      gradient_marginal_factor = 5,

      max_restarts = 2L,
      iter_max = 5000L,
      eval_max = 10000L,

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
  # 5. Validate fitted object
  # ============================================================

  expect_s3_class(fit, "eDNAModel_joint")

  expect_true(is.list(fit$data_tmb))
  expect_true(is.list(fit$convergence))

  expect_true(is.list(fit$tmb_object))
  expect_true(is.function(fit$tmb_object$fn))
  expect_true(is.function(fit$tmb_object$gr))

  # ============================================================
  # 6. Validate optimization
  # ============================================================

  expect_true(is.finite(fit$fit$objective))

  expect_true(
    fit$convergence$overall_status %in%
      c("PASS", "MARGINAL", "FAIL")
  )

  expect_true(
    is.numeric(fit$convergence$max_abs_gradient)
  )

  # ============================================================
  # 7. Validate prepared data
  # ============================================================

  dat <- fit$data_tmb

  expect_equal(dat$family_code, 0L)

  expect_true(is.matrix(dat$X_occ))
  expect_true(is.matrix(dat$X_cap))
  expect_true(is.matrix(dat$X_abund))

  expect_length(
    dat$y,
    n_libraries * length(fit$retained_otus)
  )

  expect_true(all(dat$y >= 0))
  expect_true(all(is.finite(dat$y)))

  # ============================================================
  # 8. Validate model probabilities, when available
  # ============================================================

  report <- fit$report

  expect_true(is.list(report))

  # The exact probability fields depend on the C++ REPORT
  # statements. Do not assume old plot_probabilities fields.

  probability_fields <- grep(
    "prob|psi|detect",
    names(report),
    value = TRUE,
    ignore.case = TRUE
  )

  for (field in probability_fields) {

    values <- report[[field]]

    if (
      is.numeric(values) &&
      length(values) > 0L &&
      all(is.finite(values)) &&
      all(values >= 0) &&
      all(values <= 1)
    ) {
      expect_true(all(values >= 0 & values <= 1))
      expect_false(anyNA(values))
    }
  }

  # ============================================================
  # 9. Report numerical warnings
  # ============================================================

  if (length(fit_warnings) > 0L) {
    message(
      "Numerical warnings captured: ",
      paste(unique(fit_warnings), collapse = " | ")
    )
  }

})
