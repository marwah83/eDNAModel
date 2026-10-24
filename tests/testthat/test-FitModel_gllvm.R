
test_that(
  "FitModel_gllvm returns valid outputs and retains capture failures",
  {

    skip_if_not_installed("gllvm")
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("phyloseq")

    set.seed(123)

    # ==========================================================
    # 1. Generate hierarchical synthetic data
    # ==========================================================

    n_sites <- 20
    n_rep <- 4
    n_otu <- 5

    n_samples <- n_sites * n_rep

    site_ids <- paste0("Loc", seq_len(n_sites))
    otu_ids <- paste0("OTU", seq_len(n_otu))

    # ----------------------------------------------------------
    # Sample metadata
    # ----------------------------------------------------------

    sample_df <- data.frame(
      Site = rep(site_ids, each = n_rep),
      Name = paste0("Sample", seq_len(n_samples)),
      Replicate = rep(
        seq_len(n_rep),
        times = n_sites
      ),
      row.names = paste0(
        "S",
        seq_len(n_samples)
      ),
      stringsAsFactors = FALSE
    )

    # ----------------------------------------------------------
    # Site-level occupancy probabilities
    # ----------------------------------------------------------

    psi_otu <- c(
      0.55,
      0.65,
      0.70,
      0.75,
      0.80
    )

    # ----------------------------------------------------------
    # Replicate-level detection probabilities
    # ----------------------------------------------------------

    p_otu <- c(
      0.35,
      0.45,
      0.50,
      0.55,
      0.60
    )

    # ----------------------------------------------------------
    # OTU-specific abundance parameters
    # ----------------------------------------------------------

    lambda_otu <- c(
      1.5,
      2.0,
      2.5,
      3.0,
      3.5
    )

    species_mat <- matrix(
      0L,
      nrow = n_otu,
      ncol = n_samples,
      dimnames = list(
        otu_ids,
        rownames(sample_df)
      )
    )

    # ----------------------------------------------------------
    # Generate occupancy, detection and abundance
    # ----------------------------------------------------------

    for (k in seq_len(n_otu)) {

      occupied <- rbinom(
        n_sites,
        size = 1,
        prob = psi_otu[k]
      )

      for (s in seq_len(n_sites)) {

        idx <- which(
          sample_df$Site == site_ids[s]
        )

        if (occupied[s] == 1) {

          detected <- rbinom(
            n_rep,
            size = 1,
            prob = p_otu[k]
          )

          positive_counts <- 1L + rpois(
            n_rep,
            lambda = lambda_otu[k]
          )

          species_mat[k, idx] <-
            detected * positive_counts
        }
      }
    }

    # ----------------------------------------------------------
    # Guarantee informative detection histories
    #
    # Each OTU has at least one site containing both
    # detections and non-detections.
    # ----------------------------------------------------------

    for (k in seq_len(n_otu)) {

      idx <- which(
        sample_df$Site == site_ids[k]
      )

      species_mat[k, idx] <- c(
        2L,
        0L,
        3L,
        0L
      )
    }

    # ==========================================================
    # 2. Construct phyloseq object
    # ==========================================================

    otu_tab <- phyloseq::otu_table(
      species_mat,
      taxa_are_rows = TRUE
    )

    physeq <- phyloseq::phyloseq(
      otu_tab,
      phyloseq::sample_data(sample_df)
    )

    # ==========================================================
    # 3. Check synthetic data before fitting
    # ==========================================================

    expect_equal(
      phyloseq::nsamples(physeq),
      n_samples
    )

    expect_equal(
      phyloseq::ntaxa(physeq),
      n_otu
    )

    expect_true(
      any(species_mat == 0)
    )

    expect_true(
      any(species_mat > 0)
    )

    expect_true(
      all(rowSums(species_mat) > 0)
    )

    # ==========================================================
    # 4. Check replicated detection histories
    # ==========================================================

    observed <- data.frame(
      Site = rep(
        sample_df$Site,
        each = n_otu
      ),
      OTU = rep(
        otu_ids,
        times = n_samples
      ),
      detected = as.integer(
        as.vector(species_mat > 0)
      )
    )

    histories <- observed %>%
      dplyr::group_by(Site, OTU) %>%
      dplyr::summarise(
        n_detected = sum(detected),
        n_failed = sum(detected == 0),
        .groups = "drop"
      )

    confirmed <- histories %>%
      dplyr::filter(n_detected > 0)

    expect_gt(
      nrow(confirmed),
      0
    )

    expect_gt(
      sum(confirmed$n_failed),
      0
    )

    expect_true(
      any(
        confirmed$n_detected > 0 &
          confirmed$n_failed > 0
      )
    )

    # ==========================================================
    # 5. Run FitModel_gllvm
    # ==========================================================

    out <- suppressWarnings(
      FitModel_gllvm(
        phyloseq = physeq,

        site_col = "Site",

        otu_col = "OTU",
        count_col = "y",

        sample_col = "Name",
        replicate_col = "Replicate",

        abundance_rhs =
          y ~ (1 | OTU),

        capture_formula =
          a_sim ~ 1 + (1 | OTU),

        occupancy_covars = NULL,

        abundance_family = "poisson",

        min_species_sum = 1,
        min_detection_replicates = 1,

        n_iter = 10,
        burn_in = 2,

        num_lv_c = 1,

        gllvm_auto_exclude = TRUE,
        gllvm_n_init = 5,

        verbose = FALSE
      )
    )

    # ==========================================================
    # 6. Expected output components
    # ==========================================================

    expected_components <- c(
      "summary",
      "capture",
      "capture_site",

      "psi_list",
      "capture_list",
      "lambda_list",
      "p_detect_list",

      "occupancy_models",
      "capture_models",
      "abundance_models",

      "reduced_data",
      "sample_data",
      "long_df",

      "lv_sites",
      "lv_species",
      "mean_lv_sites",
      "mean_lv_species",

      "filter_summary",
      "diagnostic_AIC",

      "gllvm_diagnostics",
      "gllvm_excluded_otus",
      "retained_iterations",

      "note"
    )

    expect_true(
      all(
        expected_components %in% names(out)
      )
    )

    # ==========================================================
    # 7. Basic output checks
    # ==========================================================

    expect_s3_class(
      out$summary,
      "data.frame"
    )

    expect_s3_class(
      out$capture,
      "data.frame"
    )

    expect_s3_class(
      out$capture_site,
      "data.frame"
    )

    expect_gt(
      nrow(out$summary),
      0
    )

    # ==========================================================
    # 8. Latent-variable outputs
    # ==========================================================

    expect_s3_class(
      out$lv_sites,
      "data.frame"
    )

    expect_s3_class(
      out$lv_species,
      "data.frame"
    )

    expect_s3_class(
      out$mean_lv_sites,
      "data.frame"
    )

    expect_s3_class(
      out$mean_lv_species,
      "data.frame"
    )

    # ==========================================================
    # 9. Iteration consistency
    # ==========================================================

    n_keep <- length(out$psi_list)

    expect_gt(
      n_keep,
      0
    )

    expect_lte(
      n_keep,
      8
    )

    iteration_components <- c(
      "capture_list",
      "lambda_list",
      "p_detect_list",
      "capture_models",
      "abundance_models"
    )

    for (nm in iteration_components) {

      expect_equal(
        length(out[[nm]]),
        n_keep
      )
    }

    # ==========================================================
    # 10. Required summary columns
    #
    # FIX:
    # Check that required columns are present without
    # requiring an exact match of all column names.
    # ==========================================================

    required_summary_cols <- c(
      "psi_mean",
      "lambda_mean",
      "p_detect_mean"
    )

    expect_true(
      all(
        required_summary_cols %in%
          names(out$summary)
      ),
      info = "Missing required columns in out$summary"
    )

    capture_cols <- c(
      "capture_mean",
      "capture_median",
      "capture_lwr",
      "capture_upr"
    )

    expect_true(
      all(
        capture_cols %in%
          names(out$capture)
      )
    )

    expect_true(
      all(
        capture_cols %in%
          names(out$capture_site)
      )
    )

    # ==========================================================
    # 11. Probability sanity checks
    # ==========================================================

    check_probability <- function(x) {

      expect_true(
        is.numeric(x)
      )

      expect_true(
        all(is.finite(x))
      )

      expect_true(
        all(x >= 0 & x <= 1)
      )
    }

    check_probability(
      out$summary$psi_mean
    )

    check_probability(
      out$summary$p_detect_mean
    )

    expect_true(
      all(
        is.finite(out$summary$lambda_mean) &
          out$summary$lambda_mean >= 0
      )
    )

    # ==========================================================
    # 12. Capture probability checks
    # ==========================================================

    check_probability(
      out$capture$capture_mean
    )

    check_probability(
      out$capture_site$capture_mean
    )

    # ==========================================================
    # 13. Diagnostic AIC
    # ==========================================================

    expect_s3_class(
      out$diagnostic_AIC,
      "data.frame"
    )

    expect_gte(
      nrow(out$diagnostic_AIC),
      1
    )

    # ==========================================================
    # 14. Capture model diagnostics
    # ==========================================================

    capture_links <- vapply(
      out$capture_models,
      function(m) {
        stats::family(m)$link
      },
      character(1)
    )

    expect_true(
      all(capture_links == "logit")
    )

    # ----------------------------------------------------------
    # Count capture failures in each retained iteration
    # ----------------------------------------------------------

    capture_counts <- lapply(
      seq_along(out$capture_models),
      function(i) {

        m <- out$capture_models[[i]]

        y_capture <- stats::model.response(
          stats::model.frame(m)
        )

        data.frame(
          iteration = i,
          n_zero = sum(y_capture == 0),
          n_one = sum(y_capture == 1),
          prop_one = mean(y_capture == 1)
        )
      }
    )

    capture_counts <- dplyr::bind_rows(
      capture_counts
    )

    # At least one iteration should contain capture failures.
    #
    # This is a minimal regression check, not proof that
    # the capture component is statistically identifiable.

    expect_true(
      any(capture_counts$n_zero > 0),
      info = paste(
        "No capture failures entered any retained",
        "capture model despite informative",
        "replicated detection histories."
      )
    )

    expect_true(
      all(capture_counts$n_one > 0)
    )

    # ----------------------------------------------------------
    # Report diagnostic information
    # ----------------------------------------------------------

    message(
      "Capture-model diagnostic summary:"
    )

    print(capture_counts)

    message(
      "Proportion of capture models with no failures: ",
      round(
        mean(capture_counts$n_zero == 0),
        3
      )
    )

    # ==========================================================
    # 15. Prediction-scale checks
    # ==========================================================

    # ----------------------------------------------------------
    # Capture GLMM: logit link
    # ----------------------------------------------------------

    capture_draws <- dplyr::bind_rows(
      out$capture_list
    )

    expect_true(
      "eta" %in% names(capture_draws)
    )

    capture_eta <- capture_draws$eta

    expect_true(
      all(is.finite(capture_eta))
    )

    check_probability(
      stats::plogis(capture_eta)
    )

    # ----------------------------------------------------------
    # Detection from p_detect_list:
    # inverse complementary log-log transformation
    # ----------------------------------------------------------

    p_detect_draws <- dplyr::bind_rows(
      out$p_detect_list
    )

    expect_true(
      "eta" %in% names(p_detect_draws)
    )

    p_eta <- p_detect_draws$eta

    expect_true(
      all(is.finite(p_eta))
    )

    p_detect <- -expm1(
      -exp(p_eta)
    )

    check_probability(
      p_detect
    )

  }
)
