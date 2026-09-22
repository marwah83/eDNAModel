test_that(
  "FitModel_gllvm runs and returns expected structure with synthetic data",
  {

    skip_if_not_installed("gllvm")
    skip_if_not_installed("glmmTMB")
    skip_if_not_installed("phyloseq")

    set.seed(123)

    # ==========================================================
    # 1. Synthetic data
    # ==========================================================

    n_sites <- 12
    n_rep   <- 3
    n_otu   <- 5

    n_samples <- n_sites * n_rep

    # ----------------------------------------------------------
    # Sample metadata
    # ----------------------------------------------------------

    sample_df <- data.frame(

      Site = rep(
        paste0("Loc", seq_len(n_sites)),
        each = n_rep
      ),

      Name = paste0(
        "Sample",
        seq_len(n_samples)
      ),

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
    # Generate heterogeneous OTU counts
    # ----------------------------------------------------------

    # Different baseline abundance among OTUs
    lambda_otu <- c(
      0.8,
      1.2,
      1.8,
      2.5,
      3.5
    )

    species_mat <- matrix(
      0,
      nrow = n_otu,
      ncol = n_samples
    )

    for (k in seq_len(n_otu)) {

      species_mat[k, ] <- rpois(
        n_samples,
        lambda = lambda_otu[k]
      )

    }

    # ----------------------------------------------------------
    # Introduce additional zeros
    #
    # This makes the data more realistic for an
    # occupancy/detection model.
    # ----------------------------------------------------------

    zero_mask <- matrix(
      runif(n_otu * n_samples) < 0.35,
      nrow = n_otu,
      ncol = n_samples
    )

    species_mat[zero_mask] <- 0

    # ----------------------------------------------------------
    # Guarantee every OTU has enough information
    # ----------------------------------------------------------

    for (k in seq_len(n_otu)) {

      if (sum(species_mat[k, ]) == 0) {

        species_mat[
          k,
          sample(
            seq_len(n_samples),
            3
          )
        ] <- c(1, 2, 1)

      }

      if (sum(species_mat[k, ] > 0) < 3) {

        idx <- sample(
          which(species_mat[k, ] == 0),
          3
        )

        species_mat[k, idx] <- c(
          1,
          2,
          1
        )
      }
    }

    rownames(species_mat) <- paste0(
      "OTU",
      seq_len(n_otu)
    )

    colnames(species_mat) <- rownames(
      sample_df
    )

    # ==========================================================
    # 2. Construct phyloseq object
    # ==========================================================

    otu_tab <- phyloseq::otu_table(
      species_mat,
      taxa_are_rows = TRUE
    )

    physeq <- phyloseq::phyloseq(

      otu_tab,

      phyloseq::sample_data(
        sample_df
      )
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
    # 4. Run FitModel_gllvm
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

        abundance_family =
          "poisson",

        min_species_sum =
          1,

        min_detection_replicates =
          1,

        # ----------------------------------------------
        # More than 3 iterations for a stochastic fit
        # ----------------------------------------------

        n_iter =
          10,

        burn_in =
          2,

        num_lv_c =
          1,

        verbose =
          FALSE
      )
    )

    # ==========================================================
    # 5. Expected output components
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
      "note"
    )

    expect_true(
      all(
        expected_components %in%
          names(out)
      )
    )

    # ==========================================================
    # 6. Basic output checks
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
    # 7. Latent-variable outputs
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
    # 8. Iteration outputs
    # ==========================================================
    #
    # IMPORTANT:
    #
    # Do NOT require exactly n_iter - burn_in here.
    #
    # FitModel_gllvm can discard failed iterations.
    # Therefore the correct unit-test condition is that
    # at least one successful post-burn-in iteration exists.
    # ==========================================================

    expect_gt(
      length(out$psi_list),
      0
    )

    expect_gt(
      length(out$capture_list),
      0
    )

    expect_gt(
      length(out$lambda_list),
      0
    )

    expect_gt(
      length(out$p_detect_list),
      0
    )

    # The retained iteration lists should agree in length

    expect_equal(
      length(out$psi_list),
      length(out$capture_list)
    )

    expect_equal(
      length(out$psi_list),
      length(out$lambda_list)
    )

    expect_equal(
      length(out$psi_list),
      length(out$p_detect_list)
    )

    # ==========================================================
    # 9. Summary columns
    # ==========================================================

    expect_true(
      all(
        c(
          "psi_mean",
          "lambda_mean",
          "p_detect_mean"
        ) %in%
          names(out$summary)
      )
    )

    expect_true(
      all(
        c(
          "capture_mean",
          "capture_median",
          "capture_lwr",
          "capture_upr"
        ) %in%
          names(out$capture)
      )
    )

    expect_true(
      all(
        c(
          "capture_mean",
          "capture_median",
          "capture_lwr",
          "capture_upr"
        ) %in%
          names(out$capture_site)
      )
    )

    # ==========================================================
    # 10. Probability sanity checks
    # ==========================================================

    expect_true(
      all(
        out$summary$psi_mean >= 0 &
          out$summary$psi_mean <= 1,
        na.rm = TRUE
      )
    )

    expect_true(
      all(
        out$summary$p_detect_mean >= 0 &
          out$summary$p_detect_mean <= 1,
        na.rm = TRUE
      )
    )

    expect_true(
      all(
        out$summary$lambda_mean >= 0,
        na.rm = TRUE
      )
    )

    # ==========================================================
    # 11. Capture sanity checks
    # ==========================================================

    expect_true(
      all(
        out$capture$capture_mean >= 0 &
          out$capture$capture_mean <= 1,
        na.rm = TRUE
      )
    )

    expect_true(
      all(
        out$capture_site$capture_mean >= 0 &
          out$capture_site$capture_mean <= 1,
        na.rm = TRUE
      )
    )

    # ==========================================================
    # 12. Diagnostic AIC
    # ==========================================================

    expect_s3_class(
      out$diagnostic_AIC,
      "data.frame"
    )

    expect_gte(
      nrow(out$diagnostic_AIC),
      1
    )

  }
)
