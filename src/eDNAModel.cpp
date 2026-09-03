Pasted text(20260903-153910).txt
Document
can you update the FitModel_Joint with the random_Zi_otu

Yes. Below is the full corrected FitModel_joint() based directly on the function you uploaded. Your original function already correctly separates occupancy at site × OTU, capture at sample × OTU, and abundance at PCR/read level.

This version adds random_zi_otu = TRUE, passes use_zi_otu to C++, adds b_zi_otu and log_sd_zi_otu, integrates b_zi_otu through TMB, maps it out when appropriate, retains your PCR-replication requirement for ZIP/ZINB, and updates the output descriptions.

Important: this R function assumes you are using the updated C++ code containing DATA_INTEGER(use_zi_otu), PARAMETER_VECTOR(b_zi_otu), and PARAMETER(log_sd_zi_otu).

FitModel_joint <- function(
    phyloseq,
    site_col,
    sample_col = "Name",
    replicate_col = NULL,
    otu_col = "OTU",
    count_col = "y",

    occupancy_formula = ~ 1,
    capture_formula = ~ 1,
    abundance_formula = ~ 1,

    abundance_offset = NULL,

    abundance_family = c(
        "poisson",
        "nbinom",
        "zip",
        "zinb"
    ),

    min_species_sum = 10,
    min_detection_replicates = 1,

    random_occ_otu = TRUE,
    random_capture_otu = TRUE,
    random_abund_otu = TRUE,
    random_sample = TRUE,
    random_sample_otu = FALSE,

    # NEW:
    # OTU-specific random intercept for ZIP/ZINB
    # structural-zero probability.
    random_zi_otu = TRUE,

    DLL = "eDNAModel",

    verbose = TRUE
) {

    # ============================================================
    # 0. Match abundance family
    # ============================================================

    abundance_family <-
        match.arg(abundance_family)


    # ------------------------------------------------------------
    # random_zi_otu only makes sense for ZIP/ZINB
    # ------------------------------------------------------------

    if (
        random_zi_otu &&
        !abundance_family %in% c("zip", "zinb")
    ) {

        if (verbose) {
            message(
                "random_zi_otu is ignored because abundance_family = '",
                abundance_family,
                "'."
            )
        }

        random_zi_otu <- FALSE
    }


    # ============================================================
    # 1. Prepare long data
    # ============================================================

    prep <- prepare_long_data(
        physeq_obj = phyloseq,
        site_col = site_col,
        nested_cols = unique(
            stats::na.omit(
                c(
                    sample_col,
                    replicate_col
                )
            )
        )
    )

    dat <- prep$long_df


    # ============================================================
    # 2. Validate required columns
    # ============================================================

    required <- unique(
        c(
            site_col,
            sample_col,
            otu_col,
            count_col
        )
    )

    if (!is.null(replicate_col)) {

        required <- unique(
            c(
                required,
                replicate_col
            )
        )
    }

    missing_cols <-
        setdiff(
            required,
            names(dat)
        )

    if (length(missing_cols) > 0) {

        stop(
            "Missing columns: ",
            paste(
                missing_cols,
                collapse = ", "
            )
        )
    }


    # ============================================================
    # 3. Type handling
    # ============================================================

    dat[[site_col]] <-
        as.character(
            dat[[site_col]]
        )

    dat[[sample_col]] <-
        as.character(
            dat[[sample_col]]
        )

    dat[[otu_col]] <-
        as.character(
            dat[[otu_col]]
        )

    dat[[count_col]] <-
        as.numeric(
            dat[[count_col]]
        )

    if (!is.null(replicate_col)) {

        dat[[replicate_col]] <-
            as.character(
                dat[[replicate_col]]
            )
    }


    # Remove missing/non-finite counts

    dat <- dat[
        !is.na(dat[[count_col]]) &
            is.finite(dat[[count_col]]),
        ,
        drop = FALSE
    ]


    if (
        any(
            dat[[count_col]] < 0
        )
    ) {

        stop(
            "Counts must be non-negative."
        )
    }


    # ============================================================
    # 4. ZIP / ZINB IDENTIFIABILITY CHECK
    # ============================================================
    #
    # ZIP/ZINB introduce a structural-zero process at the
    # PCR/read level.
    #
    # Hierarchy:
    #
    #   Z_im
    #       -> A_ijm
    #           -> Y_ijkm
    #
    # A_ijm already provides a sample-level zero-generating
    # process.
    #
    # Therefore, without repeated PCR observations within a
    # biological sample, the capture-zero process and the
    # read-level zero-inflation process cannot be cleanly
    # separated.
    #
    # ZIP/ZINB are therefore allowed only when explicit
    # PCR replication is present.
    # ============================================================

    replication_summary <- NULL


    if (
        abundance_family %in%
        c(
            "zip",
            "zinb"
        )
    ) {

        # --------------------------------------------------------
        # replicate_col must be supplied
        # --------------------------------------------------------

        if (is.null(replicate_col)) {

            stop(
                abundance_family,
                " requires PCR-level replication, but ",
                "`replicate_col` is NULL. ",
                "Without repeated PCR/read observations within ",
                "biological samples, read-level zero inflation ",
                "is confounded with the sample-level capture process. ",
                "Use abundance_family = 'poisson' or 'nbinom', ",
                "or provide a PCR replicate column."
            )
        }


        # --------------------------------------------------------
        # Replicate IDs must not be missing
        # --------------------------------------------------------

        if (
            any(
                is.na(
                    dat[[replicate_col]]
                ) |
                    dat[[replicate_col]] == ""
            )
        ) {

            stop(
                abundance_family,
                " requires valid PCR replicate identifiers. ",
                "Missing or empty values were found in `",
                replicate_col,
                "`."
            )
        }


        # --------------------------------------------------------
        # Count distinct PCR replicates within each biological
        # sample.
        #
        # distinct() is used because long-format data contain
        # one row per OTU as well as replicate.
        # --------------------------------------------------------

        replication_summary <-
            dat |>
            dplyr::distinct(
                dplyr::across(
                    dplyr::all_of(
                        c(
                            site_col,
                            sample_col,
                            replicate_col
                        )
                    )
                )
            ) |>
            dplyr::group_by(
                dplyr::across(
                    dplyr::all_of(
                        c(
                            site_col,
                            sample_col
                        )
                    )
                )
            ) |>
            dplyr::summarise(
                n_pcr_replicates =
                    dplyr::n_distinct(
                        .data[[replicate_col]]
                    ),
                .groups = "drop"
            )


        # --------------------------------------------------------
        # Require >= 2 PCR replicates for every biological sample
        # --------------------------------------------------------

        insufficient_replication <-
            replication_summary |>
            dplyr::filter(
                .data$n_pcr_replicates < 2
            )


        if (
            nrow(
                insufficient_replication
            ) > 0
        ) {

            stop(
                abundance_family,
                " requires at least two PCR replicates per ",
                "biological sample. ",
                nrow(insufficient_replication),
                " biological sample(s) have fewer than two ",
                "PCR replicates. ",
                "Without replication, the ZIP/ZINB structural-zero ",
                "probability cannot be cleanly separated from ",
                "sample-level capture failure. ",
                "Use abundance_family = 'poisson' or 'nbinom', ",
                "or provide replicated PCR observations."
            )
        }


        if (verbose) {

            message(
                "PCR replication check passed for ",
                toupper(abundance_family),
                "."
            )

            message(
                "PCR replicates per biological sample: min = ",
                min(
                    replication_summary$n_pcr_replicates
                ),
                ", median = ",
                stats::median(
                    replication_summary$n_pcr_replicates
                ),
                ", max = ",
                max(
                    replication_summary$n_pcr_replicates
                )
            )
        }
    }


    # ============================================================
    # 5. OTU filtering
    # ============================================================

    otu_stats <- dat |>
        dplyr::group_by(
            .data[[otu_col]]
        ) |>
        dplyr::summarise(

            total_count =
                sum(
                    .data[[count_col]],
                    na.rm = TRUE
                ),

            detected_replicates =
                sum(
                    .data[[count_col]] > 0,
                    na.rm = TRUE
                ),

            .groups = "drop"
        )


    retained_otus <-
        otu_stats |>
        dplyr::filter(

            .data$total_count >=
                min_species_sum,

            .data$detected_replicates >=
                min_detection_replicates
        ) |>
        dplyr::pull(
            .data[[otu_col]]
        )


    dat <-
        dat |>
        dplyr::filter(
            .data[[otu_col]] %in%
                retained_otus
        )


    if (nrow(dat) == 0) {

        stop(
            "No OTUs remain after filtering."
        )
    }


    # ============================================================
    # 6. Explicit factor indices
    # ============================================================

    # Recreate factors AFTER filtering.
    #
    # This is important because all OTU-indexed random effects:
    #
    #   b_occ_otu
    #   b_cap_otu
    #   b_abund_otu
    #   b_zi_otu
    #
    # must use exactly the same OTU indexing.

    dat$.site <-
        factor(
            dat[[site_col]]
        )

    dat$.otu <-
        factor(
            dat[[otu_col]]
        )

    dat$.sample <-
        factor(
            dat[[sample_col]]
        )


    # ------------------------------------------------------------
    # SITE x OTU
    #
    # Occupancy level
    # ------------------------------------------------------------

    dat$.site_otu <-
        interaction(
            dat$.site,
            dat$.otu,
            drop = TRUE,
            lex.order = TRUE
        )


    # ------------------------------------------------------------
    # SITE x SAMPLE x OTU
    #
    # Capture level
    # ------------------------------------------------------------

    dat$.sample_otu <-
        interaction(
            dat$.site,
            dat$.sample,
            dat$.otu,
            drop = TRUE,
            lex.order = TRUE
        )


    # ============================================================
    # 7. Build SAMPLE x OTU table
    # ============================================================

    sample_df <-
        dat |>
        dplyr::group_by(
            .data$.sample_otu
        ) |>
        dplyr::summarise(

            .site_otu =
                dplyr::first(
                    .data$.site_otu
                ),

            .otu =
                dplyr::first(
                    .data$.otu
                ),

            .sample =
                dplyr::first(
                    .data$.sample
                ),

            sample_positive =
                as.integer(
                    any(
                        .data[[count_col]] > 0
                    )
                ),

            .groups = "drop"
        )


    # ============================================================
    # 8. Carry capture-formula variables
    # ============================================================

    cap_vars <-
        all.vars(
            capture_formula
        )

    cap_vars <-
        intersect(
            cap_vars,
            names(dat)
        )


    if (length(cap_vars) > 0) {

        cap_covars <-
            dat |>
            dplyr::select(
                ".sample_otu",
                dplyr::all_of(
                    cap_vars
                )
            ) |>
            dplyr::group_by(
                .data$.sample_otu
            ) |>
            dplyr::summarise(
                dplyr::across(
                    dplyr::all_of(
                        cap_vars
                    ),
                    ~ dplyr::first(.x)
                ),
                .groups = "drop"
            )


        sample_df <-
            dplyr::left_join(
                sample_df,
                cap_covars,
                by = ".sample_otu"
            )
    }


    # ============================================================
    # 9. Build SITE x OTU table
    # ============================================================

    site_df <-
        dat |>
        dplyr::group_by(
            .data$.site_otu
        ) |>
        dplyr::summarise(

            .otu =
                dplyr::first(
                    .data$.otu
                ),

            site_positive =
                as.integer(
                    any(
                        .data[[count_col]] > 0
                    )
                ),

            .groups = "drop"
        )


    # ============================================================
    # 10. Carry occupancy-formula variables
    # ============================================================

    occ_vars <-
        all.vars(
            occupancy_formula
        )

    occ_vars <-
        intersect(
            occ_vars,
            names(dat)
        )


    if (length(occ_vars) > 0) {

        occ_covars <-
            dat |>
            dplyr::select(
                ".site_otu",
                dplyr::all_of(
                    occ_vars
                )
            ) |>
            dplyr::group_by(
                .data$.site_otu
            ) |>
            dplyr::summarise(
                dplyr::across(
                    dplyr::all_of(
                        occ_vars
                    ),
                    ~ dplyr::first(.x)
                ),
                .groups = "drop"
            )


        site_df <-
            dplyr::left_join(
                site_df,
                occ_covars,
                by = ".site_otu"
            )
    }


    # ============================================================
    # 11. Design matrices
    # ============================================================
    #
    # X_occ:
    #   one row per SITE x OTU
    #
    # X_cap:
    #   one row per BIOLOGICAL SAMPLE x OTU
    #
    # X_abund:
    #   one row per PCR/read observation
    #
    # Therefore the three ecological/observation processes live
    # at their appropriate hierarchical levels.
    # ============================================================

    X_occ <-
        stats::model.matrix(
            occupancy_formula,
            data = site_df
        )

    X_cap <-
        stats::model.matrix(
            capture_formula,
            data = sample_df
        )

    X_abund <-
        stats::model.matrix(
            abundance_formula,
            data = dat
        )


    # ============================================================
    # 12. Abundance offset
    # ============================================================

    if (is.null(abundance_offset)) {

        offset_abund <-
            rep(
                0,
                nrow(dat)
            )

    } else {

        if (
            !(abundance_offset %in%
              names(dat))
        ) {

            stop(
                "abundance_offset '",
                abundance_offset,
                "' not found."
            )
        }


        off <-
            as.numeric(
                dat[[abundance_offset]]
            )


        if (
            any(
                !is.finite(off) |
                    off <= 0
            )
        ) {

            stop(
                "abundance_offset must contain positive finite values."
            )
        }


        offset_abund <-
            log(off)
    }


    # ============================================================
    # 13. Explicit group levels
    # ============================================================

    site_group_levels <-
        levels(
            dat$.site_otu
        )

    sample_group_levels <-
        levels(
            dat$.sample_otu
        )


    # ============================================================
    # 14. Mapping indices
    # ============================================================

    # PCR/read row -> SAMPLE x OTU group

    row_sample_group <-
        match(
            as.character(
                dat$.sample_otu
            ),
            sample_group_levels
        ) -
        1L


    # SAMPLE x OTU -> SITE x OTU

    sample_site_group <-
        match(
            as.character(
                sample_df$.site_otu
            ),
            site_group_levels
        ) -
        1L


    # SITE x OTU -> OTU

    occ_otu <-
        as.integer(
            site_df$.otu
        ) -
        1L


    # SAMPLE x OTU -> OTU

    sample_otu <-
        as.integer(
            sample_df$.otu
        ) -
        1L


    # PCR/read row -> OTU
    #
    # IMPORTANT:
    # This same index is now used by:
    #
    #   b_abund_otu
    #   b_zi_otu

    row_otu <-
        as.integer(
            dat$.otu
        ) -
        1L


    # PCR/read row -> biological sample

    row_sample_id <-
        as.integer(
            dat$.sample
        ) -
        1L


    # One random-effect level per sample x OTU

    row_sample_otu_re <-
        row_sample_group


    # ============================================================
    # 15. Validate mappings
    # ============================================================

    if (
        anyNA(
            row_sample_group
        )
    ) {

        stop(
            "Invalid row_sample_group mapping."
        )
    }


    if (
        anyNA(
            sample_site_group
        )
    ) {

        stop(
            "Invalid sample_site_group mapping."
        )
    }


    if (
        any(
            row_sample_group < 0 |
                row_sample_group >= nrow(sample_df)
        )
    ) {

        stop(
            "row_sample_group contains invalid indices."
        )
    }


    if (
        any(
            sample_site_group < 0 |
                sample_site_group >= nrow(site_df)
        )
    ) {

        stop(
            "sample_site_group contains invalid indices."
        )
    }


    # ============================================================
    # 16. Family code
    # ============================================================

    family_code <-
        switch(

            abundance_family,

            poisson = 0L,

            nbinom = 1L,

            zip = 2L,

            zinb = 3L
        )


    # ============================================================
    # 17. Dimensions
    # ============================================================
    #
    # Calculate these BEFORE constructing data_tmb/parameters.
    # ============================================================

    n_otu <-
        nlevels(
            dat$.otu
        )

    n_sample <-
        nlevels(
            dat$.sample
        )

    n_sample_otu <-
        nrow(
            sample_df
        )


    # ============================================================
    # 18. TMB data list
    # ============================================================

    data_tmb <- list(

        y =
            as.numeric(
                dat[[count_col]]
            ),

        X_occ =
            X_occ,

        X_cap =
            X_cap,

        X_abund =
            X_abund,

        offset_abund =
            as.numeric(
                offset_abund
            ),

        occ_otu =
            as.integer(
                occ_otu
            ),

        sample_otu =
            as.integer(
                sample_otu
            ),

        row_sample_group =
            as.integer(
                row_sample_group
            ),

        sample_site_group =
            as.integer(
                sample_site_group
            ),

        row_otu =
            as.integer(
                row_otu
            ),

        row_sample_id =
            as.integer(
                row_sample_id
            ),

        row_sample_otu_re =
            as.integer(
                row_sample_otu_re
            ),

        sample_positive =
            as.integer(
                sample_df$sample_positive
            ),

        site_positive =
            as.integer(
                site_df$site_positive
            ),

        n_site_groups =
            nrow(
                site_df
            ),

        n_sample_groups =
            nrow(
                sample_df
            ),

        family_code =
            family_code,

        use_occ_otu =
            as.integer(
                random_occ_otu
            ),

        use_cap_otu =
            as.integer(
                random_capture_otu
            ),

        use_abund_otu =
            as.integer(
                random_abund_otu
            ),

        use_sample_re =
            as.integer(
                random_sample
            ),

        use_sample_otu_re =
            as.integer(
                random_sample_otu
            ),

        # --------------------------------------------------------
        # NEW
        # --------------------------------------------------------
        #
        # Tell C++ whether the OTU-specific zero-inflation
        # random intercept is enabled.

        use_zi_otu =
            as.integer(
                random_zi_otu
            )
    )


    # ============================================================
    # 19. Starting values
    # ============================================================

    site_positive_rate <-
        mean(
            site_df$site_positive
        )

    site_positive_rate <-
        pmin(
            pmax(
                site_positive_rate,
                0.05
            ),
            0.95
        )


    sample_positive_rate <-
        mean(
            sample_df$sample_positive
        )

    sample_positive_rate <-
        pmin(
            pmax(
                sample_positive_rate,
                0.05
            ),
            0.95
        )


    positive_counts <-
        dat[[count_col]][
            dat[[count_col]] > 0
        ]


    if (length(positive_counts) > 0) {

        mean_count_start <-
            mean(
                positive_counts
            )

    } else {

        mean_count_start <-
            1
    }


    # ------------------------------------------------------------
    # Occupancy fixed-effect starts
    # ------------------------------------------------------------

    beta_occ_start <-
        rep(
            0,
            ncol(
                X_occ
            )
        )


    occ_intercept <-
        match(
            "(Intercept)",
            colnames(
                X_occ
            )
        )


    if (!is.na(occ_intercept)) {

        beta_occ_start[
            occ_intercept
        ] <-
            stats::qlogis(
                site_positive_rate
            )
    }


    # ------------------------------------------------------------
    # Capture fixed-effect starts
    # ------------------------------------------------------------

    beta_cap_start <-
        rep(
            0,
            ncol(
                X_cap
            )
        )


    cap_intercept <-
        match(
            "(Intercept)",
            colnames(
                X_cap
            )
        )


    if (!is.na(cap_intercept)) {

        beta_cap_start[
            cap_intercept
        ] <-
            stats::qlogis(
                sample_positive_rate
            )
    }


    # ------------------------------------------------------------
    # Abundance fixed-effect starts
    # ------------------------------------------------------------

    beta_abund_start <-
        rep(
            0,
            ncol(
                X_abund
            )
        )


    abund_intercept <-
        match(
            "(Intercept)",
            colnames(
                X_abund
            )
        )


    if (!is.na(abund_intercept)) {

        beta_abund_start[
            abund_intercept
        ] <-
            log(
                pmax(
                    mean_count_start,
                    1e-4
                )
            )
    }


    # ============================================================
    # 20. Full parameter list
    # ============================================================

    parameters <- list(

        # --------------------------------------------------------
        # Fixed effects
        # --------------------------------------------------------

        beta_occ =
            beta_occ_start,

        beta_cap =
            beta_cap_start,

        beta_abund =
            beta_abund_start,


        # --------------------------------------------------------
        # OTU random effects
        # --------------------------------------------------------

        b_occ_otu =
            rep(
                0,
                n_otu
            ),

        b_cap_otu =
            rep(
                0,
                n_otu
            ),

        b_abund_otu =
            rep(
                0,
                n_otu
            ),

        # NEW:
        #
        # OTU-specific zero-inflation deviations:
        #
        #   b_zi_m ~ N(0, sigma_zi^2)

        b_zi_otu =
            rep(
                0,
                n_otu
            ),


        # --------------------------------------------------------
        # Sample-level abundance random effects
        # --------------------------------------------------------

        b_sample =
            rep(
                0,
                n_sample
            ),

        b_sample_otu =
            rep(
                0,
                n_sample_otu
            ),


        # --------------------------------------------------------
        # Log standard deviations
        # --------------------------------------------------------

        log_sd_occ_otu =
            log(0.5),

        log_sd_cap_otu =
            log(0.5),

        log_sd_abund_otu =
            log(0.5),

        # NEW:
        #
        # sigma_zi =
        # exp(log_sd_zi_otu)

        log_sd_zi_otu =
            log(0.5),

        log_sd_sample =
            log(0.5),

        log_sd_sample_otu =
            log(0.5),


        # --------------------------------------------------------
        # NB2 dispersion
        # --------------------------------------------------------

        # Currently one shared dispersion parameter across OTUs.

        log_theta =
            log(10),


        # --------------------------------------------------------
        # Zero-inflation population intercept
        # --------------------------------------------------------
        #
        # With random_zi_otu = TRUE:
        #
        #   logit(pi_m)
        #       =
        #       zi_intercept
        #       +
        #       b_zi_otu[m]
        #
        # With random_zi_otu = FALSE:
        #
        #   logit(pi)
        #       =
        #       zi_intercept

        zi_intercept =
            stats::qlogis(
                0.05
            )
    )


    # ============================================================
    # 21. Map unused parameters out of optimization
    # ============================================================

    map <- list()


    # ------------------------------------------------------------
    # Occupancy OTU random effect OFF
    # ------------------------------------------------------------

    if (!random_occ_otu) {

        map$b_occ_otu <-
            factor(
                rep(
                    NA,
                    n_otu
                )
            )

        map$log_sd_occ_otu <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # Capture OTU random effect OFF
    # ------------------------------------------------------------

    if (!random_capture_otu) {

        map$b_cap_otu <-
            factor(
                rep(
                    NA,
                    n_otu
                )
            )

        map$log_sd_cap_otu <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # Abundance OTU random effect OFF
    # ------------------------------------------------------------

    if (!random_abund_otu) {

        map$b_abund_otu <-
            factor(
                rep(
                    NA,
                    n_otu
                )
            )

        map$log_sd_abund_otu <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # Biological-sample random effect OFF
    # ------------------------------------------------------------

    if (!random_sample) {

        map$b_sample <-
            factor(
                rep(
                    NA,
                    n_sample
                )
            )

        map$log_sd_sample <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # Sample x OTU random effect OFF
    # ------------------------------------------------------------

    if (!random_sample_otu) {

        map$b_sample_otu <-
            factor(
                rep(
                    NA,
                    n_sample_otu
                )
            )

        map$log_sd_sample_otu <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # NEW:
    # OTU-specific zero-inflation random effect OFF
    # ------------------------------------------------------------

    if (!random_zi_otu) {

        map$b_zi_otu <-
            factor(
                rep(
                    NA,
                    n_otu
                )
            )

        map$log_sd_zi_otu <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # Poisson
    #
    # No NB dispersion.
    # No zero inflation.
    # ------------------------------------------------------------

    if (
        abundance_family ==
            "poisson"
    ) {

        map$log_theta <-
            factor(NA)

        map$zi_intercept <-
            factor(NA)

        map$b_zi_otu <-
            factor(
                rep(
                    NA,
                    n_otu
                )
            )

        map$log_sd_zi_otu <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # Negative binomial
    #
    # theta estimated.
    # No zero inflation.
    # ------------------------------------------------------------

    if (
        abundance_family ==
            "nbinom"
    ) {

        map$zi_intercept <-
            factor(NA)

        map$b_zi_otu <-
            factor(
                rep(
                    NA,
                    n_otu
                )
            )

        map$log_sd_zi_otu <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # ZIP
    #
    # Zero inflation estimated.
    # NB theta removed.
    #
    # If random_zi_otu = TRUE:
    #
    #   zi_intercept
    #   b_zi_otu
    #   log_sd_zi_otu
    #
    # are active.
    # ------------------------------------------------------------

    if (
        abundance_family ==
            "zip"
    ) {

        map$log_theta <-
            factor(NA)
    }


    # ------------------------------------------------------------
    # ZINB
    #
    # log_theta remains estimated.
    # zi_intercept remains estimated.
    #
    # If random_zi_otu = TRUE:
    #
    #   b_zi_otu
    #   log_sd_zi_otu
    #
    # also remain active.
    # ------------------------------------------------------------


    # ============================================================
    # 22. Gaussian random effects integrated by TMB
    # ============================================================

    random_effects <-
        character(0)


    if (random_occ_otu) {

        random_effects <-
            c(
                random_effects,
                "b_occ_otu"
            )
    }


    if (random_capture_otu) {

        random_effects <-
            c(
                random_effects,
                "b_cap_otu"
            )
    }


    if (random_abund_otu) {

        random_effects <-
            c(
                random_effects,
                "b_abund_otu"
            )
    }


    if (random_sample) {

        random_effects <-
            c(
                random_effects,
                "b_sample"
            )
    }


    if (random_sample_otu) {

        random_effects <-
            c(
                random_effects,
                "b_sample_otu"
            )
    }


    # ------------------------------------------------------------
    # NEW:
    # Integrate OTU-specific zero-inflation effects using
    # the TMB Laplace approximation.
    # ------------------------------------------------------------

    if (random_zi_otu) {

        random_effects <-
            c(
                random_effects,
                "b_zi_otu"
            )
    }


    # ============================================================
    # 23. Verbose model summary
    # ============================================================

    if (verbose) {

        message(
            "Abundance family: ",
            abundance_family
        )

        message(
            "Sites x OTUs: ",
            nrow(
                site_df
            )
        )

        message(
            "Samples x OTUs: ",
            nrow(
                sample_df
            )
        )

        message(
            "Read observations: ",
            nrow(
                dat
            )
        )

        message(
            "OTUs: ",
            n_otu
        )

        message(
            "Biological samples: ",
            n_sample
        )

        message(
            "Random effects: ",
            if (
                length(
                    random_effects
                ) == 0
            ) {

                "none"

            } else {

                paste(
                    random_effects,
                    collapse = ", "
                )
            }
        )


        if (
            abundance_family %in%
            c(
                "nbinom",
                "zinb"
            )
        ) {

            message(
                "NB2 dispersion structure: ",
                "one shared theta across OTUs."
            )
        }


        if (
            abundance_family %in%
            c(
                "zip",
                "zinb"
            )
        ) {

            if (random_zi_otu) {

                message(
                    "Zero-inflation structure: ",
                    "shared population intercept + ",
                    "OTU-specific Gaussian random intercept."
                )

            } else {

                message(
                    "Zero-inflation structure: ",
                    "one shared pi across OTUs/read observations."
                )
            }
        }
    }


    # ============================================================
    # 24. Make TMB objective
    # ============================================================

    obj <- TMB::MakeADFun(

        data =
            data_tmb,

        parameters =
            parameters,

        random =
            if (
                length(
                    random_effects
                ) == 0
            ) {

                NULL

            } else {

                random_effects
            },

        map =
            if (
                length(
                    map
                ) == 0
            ) {

                NULL

            } else {

                map
            },

        DLL =
            DLL,

        silent =
            !verbose
    )


    # ============================================================
    # 25. Check initial objective
    # ============================================================

    initial_nll <-
        obj$fn(
            obj$par
        )


    if (
        !is.finite(
            initial_nll
        )
    ) {

        stop(
            "Initial joint negative log-likelihood is not finite."
        )
    }


    if (verbose) {

        message(
            "Initial negative log-likelihood: ",
            round(
                initial_nll,
                3
            )
        )

        message(
            "Optimizing joint likelihood..."
        )
    }


    # ============================================================
    # 26. Optimize
    # ============================================================

    opt <- stats::nlminb(

        start =
            obj$par,

        objective =
            obj$fn,

        gradient =
            obj$gr,

        control = list(

            iter.max =
                2000,

            eval.max =
                4000,

            rel.tol =
                1e-10
        )
    )


    # ============================================================
    # 27. Gradient diagnostic
    # ============================================================

    final_gradient <-
        tryCatch(

            obj$gr(
                opt$par
            ),

            error = function(e)
                rep(
                    NA_real_,
                    length(
                        opt$par
                    )
                )
        )


    max_gradient <-
        if (
            all(
                is.na(
                    final_gradient
                )
            )
        ) {

            NA_real_

        } else {

            max(
                abs(
                    final_gradient
                ),
                na.rm = TRUE
            )
        }


    if (verbose) {

        message(
            "Optimizer convergence code: ",
            opt$convergence
        )

        message(
            "Final negative log-likelihood: ",
            round(
                opt$objective,
                3
            )
        )

        message(
            "Optimizer message: ",
            opt$message
        )

        message(
            "Maximum absolute gradient: ",
            signif(
                max_gradient,
                5
            )
        )
    }


    # ============================================================
    # 28. Standard errors from marginal likelihood
    # ============================================================

    invisible(
        obj$fn(
            opt$par
        )
    )


    sdr <- tryCatch(

        TMB::sdreport(
            obj,
            par.fixed =
                opt$par,
            getJointPrecision =
                length(
                    random_effects
                ) > 0
        ),

        error = function(e) {

            warning(
                "sdreport failed: ",
                e$message
            )

            NULL
        }
    )


    # ============================================================
    # 29. Model reports
    # ============================================================

    report <- tryCatch(

        obj$report(),

        error = function(e) {

            warning(
                "TMB report failed: ",
                e$message
            )

            NULL
        }
    )


    # ============================================================
    # 30. Fixed-effect / parameter summary
    # ============================================================

    if (!is.null(sdr)) {

        fixed_matrix <-
            summary(
                sdr,
                "fixed"
            )


        fixed_summary <-
            as.data.frame(
                fixed_matrix
            )


        fixed_summary$parameter <-
            rownames(
                fixed_summary
            )


        rownames(
            fixed_summary
        ) <- NULL

    } else {

        fixed_summary <-
            data.frame()
    }


    # ============================================================
    # 31. ADREPORT summary
    # ============================================================

    if (!is.null(sdr)) {

        derived_matrix <-
            tryCatch(

                summary(
                    sdr,
                    "report"
                ),

                error = function(e)
                    NULL
            )


        if (!is.null(derived_matrix)) {

            derived_summary <-
                as.data.frame(
                    derived_matrix
                )


            derived_summary$parameter <-
                rownames(
                    derived_summary
                )


            rownames(
                derived_summary
            ) <- NULL

        } else {

            derived_summary <-
                data.frame()
        }

    } else {

        derived_summary <-
            data.frame()
    }


    # ============================================================
    # 32. Hessian / standard-error diagnostics
    # ============================================================

    finite_se <-
        if (
            nrow(
                fixed_summary
            ) > 0 &&
                "Std. Error" %in%
                names(
                    fixed_summary
                )
        ) {

            all(
                is.finite(
                    fixed_summary[[
                        "Std. Error"
                    ]]
                )
            )

        } else {

            NA
        }


    positive_definite_hessian <-
        if (!is.null(sdr)) {

            isTRUE(
                sdr$pdHess
            )

        } else {

            NA
        }


    # ============================================================
    # 33. Zero-inflation description
    # ============================================================

    zero_inflation_structure <-
        if (
            abundance_family %in%
            c(
                "zip",
                "zinb"
            )
        ) {

            if (random_zi_otu) {

                paste(
                    "read-level structural-zero probability with",
                    "a shared population-level intercept and",
                    "OTU-specific Gaussian random intercepts"
                )

            } else {

                "single shared read-level structural-zero probability"
            }

        } else {

            NULL
        }


    # ============================================================
    # 34. Return
    # ============================================================

    list(

        fit =
            opt,

        tmb_object =
            obj,

        sdreport =
            sdr,

        report =
            report,

        fixed_effects =
            fixed_summary,

        derived =
            derived_summary,

        site_data =
            site_df,

        sample_data =
            sample_df,

        long_df =
            dat,

        otu_stats =
            otu_stats,

        retained_otus =
            retained_otus,

        replication_summary =
            replication_summary,

        data_tmb =
            data_tmb,

        parameters_start =
            parameters,

        parameter_map =
            map,

        random_effects =
            random_effects,

        random_zi_otu =
            random_zi_otu,

        formulas = list(

            occupancy =
                occupancy_formula,

            capture =
                capture_formula,

            abundance =
                abundance_formula
        ),

        abundance_family =
            abundance_family,

        dispersion_structure =
            if (
                abundance_family %in%
                c(
                    "nbinom",
                    "zinb"
                )
            ) {

                "single shared NB2 dispersion parameter across OTUs"

            } else {

                NULL
            },

        zero_inflation_structure =
            zero_inflation_structure,

        convergence = list(

            code =
                opt$convergence,

            message =
                opt$message,

            objective =
                opt$objective,

            max_abs_gradient =
                max_gradient,

            finite_standard_errors =
                finite_se,

            positive_definite_hessian =
                positive_definite_hessian
        ),

        note = paste(
            "Joint observed-data likelihood.",
            "Occupancy is defined at the site x OTU level,",
            "capture at the biological sample x OTU level,",
            "and abundance at the PCR/read level.",
            "Latent occupancy Z and biological-sample capture A",
            "are analytically marginalized.",
            "Enabled Gaussian random effects are independent",
            "and integrated using the TMB Laplace approximation.",
            if (
                abundance_family %in%
                c(
                    "zip",
                    "zinb"
                )
            ) {

                paste(
                    "ZIP/ZINB requires at least two PCR replicates",
                    "per biological sample to distinguish read-level",
                    "zero inflation from sample-level capture failure.",
                    if (random_zi_otu) {
                        paste(
                            "The zero-inflation probability contains",
                            "an OTU-specific Gaussian random intercept."
                        )
                    } else {
                        paste(
                            "The zero-inflation probability is shared",
                            "across OTUs."
                        )
                    }
                )

            } else {

                ""
            }
        )
    )
}

There is one important detail in this version that differs from your uploaded function: I calculate n_otu, n_sample, and n_sample_otu before constructing the parameter list, and .otu is recreated after OTU filtering. This guarantees that b_occ_otu, b_cap_otu, b_abund_otu, and the new b_zi_otu all use the same filtered OTU indexing.

With the updated C++, you can then fit:

fit_zinb <- FitModel_joint(
    phyloseq = ps,

    site_col = "Site.name",
    sample_col = "Name",
    replicate_col = "Replicate",

    otu_col = "OTU",
    count_col = "y",

    occupancy_formula = ~ 1,
    capture_formula = ~ 1,
    abundance_formula = ~ 1,

    abundance_family = "zinb",

    random_occ_otu = TRUE,
    random_capture_otu = TRUE,
    random_abund_otu = TRUE,

    random_sample = TRUE,
    random_sample_otu = FALSE,

    # NEW
    random_zi_otu = TRUE,

    verbose = TRUE
)

The key new model is therefore

$$ \operatorname{logit}(\pi_m) = \beta_{\mathrm{ZI}} + b_m^{(\mathrm{ZI})}, \qquad b_m^{(\mathrm{ZI})} \sim N(0,\sigma_{\mathrm{ZI}}^2), $$

rather than forcing every OTU to share exactly the same structural-zero probability.

Pasted text(20260903-154558).txt
Document
is this C++ correct 

Yes — the overall C++ structure is correct and it matches the updated FitModel_joint() design, including the new OTU-specific zero-inflation random effect. I checked the hierarchy, random effects, ZIP/ZINB likelihood, and both levels of marginalization. The code correctly defines occupancy at site × OTU, capture at sample × OTU, and counts at PCR/read level.

There are, however, two changes I strongly recommend before treating it as the final version.

1. The new random_zi_otu implementation is correct

You correctly added

DATA_INTEGER(use_zi_otu);

PARAMETER_VECTOR(b_zi_otu);
PARAMETER(log_sd_zi_otu);

and

Type sd_zi_otu =
    exp(log_sd_zi_otu);

The random effects are then assigned

$$ b_m^{ZI}\sim N(0,\sigma_{ZI}^2), $$

which is implemented correctly here.

You also correctly construct

$$ \operatorname{logit}(\pi_m) = \beta_{ZI}+b_m^{ZI}, $$

through

Type eta_zi =
    zi_intercept;

if (
    use_zi_otu == 1 &&
    (
        family_code == 2 ||
        family_code == 3
    )
)
{
    eta_zi +=
        b_zi_otu(m);
}

and transform it to the probability scale.

So this part does exactly what we intended.

2. Change the ZIP/ZINB probability calculations to log-scale calculations

This is the main change I recommend.

At present you calculate

pi_zi_otu(m) =
    Type(1) /
    (
        Type(1) +
        exp(-eta_zi)
    );

and later use

log(pi_zi)

and

log(Type(1) - pi_zi)

in the ZIP/ZINB likelihood.

Mathematically this is correct, but numerically it is less stable than necessary. You already wrote stable functions at the top:

log_invlogit()
log1m_invlogit()

So use them directly.

Inside the observation loop, instead of only retrieving pi_zi, calculate the ZI linear predictor:

Type eta_zi =
    zi_intercept;

if (
    use_zi_otu == 1 &&
    (
        family_code == 2 ||
        family_code == 3
    )
)
{
    eta_zi +=
        b_zi_otu(
            row_otu(r)
        );
}

Type log_pi =
    log_invlogit(
        eta_zi
    );

Type log_1m_pi =
    log1m_invlogit(
        eta_zi
    );

Then ZIP becomes:

else if (
    family_code == 2
)
{
    Type log_count_component =
        dpois(
            y(r),
            mu,
            true
        );

    if (
        y(r) == Type(0)
    )
    {
        log_count(r) =
            logspace_add(
                log_pi,
                log_1m_pi +
                    log_count_component
            );
    }
    else
    {
        log_count(r) =
            log_1m_pi +
            log_count_component;
    }
}

And ZINB becomes:

else if (
    family_code == 3
)
{
    Type variance =
        mu +
        mu * mu / theta;

    Type log_count_component =
        dnbinom2(
            y(r),
            mu,
            variance,
            true
        );

    if (
        y(r) == Type(0)
    )
    {
        log_count(r) =
            logspace_add(
                log_pi,
                log_1m_pi +
                    log_count_component
            );
    }
    else
    {
        log_count(r) =
            log_1m_pi +
            log_count_component;
    }
}

This represents exactly the same model but avoids numerical problems when \(\pi_m\) becomes extremely close to 0 or 1.

3. Your capture marginalization is correct

For a positive sample, you correctly use

$$ \log L_{ijm} = \log p_{ijm} + \sum_k\log f(y_{ijkm}). $$

For an all-zero sample, you correctly calculate

$$ L_{ijm} = (1-p_{ijm}) + p_{ijm} \prod_k f(0). $$

Your implementation uses logspace_add() for this sum, which is exactly what we want.

This is important because capture is being marginalized at the biological-sample level, not at the individual read level.

4. Your occupancy marginalization is also correct

The code first combines all sample likelihoods belonging to the same site × OTU:

$$ \log L(\mathbf y_{im}\mid Z_{im}=1) = \sum_j \log L_{ijm}. $$

Then, for an all-zero site × OTU history, it calculates

$$ L_{im} = (1-\psi_{im}) + \psi_{im} L(\mathbf y_{im}\mid Z_{im}=1). $$

That is the correct site-level occupancy marginalization.

This also answers the earlier concern from your collaborator: this C++ code does not place occupancy at the PCR/read level. There is one \(Z_{im}\) per site × OTU, and multiple samples and PCR replicates are nested beneath it.

5. One conceptual point about pi_zi_otu

The code reports

REPORT(pi_zi_otu);
ADREPORT(pi_zi_otu);

which is useful.

But be careful how you describe this result. It is

$$ \pi_m = P(\text{structural zero at read/PCR level}\mid A_{ijm}=1), $$

not simply the overall probability of observing a zero.

A zero can now arise from four different mechanisms in the ZINB model:

$$ \boxed{ \text{Observed zero} \leftarrow \begin{cases} Z_{im}=0 & \text{species absent},\\ A_{ijm}=0 & \text{capture failure},\\ \text{ZI state}=1 & \text{PCR/read structural zero},\\ Y_{ijkm}=0 & \text{NB sampling zero}. \end{cases} } $$

That is precisely why the PCR replication requirement in FitModel_joint() is important.

6. One additional guard I would add

At the beginning of the C++ objective, add:

if (
    family_code < 0 ||
    family_code > 3
)
{
    error(
        "family_code must be 0, 1, 2, or 3."
    );
}

Otherwise an invalid family code could leave

lp = Type(0);

and silently produce an incorrect likelihood rather than failing.

Bottom line

The statistical architecture of your C++ is correct:

$$ \boxed{ Z_{im} \rightarrow A_{ijm} \rightarrow Y_{ijkm} } $$

with

$$ \operatorname{logit}(\psi_{im}) = X^{\psi}_{im}\beta_\psi+b_m^\psi, $$ $$ \operatorname{logit}(p_{ijm}) = X^p_{ijm}\beta_p+b_m^p, $$ $$ \log(\lambda_{ijkm}) = X^\lambda_{ijkm}\beta_\lambda +b_m^\lambda +b_j +b_{jm} +\mathrm{offset}_{ijkm}, $$

and now

$$ \boxed{ \operatorname{logit}(\pi_m) = \beta_{ZI} + b_m^{ZI}, \qquad b_m^{ZI}\sim N(0,\sigma_{ZI}^2). } $$

I would therefore keep the model structure, but change the ZIP/ZINB computations to use log_invlogit(eta_zi) and log1m_invlogit(eta_zi) directly. That will make the implementation considerably safer numerically without changing the statistical model.

Exactly. If you just want to test the updated FitModel_joint() on your real ps object, you do not need the long simulation code.

Because we added random_zi_otu, I would run:

fit_zinb <- FitModel_joint(
    phyloseq = ps,

    site_col = "Sampling.area.Name",
    sample_col = "Name",
    replicate_col = "Replicate",

    otu_col = "OTU",
    count_col = "y",

    occupancy_formula = ~ 1,
    capture_formula = ~ 1,
    abundance_formula = ~ 1,

    abundance_family = "zinb",

    random_occ_otu = TRUE,
    random_capture_otu = TRUE,
    random_abund_otu = TRUE,

    random_sample = TRUE,
    random_sample_otu = FALSE,

    # NEW: OTU-specific zero-inflation random effect
    random_zi_otu = TRUE,

    verbose = TRUE
)

This tests the full ZINB model with

$$ \operatorname{logit}(\psi_{im}) = \beta_{\mathrm{occ}}+b^{\mathrm{occ}}_m, $$ $$ \operatorname{logit}(p_{ijm}) = \beta_{\mathrm{cap}}+b^{\mathrm{cap}}_m, $$ $$ \log(\lambda_{ijkm}) = \beta_{\mathrm{abund}} +b^{\mathrm{abund}}_m +b^{\mathrm{sample}}_j, $$

and the newly added

$$ \operatorname{logit}(\pi_m) = \beta_{\mathrm{ZI}}+b^{\mathrm{ZI}}_m. $$

After it finishes, check:

fit_zinb$convergence

then:

fit_zinb$fixed_effects

and specifically the new ZI variance:

fit_zinb$fixed_effects[
    fit_zinb$fixed_effects$parameter == "log_sd_zi_otu",
]

Convert it to the SD scale:

exp(
    fit_zinb$fixed_effects$Estimate[
        fit_zinb$fixed_effects$parameter == "log_sd_zi_otu"
    ]
)

And most importantly, look at the estimated OTU-specific zero-inflation probabilities:

fit_zinb$report$pi_zi_otu

You can put OTU names beside them:

data.frame(
    OTU = levels(fit_zinb$long_df$.otu),
    pi_ZI = as.numeric(fit_zinb$report$pi_zi_otu)
)

If the modification is working, these pi_ZI values should now differ among OTUs. In the old model, every OTU would have had the same \(\pi\); that is the key thing we are testing here.

79 |   Index nsuper_et = 0; // Number of relaxed snodes in the original etree 
         |         ^
   eDNAModel.cpp:1079:5: error: expected ')'
    1079 |     {
         |     ^
   eDNAModel.cpp:1073:8: note: to match this '('
    1073 |     if (
         |        ^
   eDNAModel.cpp:1225:13: error: use of undeclared identifier 'log_pi'
    1225 |             log_pi,
         |             ^~~~~~
   eDNAModel.cpp:1302:13: error: use of undeclared identifier 'log_pi'
    1302 |             log_pi,
         |             ^~~~~~
   3 warnings and 3 errors generated.
   make: *** [eDNAModel.o] Error 1
   ERROR: compilation failed for package ‘eDNAModel’
─  removing ‘/private/var/folders/83/0h29m01s4rg2xr7s2zqp2pyh0000gn/T/RtmpOCi9Lp/devtools_install_d0cd781e4554/eDNAModel’
Error in (function (command = NULL, args = character(), error_on_status = TRUE, …:
! System command 'R' failed
---
Exit status: 1
stdout & stderr: <printed>
---

Yes. This is a C++ syntax error in the code I gave you, not a problem with TMB or your data.

The important error is the first one:

eDNAModel.cpp:1079:5: error: expected ')'
1079 |     {
     |     ^
eDNAModel.cpp:1073:8: note: to match this '('
1073 |     if (

Because that if statement is not closed correctly, the compiler subsequently cannot interpret the scope correctly, which is why you also get:

use of undeclared identifier 'log_pi'
Fix the block around line 1073

Find the section where log_pi and log_1m_pi are defined. Replace that entire portion with this:

// ==========================================================
// Zero-inflation linear predictor for this observation
// ==========================================================

Type eta_zi =
  zi_intercept;


// Add OTU-specific zero-inflation random effect
// only for ZIP or ZINB models.

if (
  use_zi_otu == 1 &&
  (
    family_code == 2 ||
    family_code == 3
  )
)
{
  eta_zi +=
    b_zi_otu(
      row_otu(r)
    );
}


// ==========================================================
// Numerically stable zero-inflation probabilities
// ==========================================================
//
// log_pi:
//
//     log(pi_m)
//
// log_1m_pi:
//
//     log(1 - pi_m)
//
// These are calculated directly from the linear predictor
// rather than first calculating pi_m and then taking logs.
// ==========================================================

Type log_pi =
  log_invlogit(
    eta_zi
  );


Type log_1m_pi =
  log1m_invlogit(
    eta_zi
  );

The important part is that the if must be exactly:

if (
  use_zi_otu == 1 &&
  (
    family_code == 2 ||
    family_code == 3
  )
)
{

There are two closing parentheses before { because one closes

(
    family_code == 2 ||
    family_code == 3
)

and the second closes the complete

if (
    ...
)
Then log_pi will be available to ZIP

Your ZIP block should be:

else if (
  family_code == 2
)
{
  Type log_count_component =
    dpois(
      y(r),
      mu,
      true
    );

  if (
    y(r) == Type(0)
  )
  {
    lp =
      logspace_add(
        log_pi,
        log_1m_pi +
          log_count_component
      );
  }
  else
  {
    lp =
      log_1m_pi +
      log_count_component;
  }
}

And ZINB:

else if (
  family_code == 3
)
{
  Type variance =
    mu +
    mu * mu / theta;

  Type log_count_component =
    dnbinom2(
      y(r),
      mu,
      variance,
      true
    );

  if (
    y(r) == Type(0)
  )
  {
    lp =
      logspace_add(
        log_pi,
        log_1m_pi +
          log_count_component
      );
  }
  else
  {
    lp =
      log_1m_pi +
      log_count_component;
  }
}

The two later errors at lines 1225 and 1302:

use of undeclared identifier 'log_pi'

are most likely secondary errors caused by the malformed if statement. Once the parentheses are corrected, log_pi is declared inside the observation loop before both ZIP and ZINB branches and is therefore in scope.

The warnings such as:

Index nsuper_et = 0

come from Eigen/TMB headers and are not what caused compilation to fail.

After correcting the block, run from your package directory:

devtools::document()

then:

devtools::install()

Restart R, then:

library(eDNAModel)

"eDNAModel" %in% names(getLoadedDLLs())

You want:

[1] TRUE

Then test:

fit_zinb <- FitModel_joint(
    phyloseq = ps,

    site_col = "Sampling.area.Name",
    sample_col = "Name",
    replicate_col = "Replicate",

    otu_col = "OTU",
    count_col = "y",

    occupancy_formula = ~ 1,
    capture_formula = ~ 1,
    abundance_formula = ~ 1,

    abundance_family = "zinb",

    random_occ_otu = TRUE,
    random_capture_otu = TRUE,
    random_abund_otu = TRUE,
    random_sample = TRUE,
    random_sample_otu = FALSE,

    random_zi_otu = TRUE,

    verbose = TRUE
)

If you paste lines roughly 1050–1090 of your actual src/eDNAModel.cpp, I can point to the exact missing parenthesis rather than having you replace a larger section.

write the full corrected code

Yes. Below is a clean, corrected full src/eDNAModel.cpp. I have kept the hierarchy and features we discussed, but fixed the parenthesis/scope problem that caused log_pi to be undeclared.

#include <TMB.hpp>


// ============================================================
// Numerically stable logistic helper functions
// ============================================================

// log(invlogit(x))
template<class Type>
Type log_invlogit(Type x)
{
  return -logspace_add(Type(0), -x);
}


// log(1 - invlogit(x))
template<class Type>
Type log1m_invlogit(Type x)
{
  return -logspace_add(Type(0), x);
}


// ============================================================
// Joint hierarchical eDNA model
// ============================================================
//
// Hierarchy:
//
// SITE x OTU:
//
//   Z_im ~ Bernoulli(psi_im)
//
// SAMPLE x OTU:
//
//   A_ijm | Z_im = 1 ~ Bernoulli(p_ijm)
//
// PCR replicate:
//
//   Y_ijkm | A_ijm = 1 ~ Count(lambda_ijkm)
//
// where:
//
//   i = site
//   j = biological sample
//   k = PCR replicate
//   m = OTU
//
// Occupancy:
//
//   logit(psi_im)
//     = X_occ beta_occ
//       + b_occ_m
//
// Capture:
//
//   logit(p_ijm)
//     = X_cap beta_cap
//       + b_cap_m
//
// Abundance:
//
//   log(lambda_ijkm)
//     = X_abund beta_abund
//       + offset
//       + b_abund_m
//       + b_sample_j
//       + b_sampleOTU_jm
//
// Zero inflation:
//
//   logit(pi_m)
//     = zi_intercept
//       + b_zi_m
//
// The latent capture state A is marginalized analytically.
//
// The latent occupancy state Z is marginalized analytically.
//
// Gaussian random effects are integrated by TMB using the
// Laplace approximation when supplied to MakeADFun(random=...).
// ============================================================


template<class Type>
Type objective_function<Type>::operator() ()
{

  // ==========================================================
  // 1. DATA
  // ==========================================================

  DATA_VECTOR(y);

  // y:
  // PCR/read-level observed counts.


  // ----------------------------------------------------------
  // Design matrices
  // ----------------------------------------------------------

  DATA_MATRIX(X_occ);

  // One row per site x OTU.


  DATA_MATRIX(X_cap);

  // One row per sample x OTU.


  DATA_MATRIX(X_abund);

  // One row per PCR/read observation.


  // ----------------------------------------------------------
  // Abundance offset
  // ----------------------------------------------------------

  DATA_VECTOR(offset_abund);


  // ----------------------------------------------------------
  // Hierarchical mapping indices
  // ----------------------------------------------------------

  DATA_IVECTOR(occ_otu);

  // site x OTU group -> OTU


  DATA_IVECTOR(sample_otu);

  // sample x OTU group -> OTU


  DATA_IVECTOR(row_sample_group);

  // PCR/read row -> sample x OTU group


  DATA_IVECTOR(sample_site_group);

  // sample x OTU group -> site x OTU group


  DATA_IVECTOR(row_otu);

  // PCR/read row -> OTU


  DATA_IVECTOR(row_sample_id);

  // PCR/read row -> biological sample


  DATA_IVECTOR(row_sample_otu_re);

  // PCR/read row -> sample x OTU random-effect level


  // ----------------------------------------------------------
  // Positive-history indicators
  // ----------------------------------------------------------

  DATA_IVECTOR(sample_positive);

  // sample_positive(s) = 1 if at least one PCR replicate
  // within sample x OTU group s is positive.


  DATA_IVECTOR(site_positive);

  // site_positive(g) = 1 if at least one observation
  // within site x OTU group g is positive.


  // ----------------------------------------------------------
  // Number of groups
  // ----------------------------------------------------------

  DATA_INTEGER(n_site_groups);

  DATA_INTEGER(n_sample_groups);


  // ----------------------------------------------------------
  // Distribution
  // ----------------------------------------------------------

  DATA_INTEGER(family_code);

  // 0 = Poisson
  // 1 = NB2
  // 2 = ZIP
  // 3 = ZINB


  // ==========================================================
  // 2. RANDOM-EFFECT SWITCHES
  // ==========================================================

  DATA_INTEGER(use_occ_otu);

  DATA_INTEGER(use_cap_otu);

  DATA_INTEGER(use_abund_otu);

  DATA_INTEGER(use_sample_re);

  DATA_INTEGER(use_sample_otu_re);

  DATA_INTEGER(use_zi_otu);


  // ==========================================================
  // 3. VALIDATE FAMILY
  // ==========================================================

  if (
    family_code < 0 ||
    family_code > 3
  )
  {
    error(
      "family_code must be 0, 1, 2, or 3."
    );
  }


  // ==========================================================
  // 4. FIXED EFFECTS
  // ==========================================================

  PARAMETER_VECTOR(beta_occ);

  PARAMETER_VECTOR(beta_cap);

  PARAMETER_VECTOR(beta_abund);


  // ==========================================================
  // 5. RANDOM EFFECTS
  // ==========================================================

  PARAMETER_VECTOR(b_occ_otu);

  // OTU-specific occupancy random effects.


  PARAMETER_VECTOR(b_cap_otu);

  // OTU-specific capture random effects.


  PARAMETER_VECTOR(b_abund_otu);

  // OTU-specific abundance random effects.


  PARAMETER_VECTOR(b_sample);

  // Biological-sample abundance random effects.


  PARAMETER_VECTOR(b_sample_otu);

  // Sample x OTU abundance random effects.


  PARAMETER_VECTOR(b_zi_otu);

  // OTU-specific zero-inflation random effects.


  // ==========================================================
  // 6. RANDOM-EFFECT STANDARD DEVIATIONS
  // ==========================================================

  PARAMETER(log_sd_occ_otu);

  PARAMETER(log_sd_cap_otu);

  PARAMETER(log_sd_abund_otu);

  PARAMETER(log_sd_sample);

  PARAMETER(log_sd_sample_otu);

  PARAMETER(log_sd_zi_otu);


  Type sd_occ_otu =
    exp(log_sd_occ_otu);


  Type sd_cap_otu =
    exp(log_sd_cap_otu);


  Type sd_abund_otu =
    exp(log_sd_abund_otu);


  Type sd_sample =
    exp(log_sd_sample);


  Type sd_sample_otu =
    exp(log_sd_sample_otu);


  Type sd_zi_otu =
    exp(log_sd_zi_otu);


  // ==========================================================
  // 7. NB DISPERSION
  // ==========================================================

  PARAMETER(log_theta);


  Type theta =
    exp(log_theta);

  // NB2 parameterization:
  //
  // E(Y) = mu
  //
  // Var(Y) = mu + mu^2 / theta
  //
  // Currently theta is shared across OTUs.


  // ==========================================================
  // 8. ZERO-INFLATION INTERCEPT
  // ==========================================================

  PARAMETER(zi_intercept);

  // For ZIP/ZINB:
  //
  // logit(pi_m)
  //   =
  // zi_intercept + b_zi_otu(m)


  // ==========================================================
  // 9. NEGATIVE LOG-LIKELIHOOD
  // ==========================================================

  Type nll =
    Type(0);


  // ==========================================================
  // 10. RANDOM-EFFECT DISTRIBUTIONS
  // ==========================================================
  //
  // All random-effect blocks are independent.
  // ==========================================================


  // ----------------------------------------------------------
  // Occupancy OTU random effects
  // ----------------------------------------------------------

  if (
    use_occ_otu == 1
  )
  {
    for (
      int m = 0;
      m < b_occ_otu.size();
      m++
    )
    {
      nll -=
        dnorm(
          b_occ_otu(m),
          Type(0),
          sd_occ_otu,
          true
        );
    }
  }


  // ----------------------------------------------------------
  // Capture OTU random effects
  // ----------------------------------------------------------

  if (
    use_cap_otu == 1
  )
  {
    for (
      int m = 0;
      m < b_cap_otu.size();
      m++
    )
    {
      nll -=
        dnorm(
          b_cap_otu(m),
          Type(0),
          sd_cap_otu,
          true
        );
    }
  }


  // ----------------------------------------------------------
  // Abundance OTU random effects
  // ----------------------------------------------------------

  if (
    use_abund_otu == 1
  )
  {
    for (
      int m = 0;
      m < b_abund_otu.size();
      m++
    )
    {
      nll -=
        dnorm(
          b_abund_otu(m),
          Type(0),
          sd_abund_otu,
          true
        );
    }
  }


  // ----------------------------------------------------------
  // Sample random effects
  // ----------------------------------------------------------

  if (
    use_sample_re == 1
  )
  {
    for (
      int j = 0;
      j < b_sample.size();
      j++
    )
    {
      nll -=
        dnorm(
          b_sample(j),
          Type(0),
          sd_sample,
          true
        );
    }
  }


  // ----------------------------------------------------------
  // Sample x OTU random effects
  // ----------------------------------------------------------

  if (
    use_sample_otu_re == 1
  )
  {
    for (
      int q = 0;
      q < b_sample_otu.size();
      q++
    )
    {
      nll -=
        dnorm(
          b_sample_otu(q),
          Type(0),
          sd_sample_otu,
          true
        );
    }
  }


  // ----------------------------------------------------------
  // OTU-specific zero-inflation random effects
  // ----------------------------------------------------------
  //
  // Only active for ZIP/ZINB.
  // ----------------------------------------------------------

  if (
    use_zi_otu == 1 &&
    (
      family_code == 2 ||
      family_code == 3
    )
  )
  {
    for (
      int m = 0;
      m < b_zi_otu.size();
      m++
    )
    {
      nll -=
        dnorm(
          b_zi_otu(m),
          Type(0),
          sd_zi_otu,
          true
        );
    }
  }


  // ==========================================================
  // 11. OCCUPANCY LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_occ =
    X_occ * beta_occ;


  if (
    use_occ_otu == 1
  )
  {
    for (
      int g = 0;
      g < n_site_groups;
      g++
    )
    {
      eta_occ(g) +=
        b_occ_otu(
          occ_otu(g)
        );
    }
  }


  // ==========================================================
  // 12. CAPTURE LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_cap =
    X_cap * beta_cap;


  if (
    use_cap_otu == 1
  )
  {
    for (
      int s = 0;
      s < n_sample_groups;
      s++
    )
    {
      eta_cap(s) +=
        b_cap_otu(
          sample_otu(s)
        );
    }
  }


  // ==========================================================
  // 13. ABUNDANCE LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_abund =
    X_abund * beta_abund;


  for (
    int r = 0;
    r < y.size();
    r++
  )
  {

    // Add abundance offset.

    eta_abund(r) +=
      offset_abund(r);


    // OTU abundance effect.

    if (
      use_abund_otu == 1
    )
    {
      eta_abund(r) +=
        b_abund_otu(
          row_otu(r)
        );
    }


    // Biological-sample abundance effect.

    if (
      use_sample_re == 1
    )
    {
      eta_abund(r) +=
        b_sample(
          row_sample_id(r)
        );
    }


    // Sample x OTU abundance effect.

    if (
      use_sample_otu_re == 1
    )
    {
      eta_abund(r) +=
        b_sample_otu(
          row_sample_otu_re(r)
        );
    }
  }


  // ==========================================================
  // 14. ZERO-INFLATION PREDICTOR BY OTU
  // ==========================================================

  vector<Type> eta_zi_otu(
    b_zi_otu.size()
  );


  vector<Type> pi_zi_otu(
    b_zi_otu.size()
  );


  for (
    int m = 0;
    m < b_zi_otu.size();
    m++
  )
  {

    eta_zi_otu(m) =
      zi_intercept;


    if (
      use_zi_otu == 1 &&
      (
        family_code == 2 ||
        family_code == 3
      )
    )
    {
      eta_zi_otu(m) +=
        b_zi_otu(m);
    }


    pi_zi_otu(m) =
      exp(
        log_invlogit(
          eta_zi_otu(m)
        )
      );
  }


  // ==========================================================
  // 15. PCR/READ-LEVEL COUNT LIKELIHOOD
  // ==========================================================

  vector<Type> log_count(
    y.size()
  );


  for (
    int r = 0;
    r < y.size();
    r++
  )
  {

    // --------------------------------------------------------
    // Mean abundance
    // --------------------------------------------------------

    Type mu =
      exp(
        eta_abund(r)
      );


    Type lp =
      Type(0);


    // --------------------------------------------------------
    // Zero-inflation predictor
    // --------------------------------------------------------

    Type eta_zi =
      zi_intercept;


    if (
      use_zi_otu == 1 &&
      (
        family_code == 2 ||
        family_code == 3
      )
    )
    {
      eta_zi +=
        b_zi_otu(
          row_otu(r)
        );
    }


    // --------------------------------------------------------
    // Stable log zero-inflation probabilities
    // --------------------------------------------------------

    Type log_pi =
      log_invlogit(
        eta_zi
      );


    Type log_1m_pi =
      log1m_invlogit(
        eta_zi
      );


    // ========================================================
    // Poisson
    // ========================================================

    if (
      family_code == 0
    )
    {
      lp =
        dpois(
          y(r),
          mu,
          true
        );
    }


    // ========================================================
    // Negative Binomial 2
    // ========================================================

    else if (
      family_code == 1
    )
    {
      Type variance =
        mu +
        mu * mu / theta;


      lp =
        dnbinom2(
          y(r),
          mu,
          variance,
          true
        );
    }


    // ========================================================
    // Zero-Inflated Poisson
    // ========================================================
    //
    // If y = 0:
    //
    // P(Y=0)
    //   =
    // pi
    // +
    // (1-pi) Pois(0 | mu)
    //
    // If y > 0:
    //
    // P(Y=y)
    //   =
    // (1-pi) Pois(y | mu)
    // ========================================================

    else if (
      family_code == 2
    )
    {
      Type log_count_component =
        dpois(
          y(r),
          mu,
          true
        );


      if (
        y(r) == Type(0)
      )
      {
        lp =
          logspace_add(
            log_pi,
            log_1m_pi +
              log_count_component
          );
      }
      else
      {
        lp =
          log_1m_pi +
          log_count_component;
      }
    }


    // ========================================================
    // Zero-Inflated Negative Binomial 2
    // ========================================================
    //
    // If y = 0:
    //
    // P(Y=0)
    //   =
    // pi
    // +
    // (1-pi) NB(0 | mu, theta)
    //
    // If y > 0:
    //
    // P(Y=y)
    //   =
    // (1-pi) NB(y | mu, theta)
    // ========================================================

    else if (
      family_code == 3
    )
    {
      Type variance =
        mu +
        mu * mu / theta;


      Type log_count_component =
        dnbinom2(
          y(r),
          mu,
          variance,
          true
        );


      if (
        y(r) == Type(0)
      )
      {
        lp =
          logspace_add(
            log_pi,
            log_1m_pi +
              log_count_component
          );
      }
      else
      {
        lp =
          log_1m_pi +
          log_count_component;
      }
    }


    log_count(r) =
      lp;
  }


  // ==========================================================
  // 16. COMBINE PCR REPLICATES WITHIN SAMPLE x OTU
  // ==========================================================
  //
  // Conditional on A_ijm = 1:
  //
  // L_count_ijm
  //   =
  // product_k f(y_ijkm)
  //
  // Therefore:
  //
  // log L_count_ijm
  //   =
  // sum_k log f(y_ijkm)
  // ==========================================================

  vector<Type> log_sample_count(
    n_sample_groups
  );


  log_sample_count.setZero();


  for (
    int r = 0;
    r < y.size();
    r++
  )
  {
    int s =
      row_sample_group(r);


    log_sample_count(s) +=
      log_count(r);
  }


  // ==========================================================
  // 17. MARGINALIZE SAMPLE-LEVEL CAPTURE A
  // ==========================================================
  //
  // A_ijm exists once per biological sample x OTU.
  //
  // If any PCR replicate is positive:
  //
  //   A_ijm = 1
  //
  // and:
  //
  //   L_ijm
  //     =
  //   p_ijm * product_k f(y_ijkm)
  //
  //
  // If all PCR replicates are zero:
  //
  //   L_ijm
  //     =
  //   (1-p_ijm)
  //     +
  //   p_ijm * product_k f(0)
  //
  // This sums over A_ijm = 0 and A_ijm = 1.
  // ==========================================================

  vector<Type> log_sample(
    n_sample_groups
  );


  for (
    int s = 0;
    s < n_sample_groups;
    s++
  )
  {

    Type log_p =
      log_invlogit(
        eta_cap(s)
      );


    Type log_1mp =
      log1m_invlogit(
        eta_cap(s)
      );


    // --------------------------------------------------------
    // At least one positive PCR replicate
    // --------------------------------------------------------

    if (
      sample_positive(s) == 1
    )
    {
      log_sample(s) =
        log_p +
        log_sample_count(s);
    }


    // --------------------------------------------------------
    // All PCR replicates are zero
    // --------------------------------------------------------

    else
    {
      log_sample(s) =
        logspace_add(
          log_1mp,
          log_p +
            log_sample_count(s)
        );
    }
  }


  // ==========================================================
  // 18. COMBINE SAMPLES WITHIN SITE x OTU
  // ==========================================================
  //
  // Conditional on occupancy Z_im = 1:
  //
  // L(data_im | Z_im = 1)
  //
  //   =
  //
  // product_j L_ijm
  //
  // therefore:
  //
  // log L(data_im | Z_im = 1)
  //
  //   =
  //
  // sum_j log L_ijm
  // ==========================================================

  vector<Type> log_site_conditional(
    n_site_groups
  );


  log_site_conditional.setZero();


  for (
    int s = 0;
    s < n_sample_groups;
    s++
  )
  {
    int g =
      sample_site_group(s);


    log_site_conditional(g) +=
      log_sample(s);
  }


  // ==========================================================
  // 19. MARGINALIZE SITE-LEVEL OCCUPANCY Z
  // ==========================================================
  //
  // Z_im exists ONCE per site x OTU.
  //
  // This is the important distinction:
  //
  // occupancy is NOT at the PCR/read observation level.
  //
  //
  // If site x OTU has any positive observation:
  //
  //   Z_im = 1
  //
  // therefore:
  //
  //   L_im
  //     =
  //   psi_im *
  //   L(data_im | Z_im=1)
  //
  //
  // If the entire site x OTU history is zero:
  //
  //   L_im
  //     =
  //   (1-psi_im)
  //     +
  //   psi_im *
  //   L(all-zero data | Z_im=1)
  //
  // This analytically sums over Z_im.
  // ==========================================================

  vector<Type> psi(
    n_site_groups
  );


  vector<Type> log_site(
    n_site_groups
  );


  for (
    int g = 0;
    g < n_site_groups;
    g++
  )
  {

    // Occupancy probability on natural scale.

    psi(g) =
      exp(
        log_invlogit(
          eta_occ(g)
        )
      );


    // Stable log occupancy probability.

    Type log_psi =
      log_invlogit(
        eta_occ(g)
      );


    // Stable log non-occupancy probability.

    Type log_1mpsi =
      log1m_invlogit(
        eta_occ(g)
      );


    // --------------------------------------------------------
    // Positive site x OTU history
    // --------------------------------------------------------

    if (
      site_positive(g) == 1
    )
    {
      log_site(g) =
        log_psi +
        log_site_conditional(g);
    }


    // --------------------------------------------------------
    // Entire site x OTU history is zero
    // --------------------------------------------------------

    else
    {
      log_site(g) =
        logspace_add(
          log_1mpsi,
          log_psi +
            log_site_conditional(g)
        );
    }


    // Add contribution to negative log-likelihood.

    nll -=
      log_site(g);
  }


  // ==========================================================
  // 20. DERIVED CAPTURE PROBABILITIES
  // ==========================================================

  vector<Type> capture_prob(
    n_sample_groups
  );


  for (
    int s = 0;
    s < n_sample_groups;
    s++
  )
  {
    capture_prob(s) =
      exp(
        log_invlogit(
          eta_cap(s)
        )
      );
  }


  // ==========================================================
  // 21. DERIVED ABUNDANCE
  // ==========================================================

  vector<Type> lambda(
    y.size()
  );


  for (
    int r = 0;
    r < y.size();
    r++
  )
  {
    lambda(r) =
      exp(
        eta_abund(r)
      );
  }


  // ==========================================================
  // 22. REPORT
  // ==========================================================

  REPORT(
    eta_occ
  );


  REPORT(
    eta_cap
  );


  REPORT(
    eta_abund
  );


  REPORT(
    eta_zi_otu
  );


  REPORT(
    psi
  );


  REPORT(
    capture_prob
  );


  REPORT(
    lambda
  );


  REPORT(
    pi_zi_otu
  );


  REPORT(
    log_sample_count
  );


  REPORT(
    log_sample
  );


  REPORT(
    log_site_conditional
  );


  REPORT(
    log_site
  );


  // ==========================================================
  // 23. ADREPORT
  // ==========================================================
  //
  // ADREPORT allows sdreport() to obtain approximate SEs
  // for these transformed parameters using the delta method.
  // ==========================================================

  ADREPORT(
    psi
  );


  ADREPORT(
    capture_prob
  );


  if (
    family_code == 2 ||
    family_code == 3
  )
  {
    ADREPORT(
      pi_zi_otu
    );
  }


  // ==========================================================
  // 24. RETURN NEGATIVE LOG-LIKELIHOOD
  // ==========================================================

  return nll;
}

This version fixes the exact compilation issue because log_pi and log_1m_pi are declared inside the PCR observation loop but outside the ZIP/ZINB branches, so both branches can access them:

Type log_pi =
  log_invlogit(
    eta_zi
  );

Type log_1m_pi =
  log1m_invlogit(
    eta_zi
  );

There is also an important consistency requirement: your R FitModel_joint() must now create all six random-effect vectors, including b_zi_otu, and pass use_zi_otu to the TMB data list. If random_zi_otu = FALSE, the R wrapper should map b_zi_otu and log_sd_zi_otu; for Poisson/NB it should also map the ZI parameters. Otherwise the C++ can compile perfectly but MakeADFun() can still fail or leave irrelevant parameters unidentified.

#include <TMB.hpp>

// ============================================================
// Numerically stable logistic helper functions
// ============================================================
//
// TMB already provides logspace_add(a,b):
//
// logspace_add(a,b)
// = log(exp(a) + exp(b))
//
// This is preferable to explicitly exponentiating probabilities
// because likelihood terms can become extremely small.
//
// We use these helper functions repeatedly for occupancy,
// capture, and zero-inflation probabilities.
// ============================================================

template<class Type>
Type log_invlogit(Type x)
{
// ----------------------------------------------------------
// Computes:
//
// log(invlogit(x))
//
// where:
//
// invlogit(x)
// =
// 1 / (1 + exp(-x))
//
// Therefore:
//
// log(invlogit(x))
// =
// -log(1 + exp(-x))
//
// Using logspace_add makes this numerically stable.
// ----------------------------------------------------------

return -logspace_add(
Type(0),
-x
);
}




template<class Type>
Type log1m_invlogit(Type x)
{
// ----------------------------------------------------------
// Computes:
//
// log(1 - invlogit(x))
//
// Because:
//
// 1 - invlogit(x)
// =
// 1 / (1 + exp(x))
//
// therefore:
//
// log(1-p)
// =
// -log(1 + exp(x)).
// ----------------------------------------------------------

return -logspace_add(
Type(0),
x
);
}




// ============================================================
// Joint observed-data eDNA likelihood
// ============================================================
//
// Hierarchy:
//
// ------------------------------------------------------------
// OCCUPANCY LEVEL
//
// Z_im ~ Bernoulli(psi_im)
//
// i = site
// m = OTU
//
// There is ONE occupancy state per site x OTU.
//
// ------------------------------------------------------------
// CAPTURE LEVEL
//
// A_ijm | Z_im = 1
// ~ Bernoulli(p_ijm)
//
// j = biological sample
//
// ------------------------------------------------------------
// PCR / READ LEVEL
//
// Y_ijkm | A_ijm = 1
// ~ Count(lambda_ijkm)
//
// k = PCR replicate
//
// The count distribution may be:
//
// Poisson
// NB2
// ZIP
// ZINB
//
// ------------------------------------------------------------
//
// IMPORTANT:
//
// Z and A are NOT sampled or optimized directly.
//
// Instead:
//
// A is analytically marginalized
//
// and then:
//
// Z is analytically marginalized.
//
// Gaussian random effects are integrated by TMB using the
// Laplace approximation when passed through the random
// argument of MakeADFun().
// ============================================================

template<class Type>
Type objective_function<Type>::operator() ()
{

// ==========================================================
// 1. DATA
// ==========================================================

DATA_VECTOR(y);

// Observed sequencing counts.
//
// Each y(r) is one PCR/read observation.

// ----------------------------------------------------------
// Design matrices
// ----------------------------------------------------------

DATA_MATRIX(X_occ);

// Occupancy design matrix.
//
// One row per:
//
// SITE x OTU
//
// NOT one row per PCR observation.

DATA_MATRIX(X_cap);

// Capture design matrix.
//
// One row per:
//
// BIOLOGICAL SAMPLE x OTU

DATA_MATRIX(X_abund);

// Abundance design matrix.
//
// One row per PCR/read observation.

// ----------------------------------------------------------
// Abundance offset
// ----------------------------------------------------------

DATA_VECTOR(offset_abund);

// For example:
//
// log(total_reads)
//
// representing sequencing depth.

// ----------------------------------------------------------
// Hierarchical mapping indices
// ----------------------------------------------------------

DATA_IVECTOR(occ_otu);

// For each SITE x OTU occupancy group g,
// gives the OTU index m.
//
// Used by:
//
// b_occ_otu(m)

DATA_IVECTOR(sample_otu);

// For every SAMPLE x OTU group s,
// gives the OTU index m.
//
// Used by:
//
// b_cap_otu(m)

DATA_IVECTOR(row_sample_group);

// Maps:
//
// PCR/read row r
//
// to:
//
// SAMPLE x OTU group s.

DATA_IVECTOR(sample_site_group);

// Maps:
//
// SAMPLE x OTU group s
//
// to:
//
// SITE x OTU group g.
//
// This ensures multiple samples from a site share
// the same occupancy state.

DATA_IVECTOR(row_otu);

// Maps each PCR/read observation r to its OTU m.
//
// Used for:
//
// abundance OTU random effects
//
// and now also:
//
// zero-inflation OTU random effects.

DATA_IVECTOR(row_sample_id);

// Maps every PCR/read observation to its biological sample.

DATA_IVECTOR(row_sample_otu_re);

// Maps every PCR/read observation to a unique
// sample x OTU random-effect level.

// ----------------------------------------------------------
// Positive-history indicators
// ----------------------------------------------------------

DATA_IVECTOR(sample_positive);

// sample_positive(s) = 1 if at least one PCR replicate
// for sample x OTU group s is positive.
//
// If positive:
//
// A_ijm MUST equal 1.

DATA_IVECTOR(site_positive);

// site_positive(g) = 1 if at least one observation anywhere
// inside site x OTU group g is positive.
//
// If positive:
//
// Z_im MUST equal 1.

// ----------------------------------------------------------
// Number of hierarchical groups
// ----------------------------------------------------------

DATA_INTEGER(n_site_groups);

// Number of SITE x OTU combinations.

DATA_INTEGER(n_sample_groups);

// Number of SAMPLE x OTU combinations.

// ----------------------------------------------------------
// Count family
// ----------------------------------------------------------
//
// 0 = Poisson
// 1 = NB2
// 2 = ZIP
// 3 = ZINB
// ----------------------------------------------------------

DATA_INTEGER(family_code);

// ==========================================================
// 2. RANDOM-EFFECT SWITCHES
// ==========================================================

DATA_INTEGER(use_occ_otu);

// OTU random effect for occupancy.

DATA_INTEGER(use_cap_otu);

// OTU random effect for capture.

DATA_INTEGER(use_abund_otu);

// OTU random effect for abundance.

DATA_INTEGER(use_sample_re);

// Biological-sample abundance random effect.

DATA_INTEGER(use_sample_otu_re);

// Sample x OTU abundance random effect.

DATA_INTEGER(use_zi_otu);

// NEW:
//
// OTU-specific random effect for ZIP/ZINB structural zeros.
//
// 1 = estimate taxon-specific PCR dropout heterogeneity.
// 0 = all OTUs share one zero-inflation probability.




// ==========================================================
// 3. FIXED-EFFECT PARAMETERS
// ==========================================================

PARAMETER_VECTOR(beta_occ);

// Occupancy fixed effects.

PARAMETER_VECTOR(beta_cap);

// Capture fixed effects.

PARAMETER_VECTOR(beta_abund);

// Abundance fixed effects.




// ==========================================================
// 4. RANDOM EFFECTS
// ==========================================================

PARAMETER_VECTOR(b_occ_otu);

// OTU-specific occupancy deviations.

PARAMETER_VECTOR(b_cap_otu);

// OTU-specific capture deviations.

PARAMETER_VECTOR(b_abund_otu);

// OTU-specific abundance deviations.

PARAMETER_VECTOR(b_sample);

// Biological-sample abundance random effects.

PARAMETER_VECTOR(b_sample_otu);

// Sample x OTU abundance random effects.

PARAMETER_VECTOR(b_zi_otu);

// NEW:
//
// OTU-specific deviations in structural-zero probability.
//
// Mathematically:
//
// b_zi_otu(m)
// ~ Normal(0, sigma_zi^2)
//
// and:
//
// logit(pi_m)
// =
// zi_intercept
// +
// b_zi_otu(m).




// ==========================================================
// 5. RANDOM-EFFECT STANDARD DEVIATIONS
// ==========================================================

PARAMETER(log_sd_occ_otu);

PARAMETER(log_sd_cap_otu);

PARAMETER(log_sd_abund_otu);

PARAMETER(log_sd_sample);

PARAMETER(log_sd_sample_otu);

PARAMETER(log_sd_zi_otu);

// NEW:
//
// log standard deviation of OTU-specific
// zero-inflation random effects.

Type sd_occ_otu =
exp(
log_sd_occ_otu
);

Type sd_cap_otu =
exp(
log_sd_cap_otu
);

Type sd_abund_otu =
exp(
log_sd_abund_otu
);

Type sd_sample =
exp(
log_sd_sample
);

Type sd_sample_otu =
exp(
log_sd_sample_otu
);

Type sd_zi_otu =
exp(
log_sd_zi_otu
);

// exp() guarantees positive standard deviations.




// ==========================================================
// 6. ABUNDANCE-DISTRIBUTION PARAMETERS
// ==========================================================

PARAMETER(log_theta);

// Shared NB2 dispersion parameter.

Type theta =
exp(
log_theta
);

// One theta is currently shared across all OTUs:
//
// Var(Y)
// =
// mu + mu^2 / theta.

PARAMETER(zi_intercept);

// Global zero-inflation intercept.
//
// This now represents the population-average PCR dropout
// tendency across OTUs.
//
// OTUs can deviate from this through b_zi_otu.




// ==========================================================
// 7. NEGATIVE LOG-LIKELIHOOD
// ==========================================================

Type nll =
Type(0);

// TMB minimizes nll:
//
// nll = -log L.




// ==========================================================
// 8. GAUSSIAN RANDOM-EFFECT DENSITIES
// ==========================================================

if (use_occ_otu == 1)
{

for (
  int m = 0;
  m < b_occ_otu.size();
  m++
)
{

  nll -= dnorm(
    b_occ_otu(m),
    Type(0),
    sd_occ_otu,
    true
  );

  // b_occ_otu(m)
  // ~ N(0, sigma_occ^2).
}

}




if (use_cap_otu == 1)
{

for (
  int m = 0;
  m < b_cap_otu.size();
  m++
)
{

  nll -= dnorm(
    b_cap_otu(m),
    Type(0),
    sd_cap_otu,
    true
  );
}

}




if (use_abund_otu == 1)
{

for (
  int m = 0;
  m < b_abund_otu.size();
  m++
)
{

  nll -= dnorm(
    b_abund_otu(m),
    Type(0),
    sd_abund_otu,
    true
  );
}

}




if (use_sample_re == 1)
{

for (
  int j = 0;
  j < b_sample.size();
  j++
)
{

  nll -= dnorm(
    b_sample(j),
    Type(0),
    sd_sample,
    true
  );
}

}




if (use_sample_otu_re == 1)
{

for (
  int q = 0;
  q < b_sample_otu.size();
  q++
)
{

  nll -= dnorm(
    b_sample_otu(q),
    Type(0),
    sd_sample_otu,
    true
  );
}

}




// ----------------------------------------------------------
// NEW: zero-inflation OTU random effects
// ----------------------------------------------------------

if (
use_zi_otu == 1 &&
(
family_code == 2 ||
family_code == 3
)
)
{

for (
  int m = 0;
  m < b_zi_otu.size();
  m++
)
{

  nll -= dnorm(
    b_zi_otu(m),
    Type(0),
    sd_zi_otu,
    true
  );

  // Assumption:
  //
  //     b_zi_m
  //     ~ Normal(
  //         0,
  //         sigma_zi^2
  //       )
  //
  // This allows each OTU to have a different
  // PCR-level dropout probability.
}

}




// ==========================================================
// 9. OCCUPANCY LINEAR PREDICTOR
// ==========================================================

vector<Type> eta_occ =
X_occ *
beta_occ;

// Fixed-effect component:
//
// eta_occ_im
// =
// X_occ_im beta_occ.

if (use_occ_otu == 1)
{

for (
  int g = 0;
  g < n_site_groups;
  g++
)
{

  eta_occ(g) +=
    b_occ_otu(
      occ_otu(g)
    );

  // Full occupancy model:
  //
  // logit(psi_im)
  // =
  // X_occ_im beta_occ
  // +
  // b_occ_m.
}

}




// ==========================================================
// 10. CAPTURE LINEAR PREDICTOR
// ==========================================================

vector<Type> eta_cap =
X_cap *
beta_cap;

if (use_cap_otu == 1)
{

for (
  int s = 0;
  s < n_sample_groups;
  s++
)
{

  eta_cap(s) +=
    b_cap_otu(
      sample_otu(s)
    );

  // logit(p_ijm)
  // =
  // X_cap beta_cap
  // +
  // b_cap_m.
}

}




// ==========================================================
// 11. ABUNDANCE LINEAR PREDICTOR
// ==========================================================

vector<Type> eta_abund =
X_abund *
beta_abund;

for (
int r = 0;
r < y.size();
r++
)
{

// --------------------------------------------------------
// Offset
// --------------------------------------------------------

eta_abund(r) +=
  offset_abund(r);


// --------------------------------------------------------
// OTU abundance random effect
// --------------------------------------------------------

if (use_abund_otu == 1)
{

  eta_abund(r) +=
    b_abund_otu(
      row_otu(r)
    );
}


// --------------------------------------------------------
// Sample abundance random effect
// --------------------------------------------------------

if (use_sample_re == 1)
{

  eta_abund(r) +=
    b_sample(
      row_sample_id(r)
    );
}


// --------------------------------------------------------
// Sample x OTU abundance random effect
// --------------------------------------------------------

if (use_sample_otu_re == 1)
{

  eta_abund(r) +=
    b_sample_otu(
      row_sample_otu_re(r)
    );
}

}




// ==========================================================
// 12. ZERO-INFLATION LINEAR PREDICTOR
// ==========================================================
//
// NEW MODEL:
//
// logit(pi_m)
// =
// zi_intercept
// +
// b_zi_otu(m)
//
// where:
//
// b_zi_otu(m)
// ~ N(0, sigma_zi^2).
//
// We construct one pi per OTU.
// ==========================================================

vector<Type> pi_zi_otu(
b_zi_otu.size()
);

for (
int m = 0;
m < b_zi_otu.size();
m++
)
{

Type eta_zi =
  zi_intercept;


if (
  use_zi_otu == 1 &&
  (
    family_code == 2 ||
    family_code == 3
  )
)
{

  eta_zi +=
    b_zi_otu(m);
}


pi_zi_otu(m) =
  Type(1) /
  (
    Type(1) +
    exp(
      -eta_zi
    )
  );

// pi_zi_otu(m)
//
// is the OTU-specific structural-zero probability.

}




// ==========================================================
// 13. READ-LEVEL COUNT LOG LIKELIHOOD
// ==========================================================

vector<Type> log_count(
y.size()
);

for (
int r = 0;
r < y.size();
r++
)
{

// --------------------------------------------------------
// Mean abundance
// --------------------------------------------------------

Type mu =
  exp(
    eta_abund(r)
  );

// Because abundance uses a log link:
//
//     lambda_ijkm
//       =
//       exp(eta_abund_ijkm).


Type lp =
  Type(0);


// --------------------------------------------------------
// Get OTU-specific zero-inflation probability
// --------------------------------------------------------

Type pi_zi =
  pi_zi_otu(
    row_otu(r)
  );

// Each observation uses the dropout probability
// corresponding to its OTU.



// ========================================================
// POISSON
// ========================================================

if (
  family_code == 0
)
{

  lp =
    dpois(
      y(r),
      mu,
      true
    );

  // log P(Y=y | mu).
}



// ========================================================
// NEGATIVE BINOMIAL 2
// ========================================================

else if (
  family_code == 1
)
{

  Type variance =
    mu +
    mu *
    mu /
    theta;


  lp =
    dnbinom2(
      y(r),
      mu,
      variance,
      true
    );

  // NB2:
  //
  // E(Y) = mu
  //
  // Var(Y)
  // =
  // mu + mu^2/theta.
}



// ========================================================
// ZERO-INFLATED POISSON
// ========================================================

else if (
  family_code == 2
)
{

  Type log_count_component =
    dpois(
      y(r),
      mu,
      true
    );


  // ------------------------------------------------------
  // Observed zero
  // ------------------------------------------------------

  if (
    y(r) ==
    Type(0)
  )
  {

    lp =
      logspace_add(

        log(
          pi_zi
        ),

        log(
          Type(1) -
          pi_zi
        ) +
        log_count_component

      );

    // For y = 0:
    //
    // P(Y=0)
    //
    // =
    //
    // pi_m
    //
    // +
    //
    // (1-pi_m)
    // P_Poisson(Y=0 | mu).
    //
    // pi_m is now OTU-specific.
  }


  // ------------------------------------------------------
  // Positive observation
  // ------------------------------------------------------

  else
  {

    lp =
      log(
        Type(1) -
        pi_zi
      ) +
      log_count_component;

    // A positive count cannot originate from the
    // structural-zero component.
  }
}



// ========================================================
// ZERO-INFLATED NEGATIVE BINOMIAL 2
// ========================================================

else if (
  family_code == 3
)
{

  Type variance =
    mu +
    mu *
    mu /
    theta;


  Type log_count_component =
    dnbinom2(
      y(r),
      mu,
      variance,
      true
    );


  if (
    y(r) ==
    Type(0)
  )
  {

    lp =
      logspace_add(

        log(
          pi_zi
        ),

        log(
          Type(1) -
          pi_zi
        ) +
        log_count_component

      );

    // ZINB zero:
    //
    // P(Y=0)
    //
    // =
    //
    // pi_m
    //
    // +
    //
    // (1-pi_m)
    // P_NB(Y=0 | mu, theta).
  }


  else
  {

    lp =
      log(
        Type(1) -
        pi_zi
      ) +
      log_count_component;
  }
}


// --------------------------------------------------------
// Store observation-level count log likelihood
// --------------------------------------------------------

log_count(r) =
  lp;

}




// ==========================================================
// 14. SAMPLE x OTU CONDITIONAL COUNT LIKELIHOOD
// ==========================================================
//
// First combine PCR replicates conditional on A = 1:
//
// L_count_ijm
//
// =
//
// product_k
// f(y_ijkm | A_ijm = 1)
//
// On log scale:
//
// log L_count_ijm
//
// =
//
// sum_k
// log f(y_ijkm).
// ==========================================================

vector<Type> log_sample_count(
n_sample_groups
);

log_sample_count.setZero();

for (
int r = 0;
r < y.size();
r++
)
{

int s =
  row_sample_group(r);


log_sample_count(s) +=
  log_count(r);

}




// ==========================================================
// 15. CAPTURE MARGINALIZATION
// ==========================================================

vector<Type> log_sample(
n_sample_groups
);

for (
int s = 0;
s < n_sample_groups;
s++
)
{

Type log_p =
  log_invlogit(
    eta_cap(s)
  );


Type log_1mp =
  log1m_invlogit(
    eta_cap(s)
  );


// --------------------------------------------------------
// Positive biological sample
// --------------------------------------------------------
//
// If at least one PCR observation is positive:
//
//     A_ijm = 1
//
// must hold.
//
// Therefore:
//
// L_ijm
//
// =
//
// p_ijm
// *
// product_k f(y_ijkm).
// --------------------------------------------------------


if (
  sample_positive(s) == 1
)
{

  log_sample(s) =
    log_p +
    log_sample_count(s);
}


// --------------------------------------------------------
// Entire biological sample is zero
// --------------------------------------------------------
//
// There are TWO possibilities:
//
// 1. Capture failed:
//
//        A = 0
//
//        probability = 1-p
//
// 2. Capture succeeded:
//
//        A = 1
//
//    but all PCR observations are zero.
//
// Therefore:
//
// L_ijm
//
// =
//
// (1-p_ijm)
//
// +
//
// p_ijm
// product_k f(0).
//
// f(0) already contains any ZIP/ZINB structural-zero
// probability pi_m.
// --------------------------------------------------------


else
{

  log_sample(s) =
    logspace_add(

      log_1mp,

      log_p +
      log_sample_count(s)

    );
}

}




// ==========================================================
// 16. SITE x OTU CONDITIONAL LIKELIHOOD GIVEN Z = 1
// ==========================================================
//
// Combine all biological sample likelihoods belonging
// to the same site x OTU:
//
// L(data_im | Z_im=1)
//
// =
//
// product_j L_ijm.
// ==========================================================

vector<Type> log_site_conditional(
n_site_groups
);

log_site_conditional.setZero();

for (
int s = 0;
s < n_sample_groups;
s++
)
{

int g =
  sample_site_group(s);


log_site_conditional(g) +=
  log_sample(s);

}




// ==========================================================
// 17. OCCUPANCY MARGINALIZATION
// ==========================================================

vector<Type> psi(
n_site_groups
);

vector<Type> log_site(
n_site_groups
);

for (
int g = 0;
g < n_site_groups;
g++
)
{

// --------------------------------------------------------
// Occupancy probability
// --------------------------------------------------------

psi(g) =
  Type(1) /
  (
    Type(1) +
    exp(
      -eta_occ(g)
    )
  );


Type log_psi =
  log_invlogit(
    eta_occ(g)
  );


Type log_1mpsi =
  log1m_invlogit(
    eta_occ(g)
  );


// --------------------------------------------------------
// Positive site x OTU history
// --------------------------------------------------------
//
// If any observation is positive:
//
//     Z_im = 1
//
// is required.
//
// Therefore:
//
// L_im
//
// =
//
// psi_im
//
// *
//
// L(data_im | Z_im=1).
// --------------------------------------------------------


if (
  site_positive(g) == 1
)
{

  log_site(g) =
    log_psi +
    log_site_conditional(g);
}


// --------------------------------------------------------
// Entire site x OTU history is zero
// --------------------------------------------------------
//
// Two explanations are possible:
//
// 1. Species absent:
//
//        Z_im = 0
//
//        probability = 1-psi
//
// 2. Species present:
//
//        Z_im = 1
//
//    but all biological samples/PCR reads are zero.
//
// Therefore:
//
// L_im
//
// =
//
// (1-psi_im)
//
// +
//
// psi_im
// L(all-zero data | Z_im=1).
// --------------------------------------------------------


else
{

  log_site(g) =
    logspace_add(

      log_1mpsi,

      log_psi +
      log_site_conditional(g)

    );
}


// --------------------------------------------------------
// Add site x OTU contribution to total NLL
// --------------------------------------------------------

nll -=
  log_site(g);

}




// ==========================================================
// 18. DERIVED QUANTITIES
// ==========================================================

vector<Type> capture_prob(
n_sample_groups
);

for (
int s = 0;
s < n_sample_groups;
s++
)
{

capture_prob(s) =
  Type(1) /
  (
    Type(1) +
    exp(
      -eta_cap(s)
    )
  );

// p_ijm =
// invlogit(eta_cap_ijm).

}




vector<Type> lambda(
y.size()
);

for (
int r = 0;
r < y.size();
r++
)
{

lambda(r) =
  exp(
    eta_abund(r)
  );

// lambda_ijkm =
// exp(eta_abund_ijkm).

}




// ==========================================================
// 19. REPORT RESULTS TO R
// ==========================================================

REPORT(
eta_occ
);

REPORT(
eta_cap
);

REPORT(
eta_abund
);

REPORT(
psi
);

REPORT(
capture_prob
);

REPORT(
lambda
);

REPORT(
pi_zi_otu
);

// NEW:
//
// Returns estimated structural-zero probability for every OTU.
//
// Example in R:
//
// fit$report$pi_zi_otu

REPORT(
log_sample
);

REPORT(
log_site
);




// ==========================================================
// 20. DELTA-METHOD DERIVED-PARAMETER STANDARD ERRORS
// ==========================================================

ADREPORT(
psi
);

// Approximate SEs for occupancy probabilities.

ADREPORT(
capture_prob
);

// Approximate SEs for capture probabilities.

ADREPORT(
pi_zi_otu
);

// NEW:
//
// Approximate SEs for OTU-specific zero-inflation
// probabilities using TMB automatic differentiation.

// ==========================================================
// Return negative log likelihood
// ==========================================================

return nll;
}

Close
