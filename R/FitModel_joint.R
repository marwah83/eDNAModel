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
        
        DLL = "eDNAModel",
        
        verbose = TRUE
) {
    
    abundance_family <-
        match.arg(
            abundance_family
        )
    
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
    # 2. Validate columns
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
    # 4. OTU filtering
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
    
    
    retained_otus <- otu_stats |>
        dplyr::filter(
            
            .data$total_count >=
                min_species_sum,
            
            .data$detected_replicates >=
                min_detection_replicates
            
        ) |>
        dplyr::pull(
            .data[[otu_col]]
        )
    
    
    dat <- dat |>
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
    # 5. Explicit factor indices
    # ============================================================
    
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
    
    
    dat$.site_otu <-
        interaction(
            dat$.site,
            dat$.otu,
            drop = TRUE,
            lex.order = TRUE
        )
    
    
    dat$.sample_otu <-
        interaction(
            dat$.site,
            dat$.sample,
            dat$.otu,
            drop = TRUE,
            lex.order = TRUE
        )
    
    
    # ============================================================
    # 6. Build sample x OTU table
    # ============================================================
    
    sample_df <- dat |>
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
    # 7. Carry capture-formula variables
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
        
        cap_covars <- dat |>
            dplyr::select(
                .data$.sample_otu,
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
    # 8. Build site x OTU table
    # ============================================================
    
    site_df <- dat |>
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
    # 9. Carry occupancy-formula variables
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
        
        occ_covars <- dat |>
            dplyr::select(
                .data$.site_otu,
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
    # 10. Design matrices
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
    # 11. Abundance offset
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
    # 12. Rebuild explicit group levels
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
    # 13. Indices
    # ============================================================
    
    row_sample_group <-
        match(
            as.character(
                dat$.sample_otu
            ),
            sample_group_levels
        ) -
        1L
    
    
    sample_site_group <-
        match(
            as.character(
                sample_df$.site_otu
            ),
            site_group_levels
        ) -
        1L
    
    
    occ_otu <-
        as.integer(
            site_df$.otu
        ) -
        1L
    
    
    sample_otu <-
        as.integer(
            sample_df$.otu
        ) -
        1L
    
    
    row_otu <-
        as.integer(
            dat$.otu
        ) -
        1L
    
    
    row_sample_id <-
        as.integer(
            dat$.sample
        ) -
        1L
    
    
    # One RE per sample x OTU
    row_sample_otu_re <-
        row_sample_group
    
    
    # ============================================================
    # 14. Validate mapping
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
    
    
    # ============================================================
    # 15. Family code
    # ============================================================
    
    family_code <- switch(
        
        abundance_family,
        
        poisson = 0L,
        
        nbinom = 1L,
        
        zip = 2L,
        
        zinb = 3L
    )
    
    
    # ============================================================
    # 16. TMB data list
    # ============================================================
    
    data_tmb <- list(
        
        y =
            dat[[count_col]],
        
        X_occ =
            X_occ,
        
        X_cap =
            X_cap,
        
        X_abund =
            X_abund,
        
        offset_abund =
            offset_abund,
        
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
            )
    )
    
    
    # ============================================================
    # 17. Dimensions
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
    # 18. Better starting values
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
        
        mean_count_start <- 1
    }
    
    
    beta_occ_start <-
        rep(
            0,
            ncol(
                X_occ
            )
        )
    
    beta_occ_start[1] <-
        stats::qlogis(
            site_positive_rate
        )
    
    
    beta_cap_start <-
        rep(
            0,
            ncol(
                X_cap
            )
        )
    
    beta_cap_start[1] <-
        stats::qlogis(
            sample_positive_rate
        )
    
    
    beta_abund_start <-
        rep(
            0,
            ncol(
                X_abund
            )
        )
    
    beta_abund_start[1] <-
        log(
            pmax(
                mean_count_start,
                1e-4
            )
        )
    
    
    parameters <- list(
        
        beta_occ =
            beta_occ_start,
        
        beta_cap =
            beta_cap_start,
        
        beta_abund =
            beta_abund_start,
        
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
        
        log_sd_occ_otu =
            log(0.5),
        
        log_sd_cap_otu =
            log(0.5),
        
        log_sd_abund_otu =
            log(0.5),
        
        log_sd_sample =
            log(0.5),
        
        log_sd_sample_otu =
            log(0.5),
        
        log_theta =
            log(10),
        
        zi_intercept =
            stats::qlogis(
                0.05
            )
    )
    
    
    # ============================================================
    # 19. CRITICAL:
    # MAP UNUSED PARAMETERS OUT OF OPTIMIZATION
    # ============================================================
    
    map <- list()
    
    
    # ------------------------------------------------------------
    # Occupancy OTU effect OFF
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
    # Capture OTU effect OFF
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
    # Abundance OTU effect OFF
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
    # Sample effect OFF
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
    # Sample x OTU effect OFF
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
    # Poisson has neither theta nor zero-inflation
    # ------------------------------------------------------------
    
    if (
        abundance_family ==
        "poisson"
    ) {
        
        map$log_theta <-
            factor(NA)
        
        map$zi_intercept <-
            factor(NA)
    }
    
    
    # ------------------------------------------------------------
    # NB requires theta but not zero inflation
    # ------------------------------------------------------------
    
    if (
        abundance_family ==
        "nbinom"
    ) {
        
        map$zi_intercept <-
            factor(NA)
    }
    
    
    # ------------------------------------------------------------
    # ZIP requires zero inflation but not theta
    # ------------------------------------------------------------
    
    if (
        abundance_family ==
        "zip"
    ) {
        
        map$log_theta <-
            factor(NA)
    }
    
    
    # ZINB:
    # log_theta and zi_intercept both estimated
    
    
    # ============================================================
    # 20. Random parameters to integrate using Laplace
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
    
    
    # ============================================================
    # 21. Optional debugging export
    # ============================================================
    
    if (verbose) {
        
        message(
            "Sites x OTUs: ",
            nrow(site_df)
        )
        
        message(
            "Samples x OTUs: ",
            nrow(sample_df)
        )
        
        message(
            "Read observations: ",
            nrow(dat)
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
                length(random_effects) == 0
            ) {
                "none"
            } else {
                paste(
                    random_effects,
                    collapse = ", "
                )
            }
        )
    }
    
    
    # ============================================================
    # 22. Make TMB objective
    # ============================================================
    
    obj <- TMB::MakeADFun(
        
        data =
            data_tmb,
        
        parameters =
            parameters,
        
        random =
            if (
                length(random_effects) == 0
            ) {
                NULL
            } else {
                random_effects
            },
        
        map =
            if (
                length(map) == 0
            ) {
                NULL
            } else {
                map
            },
        
        DLL =
            "eDNAModel",
        
        silent =
            !verbose
    )
    
    
    # ============================================================
    # 23. Check initial objective
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
    # 24. Optimize
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
    }
    
    
    # ============================================================
    # 25. Standard errors from marginal likelihood
    # ============================================================
    
    # Make sure random effects are at their conditional mode
    # for the final outer parameter estimates.
    invisible(
        obj$fn(opt$par)
    )
    
    sdr <- tryCatch(
        
        TMB::sdreport(
            obj,
            par.fixed = opt$par,
            getJointPrecision = length(random_effects) > 0
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
    # 26. Model reports
    # ============================================================
    
    # IMPORTANT:
    # With TMB random effects, opt$par contains only the
    # outer/fixed parameters. obj$report() should therefore
    # use TMB's internally stored complete parameter vector.
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
    # 27. Fixed-effect / parameter summary
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
    # 28. ADREPORT summary
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
    # 29. Hessian / SE checks
    # ============================================================
    
    finite_se <- if (
        nrow(fixed_summary) > 0 &&
        "Std. Error" %in%
        names(fixed_summary)
    ) {
        
        all(
            is.finite(
                fixed_summary[["Std. Error"]]
            )
        )
        
    } else {
        
        NA
    }
    
    
    # ============================================================
    # 30. Return
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
        
        data_tmb =
            data_tmb,
        
        parameters_start =
            parameters,
        
        parameter_map =
            map,
        
        random_effects =
            random_effects,
        
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
        
        convergence = list(
            
            code =
                opt$convergence,
            
            message =
                opt$message,
            
            objective =
                opt$objective,
            
            finite_standard_errors =
                finite_se
        ),
        
        note = paste(
            "Joint observed-data likelihood.",
            "Latent occupancy Z and biological-sample capture A",
            "are analytically marginalized.",
            "Enabled Gaussian random effects are integrated using",
            "the TMB Laplace approximation."
        )
    )
}
