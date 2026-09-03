#' Fit a Joint Hierarchical eDNA Occupancy--Capture--Abundance Model
#'
#' @description
#' Fits a joint hierarchical model for environmental DNA (eDNA) data that
#' separates ecological occupancy, biological-sample capture, and sequencing
#' read abundance. The discrete latent occupancy and capture states are
#' analytically marginalized from the observed-data likelihood, while optional
#' Gaussian random effects are integrated using the Laplace approximation
#' implemented in \pkg{TMB}.
#'
#' The model contains three linked components:
#'
#' \deqn{
#' Z_{im} \sim \mathrm{Bernoulli}(\psi_{im}),
#' }
#'
#' \deqn{
#' \mathrm{logit}(\psi_{im})
#' =
#' \mathbf{x}^{(\psi)\top}_{i}\boldsymbol{\beta}_{\psi}
#' +
#' b^{(\psi)}_{m},
#' }
#'
#' \deqn{
#' A_{ijm}\mid Z_{im}=1
#' \sim
#' \mathrm{Bernoulli}(p_{ijm}),
#' }
#'
#' \deqn{
#' \mathrm{logit}(p_{ijm})
#' =
#' \mathbf{x}^{(p)\top}_{ij}\boldsymbol{\beta}_{p}
#' +
#' b^{(p)}_{m},
#' }
#'
#' and
#'
#' \deqn{
#' \log(\lambda_{ijkm})
#' =
#' \mathbf{x}^{(\lambda)\top}_{ijk}
#' \boldsymbol{\beta}_{\lambda}
#' +
#' b^{(\lambda)}_{m}
#' +
#' u_j
#' +
#' v_{jm}
#' +
#' o_{ijk}.
#' }
#'
#' Here, \eqn{i} indexes sites or site--time units, \eqn{j} biological
#' samples, \eqn{k} technical replicates, and \eqn{m} OTUs or taxa.
#'
#' @param phyloseq A \code{phyloseq} object containing an OTU/count table
#'   and corresponding sample metadata.
#'
#' @param site_col Character string giving the column in
#'   \code{sample_data(phyloseq)} defining the ecological site or
#'   site--time occupancy unit.
#'
#' @param sample_col Character string giving the biological-sample identifier.
#'   Default is \code{"Name"}.
#'
#' @param replicate_col Optional character string giving the technical or PCR
#'   replicate identifier. Default is \code{NULL}.
#'
#' @param otu_col Character string giving the OTU/taxon identifier used in the
#'   long-format model data. Default is \code{"OTU"}.
#'
#' @param count_col Character string giving the observed sequencing-read count
#'   column. Default is \code{"y"}.
#'
#' @param occupancy_formula One-sided formula specifying fixed effects in the
#'   occupancy component. Default is \code{~ 1}.
#'
#' @param capture_formula One-sided formula specifying fixed effects in the
#'   biological-sample capture component. Default is \code{~ 1}.
#'
#' @param abundance_formula One-sided formula specifying fixed effects in the
#'   sequencing-read abundance component. Default is \code{~ 1}.
#'
#' @param abundance_offset Optional character string naming a positive variable
#'   to be included as a log-offset in the abundance linear predictor.
#'   Default is \code{NULL}.
#'
#' @param abundance_family Character string specifying the conditional
#'   read-count distribution. One of \code{"poisson"}, \code{"nbinom"},
#'   \code{"zip"}, or \code{"zinb"}.
#'
#' @param min_species_sum Minimum total read count required for an OTU to be
#'   retained in the analysis. Default is \code{10}.
#'
#' @param min_detection_replicates Minimum number of positive replicate-level
#'   observations required for an OTU to be retained. Default is \code{1}.
#'
#' @param random_occ_otu Logical. If \code{TRUE}, include an OTU-specific
#'   Gaussian random intercept in the occupancy model.
#'
#' @param random_capture_otu Logical. If \code{TRUE}, include an OTU-specific
#'   Gaussian random intercept in the capture model.
#'
#' @param random_abund_otu Logical. If \code{TRUE}, include an OTU-specific
#'   Gaussian random intercept in the abundance model.
#'
#' @param random_sample Logical. If \code{TRUE}, include a biological-sample
#'   Gaussian random intercept in the abundance model.
#'
#' @param random_sample_otu Logical. If \code{TRUE}, include a sample-by-OTU
#'   Gaussian random effect in the abundance model.
#'
#' @param DLL Character string specifying the TMB dynamic library containing
#'   the joint likelihood template. For the installed \pkg{eDNAModel} package,
#'   this should normally be \code{"eDNAModel"}.
#'
#' @param verbose Logical. If \code{TRUE}, print model dimensions,
#'   optimization progress, convergence information, and post-fit calculations.
#'
#' @details
#' The observed-data likelihood is obtained by analytically summing over the
#' latent biological-sample capture state and latent site-level occupancy state.
#' Consequently, these binary latent variables are not imputed during model
#' fitting.
#'
#' For an all-zero biological sample, the contribution conditional on site
#' occupancy is
#'
#' \deqn{
#' (1-p_{ijm})
#' +
#' p_{ijm}
#' \prod_k f(0;\lambda_{ijkm}),
#' }
#'
#' whereas a biological sample containing at least one positive read has
#' likelihood
#'
#' \deqn{
#' p_{ijm}
#' \prod_k
#' f(Y_{ijkm};\lambda_{ijkm}).
#' }
#'
#' At the site--OTU level, if the complete observation history is zero,
#'
#' \deqn{
#' (1-\psi_{im})
#' +
#' \psi_{im}
#' \prod_j L_{ijm}^{(1)}
#' }
#'
#' is used. If any read count is positive, occupancy is necessarily present and
#' the contribution is
#'
#' \deqn{
#' \psi_{im}
#' \prod_j L_{ijm}^{(1)}.
#' }
#'
#' Enabled Gaussian random effects are integrated from the likelihood using
#' TMB's Laplace approximation. Fixed-effect and variance-component uncertainty
#' is obtained from the Hessian of the marginal likelihood through
#' \code{TMB::sdreport()}.
#'
#' The negative-binomial model uses the NB2 parameterization
#'
#' \deqn{
#' \mathrm{Var}(Y)
#' =
#' \mu + \frac{\mu^2}{\theta}.
#' }
#'
#' For ZIP and ZINB models, the parameter \eqn{\pi} represents an additional
#' structural-zero probability conditional on biological-sample capture.
#'
#' @return A list containing:
#'
#' \describe{
#'   \item{\code{fit}}{The \code{nlminb} optimization result.}
#'   \item{\code{tmb_object}}{The fitted TMB objective object.}
#'   \item{\code{sdreport}}{The \code{TMB::sdreport} object containing
#'     likelihood-based standard errors.}
#'   \item{\code{report}}{Derived quantities reported by the C++ TMB model.}
#'   \item{\code{fixed_effects}}{Estimated fixed parameters and standard errors.}
#'   \item{\code{derived}}{Delta-method summaries for reported quantities such
#'     as occupancy and capture probabilities.}
#'   \item{\code{site_data}}{Site-by-OTU model data.}
#'   \item{\code{sample_data}}{Biological-sample-by-OTU model data.}
#'   \item{\code{long_df}}{Long-format count data used for fitting.}
#'   \item{\code{otu_stats}}{OTU filtering statistics.}
#'   \item{\code{retained_otus}}{Names of OTUs retained after filtering.}
#'   \item{\code{random_effects}}{Names of Gaussian random effects integrated
#'     using the Laplace approximation.}
#'   \item{\code{formulas}}{Occupancy, capture, and abundance formulas.}
#'   \item{\code{abundance_family}}{Conditional count distribution used.}
#'   \item{\code{convergence}}{Optimizer convergence information and
#'     standard-error diagnostics.}
#' }
#'
#' @examples
#' \dontrun{
#' data("physeq_new", package = "eDNAModel")
#'
#' physeq_one <-
#'   physeq_new[["Marine Invasive species Trondheim"]]
#'
#' meta <- as.data.frame(
#'   phyloseq::sample_data(physeq_one)
#' )
#'
#' meta$site_month <- interaction(
#'   meta$Sampling.area.Name,
#'   meta$Samplingmonth,
#'   drop = TRUE,
#'   sep = "_"
#' )
#'
#' phyloseq::sample_data(physeq_one) <-
#'   phyloseq::sample_data(meta)
#'
#' fit <- FitModel_joint(
#'   phyloseq = physeq_one,
#'   site_col = "site_month",
#'   sample_col = "Name",
#'   replicate_col = "Replicate",
#'   otu_col = "OTU",
#'   count_col = "y",
#'   occupancy_formula = ~ Samplingmonth,
#'   capture_formula = ~ 1,
#'   abundance_formula = ~ Samplingmonth + Replicate,
#'   abundance_family = "zinb",
#'   random_occ_otu = TRUE,
#'   random_capture_otu = TRUE,
#'   random_abund_otu = TRUE,
#'   random_sample = TRUE,
#'   random_sample_otu = TRUE
#' )
#'
#' fit$convergence
#' fit$fixed_effects
#' }
#'
#' @references
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., and Bell, B. M.
#' (2016). TMB: Automatic differentiation and Laplace approximation.
#' \emph{Journal of Statistical Software}, 70(5), 1--21.
#'
#' @importFrom stats model.matrix nlminb qlogis
#' @importFrom dplyr group_by summarise filter pull select left_join across
#' @importFrom TMB MakeADFun sdreport
#' @export
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
