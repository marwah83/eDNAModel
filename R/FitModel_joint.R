#' Fit a Joint Hierarchical eDNA Occupancy--Capture--Abundance Model
#'
#' @description
#' Fits a joint hierarchical model for environmental DNA (eDNA) metabarcoding
#' count data using Template Model Builder (TMB). The model represents the
#' observation process through linked occupancy, capture/detection, and
#' sequencing-abundance components.
#'
#' The model can be fitted with Poisson, negative-binomial, zero-inflated
#' Poisson (ZIP), or zero-inflated negative-binomial (ZINB) abundance
#' distributions. OTU-level random effects can be included separately in the
#' occupancy, capture, abundance, and zero-inflation components. Additional
#' sample and sample-by-OTU random effects can also be included.
#'
#' Latent occupancy and biological-sample capture states are analytically
#' marginalized in the observed-data likelihood. Enabled Gaussian random
#' effects are integrated using the Laplace approximation implemented in TMB.
#'
#' @details
#' The hierarchical model contains three principal ecological/observation
#' components.
#'
#' **Occupancy model**
#'
#' For site \eqn{i} and OTU \eqn{k}, occupancy probability is modeled as
#'
#' \deqn{
#' \mathrm{logit}(\psi_{ik})
#' =
#' X^{(\psi)}_{ik}\beta_{\psi}
#' +
#' b^{(\psi)}_k,
#' }
#'
#' where \eqn{X^{(\psi)}} is the occupancy design matrix,
#' \eqn{\beta_{\psi}} contains fixed effects, and
#' \eqn{b^{(\psi)}_k} is an optional OTU-level random intercept.
#'
#' **Capture/detection model**
#'
#' Conditional on occupancy, the probability that OTU \eqn{k} is captured in
#' biological sample \eqn{j} is modeled as
#'
#' \deqn{
#' \mathrm{logit}(p_{jk})
#' =
#' X^{(p)}_{jk}\beta_p
#' +
#' b^{(p)}_k.
#' }
#'
#' **Abundance model**
#'
#' Conditional on capture, expected read abundance is modeled as
#'
#' \deqn{
#' \log(\lambda_{rjk})
#' =
#' X^{(\lambda)}_{rjk}\beta_{\lambda}
#' +
#' b^{(\lambda)}_k
#' +
#' b^{(s)}_j
#' +
#' b^{(s \times k)}_{jk}
#' +
#' o_{rjk},
#' }
#'
#' where \eqn{o_{rjk}} is an optional log-scale abundance offset.
#'
#' For ZIP and ZINB models, additional structural zeros are represented using
#' a zero-inflation probability
#'
#' \deqn{
#' \mathrm{logit}(\pi_k)
#' =
#' \alpha_{\pi}
#' +
#' b^{(\pi)}_k,
#' }
#'
#' where \eqn{\alpha_{\pi}} is the population-level zero-inflation intercept
#' and \eqn{b^{(\pi)}_k} is an optional OTU-level random effect. This allows
#' structural-zero probabilities to vary among OTUs rather than assuming a
#' common dropout probability for all taxa.
#'
#' @section Numerical convergence:
#'
#' Numerical convergence is assessed using several complementary diagnostics
#' rather than relying solely on the optimizer return code.
#'
#' A fit is classified as strictly numerically converged only when all of the
#' following conditions are satisfied:
#'
#' \itemize{
#'   \item \code{nlminb()} returns convergence code 0;
#'   \item all components of the final outer TMB gradient are finite;
#'   \item the maximum absolute outer gradient is no greater than
#'         \code{gradient_tol};
#'   \item \code{TMB::sdreport()} completes successfully;
#'   \item the Hessian is positive definite according to
#'         \code{sdreport()$pdHess}; and
#'   \item all standard errors returned for the fixed parameter vector are
#'         finite.
#' }
#'
#' The default gradient tolerance is \eqn{10^{-3}}. The gradient criterion is
#'
#' \deqn{
#' \max_j
#' \left|
#' \frac{\partial \tilde{\ell}(\theta)}
#'      {\partial \theta_j}
#' \right|
#' \leq
#' \mathrm{gradient\_tol},
#' }
#'
#' where \eqn{\tilde{\ell}(\theta)} denotes the TMB objective after Laplace
#' approximation over enabled random effects.
#'
#' If the optimizer returns convergence code 0 but the gradient remains above
#' the specified tolerance, the fit is not classified as strictly numerically
#' converged. The function can automatically restart optimization from the
#' previous solution. Restarts terminate when the convergence criteria are
#' satisfied, the maximum number of restarts is reached, or changes in the
#' objective and gradient indicate that optimization has stalled.
#'
#' For user-facing diagnostics, gradient results may additionally be labelled
#' \code{"PASS"}, \code{"MARGINAL"}, or \code{"FAIL"}. The
#' \code{"MARGINAL"} classification is descriptive only and does not override
#' the strict gradient convergence criterion.
#'
#' @section Hessian diagnostic:
#'
#' Positive definiteness of the Hessian is assessed using the
#' \code{pdHess} diagnostic returned by \code{TMB::sdreport()}. A
#' non-positive-definite Hessian may indicate a flat likelihood direction,
#' weak or non-identification of one or more parameters, a saddle point, or
#' other numerical difficulties.
#'
#' This diagnostic evaluates positive definiteness rather than merely positive
#' semidefiniteness. It is considered jointly with the optimizer status and
#' gradient diagnostic; \code{pdHess = TRUE} alone is not interpreted as
#' evidence of convergence.
#'
#' @section Standard-error diagnostic:
#'
#' Fixed-parameter standard errors are obtained from
#' \code{summary(sdreport_object, "fixed")}. They pass the finite-standard-error
#' diagnostic only when every returned standard error satisfies
#' \code{is.finite()}. Consequently, \code{NA}, \code{NaN}, \code{Inf}, and
#' \code{-Inf} standard errors cause this diagnostic to fail.
#'
#' This is a numerical sanity check and does not imply that finite standard
#' errors are necessarily small or scientifically informative. Large but
#' finite standard errors are handled separately as heuristic parameter
#' warnings.
#'
#' @section Heuristic parameter diagnostics:
#'
#' In addition to the formal numerical convergence criteria, the function can
#' flag potentially problematic parameter estimates. These diagnostics are
#' intended to assist interpretation and do not by themselves cause a fit to
#' be classified as non-converged.
#'
#' Two heuristic checks are currently used:
#'
#' \itemize{
#'   \item estimated log-standard-deviation parameters below
#'         \code{log_sd_warning_threshold}, which may indicate a variance
#'         component close to zero; and
#'   \item unusually large standard errors relative to the absolute parameter
#'         estimate, based on \code{se_estimate_ratio_threshold}.
#' }
#'
#' These thresholds are diagnostic guidelines rather than universal
#' statistical criteria and should be interpreted in the context of the fitted
#' model and data.
#'
#' @section Memory use and sdreport:
#'
#' For large eDNA OTU matrices, covariance calculations performed by
#' \code{TMB::sdreport()} can require substantial memory. Therefore,
#' \code{get_report_covariance = FALSE} is the default.
#'
#' Users requiring the covariance matrix of reported quantities can explicitly
#' set \code{get_report_covariance = TRUE}. This option may substantially
#' increase memory requirements and should be used with care for large
#' datasets.
#'
#' @param phyloseq A \code{phyloseq} object containing the OTU count table and
#'   associated sample metadata.
#'
#' @param site_col Character string giving the sample-data column identifying
#'   sampling sites.
#'
#' @param sample_col Character string giving the column identifying biological
#'   samples. Default is \code{"Name"}.
#'
#' @param replicate_col Optional character string identifying technical or PCR
#'   replicates. Default is \code{NULL}.
#'
#' @param otu_col Character string identifying the OTU column in the long-format
#'   data. Default is \code{"OTU"}.
#'
#' @param count_col Character string identifying the observed read-count
#'   column. Default is \code{"y"}.
#'
#' @param occupancy_formula A model formula specifying fixed effects for the
#'   occupancy component. Default is \code{~ 1}.
#'
#' @param capture_formula A model formula specifying fixed effects for the
#'   capture/detection component. Default is \code{~ 1}.
#'
#' @param abundance_formula A model formula specifying fixed effects for the
#'   abundance component. Default is \code{~ 1}.
#'
#' @param abundance_offset Optional character string identifying a positive
#'   abundance-exposure variable. The logarithm of this variable is included
#'   as an offset. Default is \code{NULL}.
#'
#' @param abundance_family Character string specifying the abundance
#'   distribution. One of \code{"poisson"}, \code{"nbinom"}, \code{"zip"}, or
#'   \code{"zinb"}.
#'
#' @param min_species_sum Minimum total read count required for an OTU to be
#'   retained. Default is 10.
#'
#' @param min_detection_replicates Minimum number of positive observations
#'   required for an OTU to be retained. Default is 1.
#'
#' @param random_occ_otu Logical. Include an OTU-level random intercept in the
#'   occupancy model. Default is \code{TRUE}.
#'
#' @param random_capture_otu Logical. Include an OTU-level random intercept in
#'   the capture model. Default is \code{TRUE}.
#'
#' @param random_abund_otu Logical. Include an OTU-level random intercept in the
#'   abundance model. Default is \code{TRUE}.
#'
#' @param random_sample Logical. Include a biological-sample random intercept
#'   in the abundance model. Default is \code{TRUE}.
#'
#' @param random_sample_otu Logical. Include a sample-by-OTU random intercept
#'   in the abundance model. Default is \code{FALSE}.
#'
#' @param random_zi_otu Logical. Include an OTU-level random intercept in the
#'   zero-inflation model. This option is used only for ZIP and ZINB models.
#'   Default is \code{TRUE}.
#'
#' @param gradient_tol Positive numeric value specifying the maximum acceptable
#'   absolute component of the final outer gradient for strict numerical
#'   convergence. Default is \code{1e-3}.
#'
#' @param gradient_marginal_factor Numeric value greater than 1 defining the
#'   optional user-facing marginal-gradient region. This affects diagnostic
#'   labels only and does not change the strict convergence criterion.
#'   Default is 2.
#'
#' @param max_restarts Non-negative integer giving the maximum number of
#'   automatic optimizer restarts. Default is 2.
#'
#' @param iter_max Maximum number of \code{nlminb()} iterations per
#'   optimization pass. Default is 5000.
#'
#' @param eval_max Maximum number of objective evaluations per optimization
#'   pass. Default is 10000.
#'
#' @param rel_tol Relative convergence tolerance supplied to \code{nlminb()}.
#'   Default is \code{1e-10}.
#'
#' @param stall_objective_tol Tolerance used to identify negligible changes in
#'   the objective between consecutive optimization passes.
#'
#' @param stall_gradient_tol Tolerance used to identify negligible changes in
#'   the maximum absolute gradient between consecutive optimization passes.
#'
#' @param n_gradient_report Number of parameters with the largest absolute
#'   gradients to report when the gradient criterion is not satisfied.
#'   Default is 2.
#'
#' @param get_report_covariance Logical. Passed to
#'   \code{TMB::sdreport(getReportCovariance = ...)}. The default is
#'   \code{FALSE} to reduce memory use for large OTU datasets.
#'
#' @param get_joint_precision Logical. Request the joint precision matrix from
#'   \code{TMB::sdreport()} when random effects are present. Default is
#'   \code{TRUE}. This may also increase computational and memory requirements
#'   for large models.
#'
#' @param log_sd_warning_threshold Numeric threshold used to flag very small
#'   estimated random-effect standard deviations on the log-SD scale. This is
#'   a heuristic diagnostic only. Default is -5.
#'
#' @param se_estimate_ratio_threshold Numeric threshold used to flag fixed
#'   parameters whose standard error is large relative to the absolute
#'   estimate. This is a heuristic diagnostic and does not determine
#'   convergence. Default is 10.
#'
#' @param estimate_zero_tol Positive numeric tolerance below which an estimate
#'   is treated as effectively zero when computing SE-to-estimate ratios.
#'   Default is \code{1e-6}.
#'
#' @param DLL Character string giving the name of the compiled TMB dynamic
#'   library. Default is \code{"eDNAModel"}.
#'
#' @param verbose Logical. If \code{TRUE}, print fitting progress and numerical
#'   diagnostics. Default is \code{TRUE}.
#'
#' @return
#' A list containing the fitted optimizer object, TMB object, \code{sdreport}
#' results, parameter summaries, processed model data, and numerical
#' diagnostics. Important components include:
#'
#' \describe{
#'   \item{\code{fit}}{The final \code{nlminb()} optimization result.}
#'   \item{\code{tmb_object}}{The fitted TMB objective object.}
#'   \item{\code{sdreport}}{The \code{TMB::sdreport()} object, when successful.}
#'   \item{\code{fixed_effects}}{Estimates and standard errors for fixed
#'     parameters.}
#'   \item{\code{derived}}{Reported/derived parameter summaries.}
#'   \item{\code{convergence}}{Detailed optimizer, gradient, Hessian,
#'     standard-error, restart, and overall numerical-convergence diagnostics.}
#'   \item{\code{diagnostics}}{User-facing diagnostic table and heuristic
#'     parameter warnings.}
#'   \item{\code{sdreport_settings}}{Settings used when calling
#'     \code{TMB::sdreport()}.}
#' }
#'
#' @note
#' The convergence diagnostics assess numerical optimization and local
#' curvature of the Laplace-approximated objective. They do not, by
#' themselves, establish the accuracy of the Laplace approximation.
#' Simulation-based parameter-recovery and coverage studies are recommended
#' when evaluating approximation accuracy, particularly for sparse and highly
#' discrete eDNA datasets.
#'
#' @examples
#' \dontrun{
#'
#' fit_zinb <- FitModel_joint(
#'     phyloseq = ps,
#'     site_col = "Sampling.area.Name",
#'     sample_col = "Name",
#'     replicate_col = "Replicate",
#'     otu_col = "OTU",
#'     count_col = "y",
#'     occupancy_formula = ~ 1,
#'     capture_formula = ~ 1,
#'     abundance_formula = ~ 1,
#'     abundance_family = "zinb",
#'     random_occ_otu = TRUE,
#'     random_capture_otu = TRUE,
#'     random_abund_otu = TRUE,
#'     random_zi_otu = TRUE,
#'     random_sample = TRUE,
#'     random_sample_otu = FALSE,
#'     gradient_tol = 1e-3,
#'     get_report_covariance = FALSE,
#'     verbose = TRUE
#' )
#'
#' fit_zinb$fixed_effects
#' fit_zinb$convergence$converged
#' fit_zinb$convergence$max_abs_gradient
#' fit_zinb$convergence$pd_hessian
#' fit_zinb$diagnostics$table
#' }
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats model.matrix nlminb qlogis
#' @importFrom dplyr group_by summarise filter select pull first across
#'   all_of left_join
#'
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

    abundance_family = c("poisson", "nbinom", "zip", "zinb"),

    min_species_sum = 10,
    min_detection_replicates = 1,

    # ----------------------------------------------------------
    # Built-in Gaussian random effects
    # ----------------------------------------------------------
    random_occ_otu = TRUE,
    random_capture_otu = TRUE,
    random_abund_otu = TRUE,
    random_sample = TRUE,
    random_sample_otu = FALSE,
    random_zi_otu = TRUE,

    # ----------------------------------------------------------
    # Numerical-convergence controls
    # ----------------------------------------------------------
    gradient_tol = 1e-3,
    gradient_marginal_factor = 2,
    max_restarts = 2L,

    iter_max = 5000L,
    eval_max = 10000L,
    rel_tol = 1e-10,

    # Stop restarting if neither the objective nor gradient improves.
    stall_objective_tol = 1e-8,
    stall_gradient_tol = 1e-8,

    n_gradient_report = 2L,

    # ----------------------------------------------------------
    # sdreport memory controls
    #
    # FALSE is deliberately the default for report covariance.
    # ----------------------------------------------------------
    get_report_covariance = FALSE,
    get_joint_precision = TRUE,

    # ----------------------------------------------------------
    # Heuristic parameter diagnostics
    # These DO NOT determine strict convergence.
    # ----------------------------------------------------------
    log_sd_warning_threshold = -5,
    se_estimate_ratio_threshold = 10,
    estimate_zero_tol = 1e-6,

    DLL = "eDNAModel",
    verbose = TRUE
) {

    # ==========================================================
    # 1. Match abundance family
    # ==========================================================

    abundance_family <- match.arg(abundance_family)


    # ==========================================================
    # 2. Validate control arguments
    # ==========================================================

    if (length(gradient_tol) != 1L ||
        !is.finite(gradient_tol) ||
        gradient_tol <= 0) {
        stop("gradient_tol must be a single positive finite number.")
    }

    if (length(gradient_marginal_factor) != 1L ||
        !is.finite(gradient_marginal_factor) ||
        gradient_marginal_factor <= 1) {
        stop("gradient_marginal_factor must be greater than 1.")
    }

    max_restarts <- as.integer(max_restarts)
    iter_max <- as.integer(iter_max)
    eval_max <- as.integer(eval_max)
    n_gradient_report <- as.integer(n_gradient_report)

    if (length(max_restarts) != 1L ||
        is.na(max_restarts) ||
        max_restarts < 0L) {
        stop("max_restarts must be a non-negative integer.")
    }

    if (length(iter_max) != 1L ||
        is.na(iter_max) ||
        iter_max < 1L) {
        stop("iter_max must be a positive integer.")
    }

    if (length(eval_max) != 1L ||
        is.na(eval_max) ||
        eval_max < 1L) {
        stop("eval_max must be a positive integer.")
    }

    if (length(n_gradient_report) != 1L ||
        is.na(n_gradient_report) ||
        n_gradient_report < 1L) {
        stop("n_gradient_report must be a positive integer.")
    }


    # ==========================================================
    # 3. Convert the phyloseq object to long format
    #
    # prepare_long_data() is assumed to be an internal
    # eDNAModel helper already defined in the package.
    # ==========================================================

    prep <- prepare_long_data(
        physeq_obj = phyloseq,
        site_col = site_col,
        nested_cols = unique(
            stats::na.omit(
                c(sample_col, replicate_col)
            )
        )
    )

    dat <- prep$long_df


    # ==========================================================
    # 4. Check required variables
    # ==========================================================

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
            c(required, replicate_col)
        )
    }

    missing_cols <- setdiff(
        required,
        names(dat)
    )

    if (length(missing_cols) > 0L) {
        stop(
            "Missing columns: ",
            paste(missing_cols, collapse = ", ")
        )
    }


    # ==========================================================
    # 5. Standardize basic variable types
    # ==========================================================

    dat[[site_col]] <- as.character(dat[[site_col]])
    dat[[sample_col]] <- as.character(dat[[sample_col]])
    dat[[otu_col]] <- as.character(dat[[otu_col]])
    dat[[count_col]] <- as.numeric(dat[[count_col]])

    if (!is.null(replicate_col)) {
        dat[[replicate_col]] <-
            as.character(dat[[replicate_col]])
    }


    # ==========================================================
    # 6. Remove observations with invalid counts
    # ==========================================================

    dat <- dat[
        !is.na(dat[[count_col]]) &
            is.finite(dat[[count_col]]),
        ,
        drop = FALSE
    ]

    if (any(dat[[count_col]] < 0)) {
        stop("Counts must be non-negative.")
    }


    # ==========================================================
    # 7. Filter rare OTUs
    # ==========================================================

    otu_stats <- dat |>
        dplyr::group_by(.data[[otu_col]]) |>
        dplyr::summarise(
            total_count =
                sum(.data[[count_col]], na.rm = TRUE),

            detected_replicates =
                sum(.data[[count_col]] > 0, na.rm = TRUE),

            .groups = "drop"
        )

    retained_otus <- otu_stats |>
        dplyr::filter(
            .data$total_count >= min_species_sum,
            .data$detected_replicates >= min_detection_replicates
        ) |>
        dplyr::pull(.data[[otu_col]])

    dat <- dat |>
        dplyr::filter(
            .data[[otu_col]] %in% retained_otus
        )

    if (nrow(dat) == 0L) {
        stop("No OTUs remain after filtering.")
    }


    # ==========================================================
    # 8. Construct hierarchical grouping factors
    # ==========================================================

    dat$.site <- factor(dat[[site_col]])
    dat$.otu <- factor(dat[[otu_col]])
    dat$.sample <- factor(dat[[sample_col]])

    # One latent occupancy state for each site x OTU.
    dat$.site_otu <- interaction(
        dat$.site,
        dat$.otu,
        drop = TRUE,
        lex.order = TRUE
    )

    # One latent capture state for each biological sample x OTU.
    dat$.sample_otu <- interaction(
        dat$.site,
        dat$.sample,
        dat$.otu,
        drop = TRUE,
        lex.order = TRUE
    )


    # ==========================================================
    # 9. Construct sample x OTU data
    # ==========================================================

    sample_df <- dat |>
        dplyr::group_by(.data$.sample_otu) |>
        dplyr::summarise(
            .site_otu = dplyr::first(.data$.site_otu),
            .otu = dplyr::first(.data$.otu),
            .sample = dplyr::first(.data$.sample),

            sample_positive =
                as.integer(any(.data[[count_col]] > 0)),

            .groups = "drop"
        )


    # ==========================================================
    # 10. Add capture covariates
    # ==========================================================

    cap_vars <- intersect(
        all.vars(capture_formula),
        names(dat)
    )

    if (length(cap_vars) > 0L) {

        cap_covars <- dat |>
            dplyr::select(
                .data$.sample_otu,
                dplyr::all_of(cap_vars)
            ) |>
            dplyr::group_by(.data$.sample_otu) |>
            dplyr::summarise(
                dplyr::across(
                    dplyr::all_of(cap_vars),
                    ~ dplyr::first(.x)
                ),
                .groups = "drop"
            )

        sample_df <- dplyr::left_join(
            sample_df,
            cap_covars,
            by = ".sample_otu"
        )
    }


    # ==========================================================
    # 11. Construct site x OTU data
    # ==========================================================

    site_df <- dat |>
        dplyr::group_by(.data$.site_otu) |>
        dplyr::summarise(
            .otu = dplyr::first(.data$.otu),

            site_positive =
                as.integer(any(.data[[count_col]] > 0)),

            .groups = "drop"
        )


    # ==========================================================
    # 12. Add occupancy covariates
    # ==========================================================

    occ_vars <- intersect(
        all.vars(occupancy_formula),
        names(dat)
    )

    if (length(occ_vars) > 0L) {

        occ_covars <- dat |>
            dplyr::select(
                .data$.site_otu,
                dplyr::all_of(occ_vars)
            ) |>
            dplyr::group_by(.data$.site_otu) |>
            dplyr::summarise(
                dplyr::across(
                    dplyr::all_of(occ_vars),
                    ~ dplyr::first(.x)
                ),
                .groups = "drop"
            )

        site_df <- dplyr::left_join(
            site_df,
            occ_covars,
            by = ".site_otu"
        )
    }


    # ==========================================================
    # 13. Construct fixed-effect design matrices
    # ==========================================================

    X_occ <- stats::model.matrix(
        occupancy_formula,
        data = site_df
    )

    X_cap <- stats::model.matrix(
        capture_formula,
        data = sample_df
    )

    X_abund <- stats::model.matrix(
        abundance_formula,
        data = dat
    )


    # ==========================================================
    # 14. Construct abundance offset
    #
    # The user supplies the exposure itself. We put it on
    # the log scale before passing it to the C++ model.
    # ==========================================================

    if (is.null(abundance_offset)) {

        offset_abund <- rep(0, nrow(dat))

    } else {

        if (!(abundance_offset %in% names(dat))) {
            stop(
                "abundance_offset '",
                abundance_offset,
                "' was not found."
            )
        }

        off <- as.numeric(
            dat[[abundance_offset]]
        )

        if (any(!is.finite(off) | off <= 0)) {
            stop(
                "abundance_offset must contain positive finite values."
            )
        }

        offset_abund <- log(off)
    }


    # ==========================================================
    # 15. Create mappings between hierarchy levels
    #
    # C++/TMB indices are zero based.
    # ==========================================================

    site_group_levels <- levels(dat$.site_otu)
    sample_group_levels <- levels(dat$.sample_otu)

    row_sample_group <- match(
        as.character(dat$.sample_otu),
        sample_group_levels
    ) - 1L

    sample_site_group <- match(
        as.character(sample_df$.site_otu),
        site_group_levels
    ) - 1L

    occ_otu <- as.integer(site_df$.otu) - 1L
    sample_otu <- as.integer(sample_df$.otu) - 1L
    row_otu <- as.integer(dat$.otu) - 1L
    row_sample_id <- as.integer(dat$.sample) - 1L

    # The sample x OTU abundance random effect uses the
    # same grouping as the sample x OTU capture state.
    row_sample_otu_re <- row_sample_group

    if (anyNA(row_sample_group)) {
        stop("Invalid row_sample_group mapping.")
    }

    if (anyNA(sample_site_group)) {
        stop("Invalid sample_site_group mapping.")
    }


    # ==========================================================
    # 16. Encode abundance distribution for C++
    # ==========================================================

    family_code <- switch(
        abundance_family,
        poisson = 0L,
        nbinom = 1L,
        zip = 2L,
        zinb = 3L
    )

    # OTU-specific zero inflation exists only for ZIP/ZINB.
    use_zi_otu <-
        isTRUE(random_zi_otu) &&
        abundance_family %in% c("zip", "zinb")


    # ==========================================================
    # 17. Construct TMB data list
    # ==========================================================

    data_tmb <- list(
        y = dat[[count_col]],

        X_occ = X_occ,
        X_cap = X_cap,
        X_abund = X_abund,

        offset_abund = offset_abund,

        occ_otu = as.integer(occ_otu),
        sample_otu = as.integer(sample_otu),

        row_sample_group =
            as.integer(row_sample_group),

        sample_site_group =
            as.integer(sample_site_group),

        row_otu =
            as.integer(row_otu),

        row_sample_id =
            as.integer(row_sample_id),

        row_sample_otu_re =
            as.integer(row_sample_otu_re),

        sample_positive =
            as.integer(sample_df$sample_positive),

        site_positive =
            as.integer(site_df$site_positive),

        n_site_groups =
            nrow(site_df),

        n_sample_groups =
            nrow(sample_df),

        family_code =
            family_code,

        use_occ_otu =
            as.integer(random_occ_otu),

        use_cap_otu =
            as.integer(random_capture_otu),

        use_abund_otu =
            as.integer(random_abund_otu),

        use_sample_re =
            as.integer(random_sample),

        use_sample_otu_re =
            as.integer(random_sample_otu),

        use_zi_otu =
            as.integer(use_zi_otu)
    )


    # ==========================================================
    # 18. Dimensions
    # ==========================================================

    n_otu <- nlevels(dat$.otu)
    n_sample <- nlevels(dat$.sample)
    n_sample_otu <- nrow(sample_df)


    # ==========================================================
    # 19. Generate sensible starting values
    # ==========================================================

    # Keep initial probabilities away from exactly 0 and 1.
    site_positive_rate <- mean(
        site_df$site_positive
    )

    site_positive_rate <- pmin(
        pmax(site_positive_rate, 0.05),
        0.95
    )

    sample_positive_rate <- mean(
        sample_df$sample_positive
    )

    sample_positive_rate <- pmin(
        pmax(sample_positive_rate, 0.05),
        0.95
    )

    positive_counts <- dat[[count_col]][
        dat[[count_col]] > 0
    ]

    mean_count_start <- if (length(positive_counts) > 0L) {
        mean(positive_counts)
    } else {
        1
    }

    beta_occ_start <- rep(
        0,
        ncol(X_occ)
    )

    beta_occ_start[1] <- stats::qlogis(
        site_positive_rate
    )

    beta_cap_start <- rep(
        0,
        ncol(X_cap)
    )

    beta_cap_start[1] <- stats::qlogis(
        sample_positive_rate
    )

    beta_abund_start <- rep(
        0,
        ncol(X_abund)
    )

    beta_abund_start[1] <- log(
        pmax(mean_count_start, 1e-4)
    )


    # ==========================================================
    # 20. Define TMB parameter list
    # ==========================================================

    parameters <- list(

        # Fixed effects
        beta_occ = beta_occ_start,
        beta_cap = beta_cap_start,
        beta_abund = beta_abund_start,

        # OTU random effects
        b_occ_otu = rep(0, n_otu),
        b_cap_otu = rep(0, n_otu),
        b_abund_otu = rep(0, n_otu),

        # Sample random effects
        b_sample = rep(0, n_sample),
        b_sample_otu = rep(0, n_sample_otu),

        # OTU-specific zero-inflation deviations
        b_zi_otu = rep(0, n_otu),

        # Random-effect log standard deviations
        log_sd_occ_otu = log(0.5),
        log_sd_cap_otu = log(0.5),
        log_sd_abund_otu = log(0.5),
        log_sd_sample = log(0.5),
        log_sd_sample_otu = log(0.5),
        log_sd_zi_otu = log(0.5),

        # Negative-binomial dispersion
        log_theta = log(10),

        # Population-level zero-inflation intercept
        zi_intercept = stats::qlogis(0.05)
    )


    # ==========================================================
    # 21. Parameter mapping
    #
    # Parameters that are not used by a particular model are
    # mapped to NA so that TMB does not estimate them.
    # ==========================================================

    map <- list()

    if (!random_occ_otu) {
        map$b_occ_otu <- factor(rep(NA, n_otu))
        map$log_sd_occ_otu <- factor(NA)
    }

    if (!random_capture_otu) {
        map$b_cap_otu <- factor(rep(NA, n_otu))
        map$log_sd_cap_otu <- factor(NA)
    }

    if (!random_abund_otu) {
        map$b_abund_otu <- factor(rep(NA, n_otu))
        map$log_sd_abund_otu <- factor(NA)
    }

    if (!random_sample) {
        map$b_sample <- factor(rep(NA, n_sample))
        map$log_sd_sample <- factor(NA)
    }

    if (!random_sample_otu) {
        map$b_sample_otu <- factor(rep(NA, n_sample_otu))
        map$log_sd_sample_otu <- factor(NA)
    }

    if (!use_zi_otu) {
        map$b_zi_otu <- factor(rep(NA, n_otu))
        map$log_sd_zi_otu <- factor(NA)
    }

    # Distribution-specific parameters.
    if (abundance_family == "poisson") {
        map$log_theta <- factor(NA)
        map$zi_intercept <- factor(NA)
    }

    if (abundance_family == "nbinom") {
        map$zi_intercept <- factor(NA)
    }

    if (abundance_family == "zip") {
        map$log_theta <- factor(NA)
    }


    # ==========================================================
    # 22. Tell TMB which parameters are random
    # ==========================================================

    random_effects <- character(0)

    if (random_occ_otu) {
        random_effects <- c(
            random_effects,
            "b_occ_otu"
        )
    }

    if (random_capture_otu) {
        random_effects <- c(
            random_effects,
            "b_cap_otu"
        )
    }

    if (random_abund_otu) {
        random_effects <- c(
            random_effects,
            "b_abund_otu"
        )
    }

    if (random_sample) {
        random_effects <- c(
            random_effects,
            "b_sample"
        )
    }

    if (random_sample_otu) {
        random_effects <- c(
            random_effects,
            "b_sample_otu"
        )
    }

    if (use_zi_otu) {
        random_effects <- c(
            random_effects,
            "b_zi_otu"
        )
    }


    # ==========================================================
    # 23. Print basic model information
    # ==========================================================

    if (verbose) {

        message("Sites x OTUs: ", nrow(site_df))
        message("Samples x OTUs: ", nrow(sample_df))
        message("Read observations: ", nrow(dat))
        message("OTUs: ", n_otu)
        message("Biological samples: ", n_sample)

        message(
            "Random effects: ",
            if (length(random_effects) == 0L) {
                "none"
            } else {
                paste(random_effects, collapse = ", ")
            }
        )

        if (
            random_zi_otu &&
            !(abundance_family %in% c("zip", "zinb"))
        ) {
            message(
                "Note: random_zi_otu = TRUE is ignored because ",
                "abundance_family = '",
                abundance_family,
                "' has no zero-inflation component."
            )
        }

        message(
            "sdreport report covariance: ",
            if (get_report_covariance) {
                "enabled"
            } else {
                "disabled (memory-efficient default)"
            }
        )
    }


    # ==========================================================
    # 24. Build the TMB objective
    # ==========================================================

    obj <- TMB::MakeADFun(
        data = data_tmb,
        parameters = parameters,

        random =
            if (length(random_effects) == 0L) {
                NULL
            } else {
                random_effects
            },

        map =
            if (length(map) == 0L) {
                NULL
            } else {
                map
            },

        DLL = DLL,
        silent = !verbose
    )


    # ==========================================================
    # 25. Check the initial objective
    # ==========================================================

    initial_nll <- obj$fn(obj$par)

    if (!is.finite(initial_nll)) {
        stop(
            "Initial joint negative log-likelihood is not finite."
        )
    }

    if (verbose) {
        message(
            "Initial negative log-likelihood: ",
            signif(initial_nll, 10)
        )

        message("Optimizing joint likelihood...")
    }


    # ==========================================================
    # 26. Helper function: evaluate the final outer gradient
    #
    # Calling fn() first ensures that the inner random-effect
    # mode corresponds to the supplied fixed parameters.
    # ==========================================================

    evaluate_gradient <- function(par) {

        invisible(
            obj$fn(par)
        )

        g <- tryCatch(
            obj$gr(par),
            error = function(e) {
                rep(NA_real_, length(par))
            }
        )

        if (length(g) == length(par)) {
            names(g) <- names(par)
        }

        gradient_finite <-
            length(g) > 0L &&
            all(is.finite(g))

        max_gradient <- if (gradient_finite) {
            max(abs(g))
        } else {
            Inf
        }

        list(
            gradient = g,
            finite = gradient_finite,
            max_abs_gradient = max_gradient,

            gradient_ok =
                gradient_finite &&
                max_gradient <= gradient_tol
        )
    }


    # ==========================================================
    # 27. Optimize with optional automatic restarts
    #
    # Restarting from the previous optimum can materially reduce
    # the outer gradient even when nlminb() has already returned
    # convergence code 0.
    # ==========================================================

    optimization_history <- list()

    current_start <- obj$par

    opt <- NULL

    previous_objective <- NA_real_
    previous_gradient <- NA_real_

    stalled <- FALSE

    max_passes <- max_restarts + 1L


    for (pass in seq_len(max_passes)) {

        if (verbose) {
            message("")
            message(
                "Optimization pass ",
                pass,
                " of ",
                max_passes,
                "..."
            )
        }

        opt <- stats::nlminb(
            start = current_start,
            objective = obj$fn,
            gradient = obj$gr,

            control = list(
                iter.max = iter_max,
                eval.max = eval_max,
                rel.tol = rel_tol
            )
        )

        grad_info <- evaluate_gradient(
            opt$par
        )

        optimizer_ok_pass <- identical(
            as.integer(opt$convergence),
            0L
        )

        # Save every optimization pass so advanced users can
        # inspect whether the restart sequence improved the fit.
        optimization_history[[pass]] <- list(
            pass = pass,
            objective = opt$objective,
            convergence_code = opt$convergence,
            optimizer_message = opt$message,
            optimizer_ok = optimizer_ok_pass,
            gradient = grad_info$gradient,
            max_abs_gradient =
                grad_info$max_abs_gradient,
            gradient_ok =
                grad_info$gradient_ok
        )

        if (verbose) {
            message(
                "  NLL: ",
                format(opt$objective, digits = 12)
            )

            message(
                "  Optimizer code: ",
                opt$convergence
            )

            message(
                "  Maximum absolute gradient: ",
                signif(
                    grad_info$max_abs_gradient,
                    8
                )
            )

            message(
                "  Gradient OK: ",
                grad_info$gradient_ok
            )
        }


        # ------------------------------------------------------
        # Both formal optimizer conditions have been satisfied.
        # ------------------------------------------------------

        if (
            optimizer_ok_pass &&
            grad_info$gradient_ok
        ) {
            if (verbose) {
                message(
                    "  Numerical optimizer criteria satisfied."
                )
            }

            break
        }


        # ------------------------------------------------------
        # Check whether repeated optimization has stalled.
        #
        # If both the objective and maximum gradient are
        # essentially unchanged, further restarts are unlikely
        # to help.
        # ------------------------------------------------------

        if (
            pass > 1L &&
            is.finite(previous_objective) &&
            is.finite(previous_gradient) &&
            is.finite(opt$objective) &&
            is.finite(grad_info$max_abs_gradient)
        ) {

            objective_change <- abs(
                opt$objective -
                    previous_objective
            )

            gradient_change <- abs(
                grad_info$max_abs_gradient -
                    previous_gradient
            )

            stalled <-
                objective_change <=
                    stall_objective_tol &&
                gradient_change <=
                    stall_gradient_tol

            if (stalled) {

                if (verbose) {
                    message(
                        "  Optimization has stalled."
                    )

                    message(
                        "  Objective change: ",
                        signif(
                            objective_change,
                            6
                        )
                    )

                    message(
                        "  Gradient change: ",
                        signif(
                            gradient_change,
                            6
                        )
                    )

                    message(
                        "  Additional restarts are unlikely ",
                        "to improve convergence."
                    )
                }

                break
            }
        }


        # No more permitted restarts.
        if (pass >= max_passes) {
            break
        }


        # Restart from the current solution.
        if (verbose) {
            message(
                "  Convergence criteria not satisfied."
            )

            message(
                "  Restarting optimization from current solution..."
            )
        }

        previous_objective <- opt$objective

        previous_gradient <-
            grad_info$max_abs_gradient

        current_start <- opt$par
    }


    # ==========================================================
    # 28. Evaluate final gradient
    # ==========================================================

    final_grad_info <- evaluate_gradient(
        opt$par
    )

    final_gradient <-
        final_grad_info$gradient

    gradient_finite <-
        final_grad_info$finite

    max_abs_gradient <-
        final_grad_info$max_abs_gradient

    gradient_ok <-
        final_grad_info$gradient_ok

    optimizer_ok <- identical(
        as.integer(opt$convergence),
        0L
    )


    # ==========================================================
    # 29. User-facing gradient classification
    #
    # IMPORTANT:
    # MARGINAL is descriptive only.
    #
    # Strict convergence still requires:
    #
    #     max_abs_gradient <= gradient_tol
    # ==========================================================

    gradient_status <- if (!gradient_finite) {

        "FAIL"

    } else if (max_abs_gradient <= gradient_tol) {

        "PASS"

    } else if (
        max_abs_gradient <=
            gradient_marginal_factor * gradient_tol
    ) {

        "MARGINAL"

    } else {

        "FAIL"
    }


    # ==========================================================
    # 30. Construct gradient diagnostic table
    # ==========================================================

    gradient_table <- data.frame(
        parameter = names(final_gradient),
        gradient = as.numeric(final_gradient),
        abs_gradient =
            abs(as.numeric(final_gradient)),
        stringsAsFactors = FALSE
    )

    gradient_table <- gradient_table[
        order(
            gradient_table$abs_gradient,
            decreasing = TRUE
        ),
        ,
        drop = FALSE
    ]

    n_show <- min(
        n_gradient_report,
        nrow(gradient_table)
    )

    largest_gradients <- if (n_show > 0L) {
        gradient_table[
            seq_len(n_show),
            ,
            drop = FALSE
        ]
    } else {
        gradient_table
    }


    # ==========================================================
    # 31. Calculate sdreport
    #
    # IMPORTANT:
    #
    # getReportCovariance is FALSE by default because the
    # covariance matrix of all ADREPORT quantities can become
    # extremely large for OTU-rich eDNA datasets.
    # ==========================================================

    sdr <- tryCatch(

        TMB::sdreport(
            obj,
            par.fixed = opt$par,

            getJointPrecision =
                isTRUE(get_joint_precision) &&
                length(random_effects) > 0L,

            getReportCovariance =
                isTRUE(get_report_covariance)
        ),

        error = function(e) {

            if (verbose) {
                message(
                    "sdreport failed: ",
                    conditionMessage(e)
                )
            }

            NULL
        }
    )

    sdreport_ok <- !is.null(sdr)


    # ==========================================================
    # 32. Positive-definite Hessian diagnostic
    #
    # pdHess is TMB's positive-definiteness diagnostic.
    #
    # This is PD, not merely PSD.
    # ==========================================================

    pd_hessian <- if (
        !is.null(sdr) &&
        !is.null(sdr$pdHess)
    ) {
        isTRUE(sdr$pdHess)
    } else {
        FALSE
    }


    # ==========================================================
    # 33. Obtain TMB report
    # ==========================================================

    report <- tryCatch(
        obj$report(),

        error = function(e) {
            warning(
                "TMB report failed: ",
                conditionMessage(e),
                call. = FALSE
            )

            NULL
        }
    )


    # ==========================================================
    # 34. Fixed-parameter summary
    # ==========================================================

    if (!is.null(sdr)) {

        fixed_matrix <- tryCatch(
            summary(sdr, "fixed"),
            error = function(e) NULL
        )

        if (!is.null(fixed_matrix)) {

            fixed_summary <- as.data.frame(
                fixed_matrix
            )

            fixed_summary$parameter <-
                rownames(fixed_summary)

            rownames(fixed_summary) <- NULL

        } else {

            fixed_summary <- data.frame()
        }

    } else {

        fixed_summary <- data.frame()
    }


    # ==========================================================
    # 35. Derived/ADREPORT parameter summary
    # ==========================================================

    if (!is.null(sdr)) {

        derived_matrix <- tryCatch(
            summary(sdr, "report"),
            error = function(e) NULL
        )

        if (!is.null(derived_matrix)) {

            derived_summary <- as.data.frame(
                derived_matrix
            )

            derived_summary$parameter <-
                rownames(derived_summary)

            rownames(derived_summary) <- NULL

        } else {

            derived_summary <- data.frame()
        }

    } else {

        derived_summary <- data.frame()
    }


    # ==========================================================
    # 36. Check that fixed-parameter SEs are finite
    #
    # This detects:
    #
    # NA
    # NaN
    # Inf
    # -Inf
    #
    # It does NOT imply that finite SEs are necessarily small
    # or scientifically informative.
    # ==========================================================

    finite_se <- if (
        nrow(fixed_summary) > 0L &&
        "Std. Error" %in% names(fixed_summary)
    ) {

        all(
            is.finite(
                fixed_summary[["Std. Error"]]
            )
        )

    } else {

        FALSE
    }


    # ==========================================================
    # 37. Heuristic diagnostic:
    #     near-zero variance components
    #
    # Example:
    #
    # log(SD) < -5
    #
    # corresponds to:
    #
    # SD < exp(-5) ~= 0.0067.
    #
    # This is a WARNING, not a convergence failure.
    # ==========================================================

    near_zero_sd <- data.frame()

    if (
        nrow(fixed_summary) > 0L &&
        all(
            c("Estimate", "parameter") %in%
                names(fixed_summary)
        )
    ) {

        log_sd_rows <- grepl(
            "^log_sd",
            fixed_summary$parameter
        )

        if (any(log_sd_rows)) {

            tmp_sd <- fixed_summary[
                log_sd_rows,
                ,
                drop = FALSE
            ]

            tmp_sd$sd_estimate <- exp(
                tmp_sd$Estimate
            )

            tmp_sd$near_zero <-
                tmp_sd$Estimate <
                log_sd_warning_threshold

            near_zero_sd <- tmp_sd[
                tmp_sd$near_zero,
                ,
                drop = FALSE
            ]
        }
    }

    variance_status <- if (
        nrow(near_zero_sd) == 0L
    ) {
        "PASS"
    } else {
        "WARNING"
    }


    # ==========================================================
    # 38. Heuristic diagnostic:
    #     large SE relative to estimate
    #
    # Do not calculate the ratio for estimates effectively
    # equal to zero because SE / |estimate| then becomes
    # inherently unstable and uninformative.
    # ==========================================================

    large_se_parameters <- data.frame()

    if (
        nrow(fixed_summary) > 0L &&
        all(
            c(
                "Estimate",
                "Std. Error",
                "parameter"
            ) %in% names(fixed_summary)
        )
    ) {

        tmp <- fixed_summary

        tmp$se_to_estimate_ratio <- NA_real_

        valid_ratio <-
            is.finite(tmp$Estimate) &
            is.finite(tmp[["Std. Error"]]) &
            abs(tmp$Estimate) >
                estimate_zero_tol

        tmp$se_to_estimate_ratio[
            valid_ratio
        ] <-
            tmp[["Std. Error"]][valid_ratio] /
            abs(tmp$Estimate[valid_ratio])

        tmp$large_relative_se <-
            !is.na(tmp$se_to_estimate_ratio) &
            tmp$se_to_estimate_ratio >
                se_estimate_ratio_threshold

        large_se_parameters <- tmp[
            tmp$large_relative_se,
            ,
            drop = FALSE
        ]
    }

    uncertainty_status <- if (
        nrow(large_se_parameters) == 0L
    ) {
        "PASS"
    } else {
        "WARNING"
    }


    # ==========================================================
    # 39. Individual formal diagnostic statuses
    # ==========================================================

    optimizer_status <- if (optimizer_ok) {
        "PASS"
    } else {
        "FAIL"
    }

    sdreport_status <- if (sdreport_ok) {
        "PASS"
    } else {
        "FAIL"
    }

    hessian_status <- if (pd_hessian) {
        "PASS"
    } else {
        "FAIL"
    }

    se_status <- if (finite_se) {
        "PASS"
    } else {
        "FAIL"
    }


    # ==========================================================
    # 40. Strict overall numerical convergence
    #
    # ALL five formal conditions must be satisfied.
    #
    # Heuristic parameter warnings are deliberately excluded.
    # ==========================================================

    converged <-
        optimizer_ok &&
        gradient_ok &&
        sdreport_ok &&
        pd_hessian &&
        finite_se


    # ==========================================================
    # 41. Overall user-facing status
    #
    # A fit may be described as MARGINAL if the only strict
    # failure is a gradient just above the requested tolerance.
    #
    # converged remains FALSE in that situation.
    # ==========================================================

    overall_status <- if (converged) {

        "PASS"

    } else if (
        optimizer_ok &&
        gradient_status == "MARGINAL" &&
        sdreport_ok &&
        pd_hessian &&
        finite_se
    ) {

        "MARGINAL"

    } else {

        "FAIL"
    }


    # ==========================================================
    # 42. Explain formal convergence failures
    # ==========================================================

    failure_reasons <- character(0)

    if (!optimizer_ok) {
        failure_reasons <- c(
            failure_reasons,
            paste0(
                "Optimizer returned convergence code ",
                opt$convergence,
                "."
            )
        )
    }

    if (!gradient_ok) {

        if (gradient_status == "MARGINAL") {

            failure_reasons <- c(
                failure_reasons,
                paste0(
                    "Maximum absolute gradient narrowly exceeds ",
                    "the specified tolerance."
                )
            )

        } else {

            failure_reasons <- c(
                failure_reasons,
                paste0(
                    "Maximum absolute gradient exceeds ",
                    "the specified tolerance."
                )
            )
        }
    }

    if (!sdreport_ok) {
        failure_reasons <- c(
            failure_reasons,
            "TMB::sdreport() was not successful."
        )
    }

    if (sdreport_ok && !pd_hessian) {
        failure_reasons <- c(
            failure_reasons,
            "The Hessian is not positive definite."
        )
    }

    if (sdreport_ok && !finite_se) {
        failure_reasons <- c(
            failure_reasons,
            paste0(
                "One or more fixed-parameter standard ",
                "errors are non-finite."
            )
        )
    }

    if (length(failure_reasons) == 0L) {
        failure_reasons <-
            "All numerical convergence criteria were satisfied."
    }


    # ==========================================================
    # 43. Construct heuristic warning messages
    # ==========================================================

    diagnostic_warnings <- character(0)

    if (nrow(near_zero_sd) > 0L) {

        diagnostic_warnings <- c(
            diagnostic_warnings,
            paste0(
                nrow(near_zero_sd),
                " random-effect log-SD parameter(s) are below ",
                log_sd_warning_threshold,
                ". This may indicate a variance component ",
                "close to zero but is not, by itself, a ",
                "numerical-convergence failure."
            )
        )
    }

    if (nrow(large_se_parameters) > 0L) {

        diagnostic_warnings <- c(
            diagnostic_warnings,
            paste0(
                nrow(large_se_parameters),
                " fixed parameter(s) have SE / |estimate| > ",
                se_estimate_ratio_threshold,
                ". This is a heuristic uncertainty flag and ",
                "does not, by itself, imply numerical failure."
            )
        )
    }


    # ==========================================================
    # 44. User-facing diagnostic table
    # ==========================================================

    diagnostic_table <- data.frame(

        Diagnostic = c(
            "Optimizer convergence",
            "Maximum absolute gradient",
            "sdreport",
            "Positive-definite Hessian",
            "Finite fixed-parameter SEs",
            "Near-zero variance components",
            "Large SE relative to estimate",
            "Overall numerical convergence"
        ),

        Value = c(
            paste0("code ", opt$convergence),

            format(
                max_abs_gradient,
                digits = 8,
                scientific = FALSE
            ),

            if (sdreport_ok) {
                "successful"
            } else {
                "failed"
            },

            as.character(pd_hessian),

            as.character(finite_se),

            as.character(
                nrow(near_zero_sd)
            ),

            as.character(
                nrow(large_se_parameters)
            ),

            as.character(converged)
        ),

        Status = c(
            optimizer_status,
            gradient_status,
            sdreport_status,
            hessian_status,
            se_status,
            variance_status,
            uncertainty_status,
            overall_status
        ),

        stringsAsFactors = FALSE
    )


    # ==========================================================
    # 45. Print final numerical diagnostics
    # ==========================================================

    if (verbose) {

        cat("\n")
        cat(
            "----- Numerical convergence diagnostics -----\n"
        )

        cat(
            "Optimizer convergence code: ",
            opt$convergence,
            "\n",
            sep = ""
        )

        cat(
            "Optimizer status:           ",
            optimizer_status,
            "\n\n",
            sep = ""
        )

        cat(
            "Maximum absolute gradient:  ",
            format(
                max_abs_gradient,
                digits = 8,
                scientific = FALSE
            ),
            "\n",
            sep = ""
        )

        cat(
            "Gradient tolerance:         ",
            format(
                gradient_tol,
                digits = 8,
                scientific = FALSE
            ),
            "\n",
            sep = ""
        )

        cat(
            "Gradient status:            ",
            gradient_status,
            "\n\n",
            sep = ""
        )

        cat(
            "sdreport status:            ",
            sdreport_status,
            "\n",
            sep = ""
        )

        cat(
            "Positive-definite Hessian:  ",
            pd_hessian,
            " [",
            hessian_status,
            "]\n",
            sep = ""
        )

        cat(
            "Finite fixed-parameter SEs: ",
            finite_se,
            " [",
            se_status,
            "]\n\n",
            sep = ""
        )

        cat(
            "Overall numerical status:   ",
            overall_status,
            "\n",
            sep = ""
        )

        cat(
            "Strict convergence:         ",
            converged,
            "\n\n",
            sep = ""
        )

        cat("Reason:\n")

        for (reason in failure_reasons) {
            cat(
                "  ",
                reason,
                "\n",
                sep = ""
            )
        }


        # Print the parameters responsible for the largest
        # remaining gradient when the gradient criterion fails.
        if (
            !gradient_ok &&
            nrow(largest_gradients) > 0L
        ) {

            cat("\n")
            cat(
                "Largest gradient parameters:\n"
            )

            for (
                i in seq_len(
                    nrow(largest_gradients)
                )
            ) {

                cat(
                    sprintf(
                        "  %-20s = %.8g\n",
                        largest_gradients$parameter[i],
                        largest_gradients$abs_gradient[i]
                    )
                )
            }
        }


        cat("\n")
        cat(
            "----- Parameter diagnostics -----\n"
        )

        cat(
            "Near-zero variance components: ",
            variance_status,
            " (",
            nrow(near_zero_sd),
            " flagged)\n",
            sep = ""
        )

        cat(
            "Large relative SEs:             ",
            uncertainty_status,
            " (",
            nrow(large_se_parameters),
            " flagged)\n",
            sep = ""
        )

        if (length(diagnostic_warnings) > 0L) {

            cat("\n")
            cat("Heuristic warnings:\n")

            for (w in diagnostic_warnings) {
                cat(
                    "  ",
                    w,
                    "\n",
                    sep = ""
                )
            }
        }

        cat(
            "---------------------------------------------\n"
        )
    }


    # ==========================================================
    # 46. Issue an R warning when strict convergence fails
    # ==========================================================

    if (!converged) {

        # Special warning when every diagnostic except the
        # strict gradient criterion passes.
        if (
            optimizer_ok &&
            gradient_status == "MARGINAL" &&
            sdreport_ok &&
            pd_hessian &&
            finite_se
        ) {

            warning(
                sprintf(
                    paste0(
                        "The fit is numerically marginal: optimizer, ",
                        "Hessian, sdreport and standard-error diagnostics ",
                        "are satisfactory, but the maximum absolute ",
                        "gradient (%.6g) exceeds gradient_tol (%.6g). ",
                        "Strict convergence is FALSE."
                    ),
                    max_abs_gradient,
                    gradient_tol
                ),
                call. = FALSE
            )

        } else {

            warning(
                paste0(
                    "Fit did not satisfy all numerical convergence ",
                    "criteria. ",
                    paste(
                        failure_reasons,
                        collapse = " "
                    )
                ),
                call. = FALSE
            )
        }
    }


    # ==========================================================
    # 47. Return fitted model
    # ==========================================================

    out <- list(

        # Core TMB objects
        fit = opt,
        tmb_object = obj,
        sdreport = sdr,
        report = report,

        # Parameter summaries
        fixed_effects = fixed_summary,
        derived = derived_summary,

        # Processed data
        site_data = site_df,
        sample_data = sample_df,
        long_df = dat,

        # OTU filtering information
        otu_stats = otu_stats,
        retained_otus = retained_otus,

        # TMB construction information
        data_tmb = data_tmb,
        parameters_start = parameters,
        parameter_map = map,
        random_effects = random_effects,

        # Model specification
        formulas = list(
            occupancy = occupancy_formula,
            capture = capture_formula,
            abundance = abundance_formula
        ),

        abundance_family =
            abundance_family,


        # ------------------------------------------------------
        # Formal numerical-convergence diagnostics
        # ------------------------------------------------------
        convergence = list(

            converged = converged,
            overall_status = overall_status,

            optimizer_ok = optimizer_ok,
            optimizer_status = optimizer_status,

            optimizer_code =
                opt$convergence,

            optimizer_message =
                opt$message,

            # Backward-compatible aliases
            code = opt$convergence,
            message = opt$message,

            objective =
                opt$objective,

            gradient =
                final_gradient,

            max_abs_gradient =
                max_abs_gradient,

            gradient_tolerance =
                gradient_tol,

            gradient_marginal_factor =
                gradient_marginal_factor,

            gradient_ok =
                gradient_ok,

            gradient_status =
                gradient_status,

            gradient_finite =
                gradient_finite,

            gradient_table =
                gradient_table,

            largest_gradients =
                largest_gradients,

            sdreport_ok =
                sdreport_ok,

            sdreport_status =
                sdreport_status,

            pd_hessian =
                pd_hessian,

            hessian_status =
                hessian_status,

            finite_standard_errors =
                finite_se,

            finite_fixed_parameter_standard_errors =
                finite_se,

            se_status =
                se_status,

            failure_reasons =
                failure_reasons,

            optimization_history =
                optimization_history,

            n_optimization_passes =
                length(
                    optimization_history
                ),

            max_restarts =
                max_restarts,

            stalled =
                stalled
        ),


        # ------------------------------------------------------
        # Heuristic parameter diagnostics
        #
        # These do not define formal convergence.
        # ------------------------------------------------------
        diagnostics = list(

            table =
                diagnostic_table,

            near_zero_sd =
                near_zero_sd,

            near_zero_sd_status =
                variance_status,

            log_sd_warning_threshold =
                log_sd_warning_threshold,

            large_se_parameters =
                large_se_parameters,

            large_se_status =
                uncertainty_status,

            se_estimate_ratio_threshold =
                se_estimate_ratio_threshold,

            warnings =
                diagnostic_warnings
        ),


        # ------------------------------------------------------
        # Record expensive sdreport settings for reproducibility.
        # ------------------------------------------------------
        sdreport_settings = list(

            getReportCovariance =
                get_report_covariance,

            getJointPrecision =
                get_joint_precision
        ),


        # ------------------------------------------------------
        # Record zero-inflation structure.
        # ------------------------------------------------------
        zero_inflation_structure = list(

            enabled =
                abundance_family %in%
                c("zip", "zinb"),

            random_zi_otu =
                use_zi_otu
        ),


        note = paste(
            "Joint observed-data likelihood.",
            "Latent occupancy Z and biological-sample capture A",
            "are analytically marginalized.",
            "Enabled Gaussian random effects are integrated using",
            "the TMB Laplace approximation.",
            "Strict numerical convergence requires optimizer code 0,",
            "a finite outer gradient below the specified tolerance,",
            "successful sdreport, a positive-definite Hessian, and",
            "finite fixed-parameter standard errors.",
            "PASS/MARGINAL/FAIL labels aid interpretation but do not",
            "alter the strict convergence definition.",
            "Near-zero variance components and large SE-to-estimate",
            "ratios are heuristic warnings only.",
            "getReportCovariance is FALSE by default to reduce",
            "memory requirements for large eDNA datasets."
        )
    )


    # Give the object a class so that summary(), print(), etc.
    # can later be implemented as S3 methods without changing
    # the underlying return structure.
    class(out) <- c(
        "eDNAModel_joint",
        "list"
    )

    return(out)
}
