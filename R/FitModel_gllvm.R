#' Fit a Multispecies 3-Level Occupancy–Detection–Abundance Model using GLLVM and GLMM
#'
#' This function fits a hierarchical three-level model for eDNA / microbiome data
#' stored as a \code{phyloseq} object. The model separates:
#'
#' \itemize{
#'   \item \strong{Site level (Z)}: true presence/absence of each taxon
#'   \item \strong{Biological sample level (A)}: capture / detection process
#'   \item \strong{PCR replicate level (Y)}: observed counts
#' }
#'
#' Occupancy is modeled at the site-by-taxon level using a binomial
#' Generalized Linear Latent Variable Model (GLLVM), while capture and
#' abundance are modeled using GLMMs via \code{glmmTMB}.
#'
#' @param phyloseq A \code{phyloseq} object containing count data and sample metadata.
#' @param site_col Character string identifying sites (occupancy level).
#' @param sample_col Character string identifying biological samples (capture level).
#' @param replicate_col Character or \code{NULL} identifying PCR replicates (abundance level).
#' @param otu_col Character string identifying taxa/species column.
#' @param count_col Character string identifying the response/count variable.
#' @param abundance_rhs Right-hand side of the abundance model, e.g.
#'   \code{(1 | taxon)} or \code{offset(log(total_reads)) + (1 | taxon)}.
#' @param capture_formula Formula for capture model (response must be \code{a_sim}).
#' @param occupancy_covars Optional character vector of site-level covariates.
#' @param min_species_sum Minimum total count required to retain a taxon.
#' @param min_detection_replicates Minimum number of detections required per taxon.
#' @param abundance_threshold Threshold defining detection (default = 1).
#' @param n_iter Number of EM-like iterations.
#' @param burn_in Number of initial iterations discarded.
#' @param abundance_family One of \code{"poisson"}, \code{"nbinom"}, \code{"zip"}, \code{"zinb"}.
#' @param num_lv_c Number of latent variables in the GLLVM occupancy model.
#' @param verbose Logical; print progress during fitting.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{Posterior summaries at the site-by-taxon level for:
#'     occupancy (\code{psi_*}), capture (\code{capture_*}),
#'     abundance (\code{lambda_*}), and detection (\code{p_detect_*}).}
#'   \item{capture}{Sample-level capture summaries.}
#'   \item{capture_site}{Capture summaries aggregated to the site level.}
#'   \item{psi_list, capture_list, lambda_list, p_detect_list}{Per-iteration linear predictors.}
#'   \item{occupancy_models}{List of fitted GLLVM occupancy models.}
#'   \item{capture_models}{List of fitted GLMM capture models.}
#'   \item{abundance_models}{List of fitted GLMM abundance models.}
#'   \item{reduced_data, sample_data, long_df}{Processed hierarchical datasets.}
#'   \item{lv_sites, lv_species}{Latent variable coordinates from the GLLVM.}
#'   \item{mean_lv_sites, mean_lv_species}{Posterior mean latent variables.}
#'   \item{filter_summary}{Filtering information for retained taxa.}
#'   \item{diagnostic_AIC}{Per-iteration model diagnostics.}
#' }
#'
#' @details
#' The model defines a three-level hierarchical structure:
#'
#' \deqn{
#' Z_{i,m} \sim \mathrm{Bernoulli}(\psi_{i,m})
#' }
#'
#' \deqn{
#' A_{i,j,m} \mid Z_{i,m} \sim \mathrm{Bernoulli}(p_{i,j,m} \cdot Z_{i,m})
#' }
#'
#' \deqn{
#' Y_{i,j,k,m} \mid A_{i,j,m} \sim \mathrm{Count}(\lambda_{i,j,k,m} \cdot A_{i,j,m})
#' }
#'
#' where:
#' \itemize{
#'   \item \(i\): site,
#'   \item \(j\): biological sample,
#'   \item \(k\): PCR replicate,
#'   \item \(m\): taxon.
#' }
#'
#' Observed occupancy is defined as:
#'
#' \deqn{
#' z^{obs}_{i,m} = I\left(\max_{j,k} Y_{i,j,k,m} > c\right)
#' }
#'
#' where \(c =\) \code{abundance_threshold}.
#'
#' The occupancy model is:
#'
#' \deqn{
#' \mathrm{logit}(\psi_{i,m}) = X_i \beta_m + \text{latent variables}
#' }
#'
#' The capture model is:
#'
#' \deqn{
#' \mathrm{logit}(p_{i,j,m}) = W_{i,j} \gamma_m
#' }
#'
#' The abundance model defines \eqn{\lambda}, which induces a detection probability:
#'
#' \deqn{
#' p_{\mathrm{detect}} = 1 - P(Y = 0 \mid A = 1)
#' }
#'
#' Detection is computed internally using the complementary log-log link:
#'
#' \deqn{
#' \eta = \log\left(-\log(P(Y = 0))\right)
#' }
#'
#' and transformed back to the natural scale.
#'
#' \strong{Important:}
#' \itemize{
#'   \item \code{p_detect} is a deterministic function of the abundance model
#'   \item No standard error is directly estimated for \code{p_detect}
#' }
#'
#' @section EM-like algorithm:
#' The algorithm follows a Monte Carlo EM (data augmentation) scheme:
#' \enumerate{
#'   \item Fit GLLVM occupancy model for \(Z\)
#'   \item Fit GLMM capture model for \(A \mid Z = 1\)
#'   \item Fit GLMM abundance model for \(Y \mid A = 1\)
#'   \item Update latent \(A\) via posterior sampling
#'   \item Update latent \(Z\) via posterior sampling
#' }
#'
#' @section Model features:
#' \itemize{
#'   \item Fully general (no hard-coded taxon or response names)
#'   \item Three-level hierarchy (Z → A → Y)
#'   \item Latent ecological structure via GLLVM
#'   \item Flexible abundance families (Poisson, NB, ZIP, ZINB)
#'   \item Supports offsets (e.g. \code{offset(log(total_reads))})
#' }
#'
#' @section Caveats:
#' \itemize{
#'   \item Approximate inference (not full joint likelihood)
#'   \item Detection depends on abundance model specification
#'   \item Large \eqn{\lambda} implies near-certain detection
#' }
#'
#' @examples
#' \dontrun{
#' fit <- FitModel_gllvm(
#'   phyloseq = ps,
#'   site_col = "site_month",
#'   sample_col = "Name",
#'   replicate_col = "Replicate",
#'   otu_col = "Taxon",
#'   count_col = "count",
#'   abundance_rhs = y ~ (1 | Taxon),
#'   capture_formula = a_sim ~ 1 + (1 | Taxon),
#'   abundance_family = "nbinom",
#'   n_iter = 10,
#'   burn_in = 2
#' )
#'
#' head(fit$summary)
#' }
#'
#' @import dplyr
#' @import tidyr
#' @importFrom gllvm gllvm
#' @importFrom glmmTMB glmmTMB nbinom2
#' @importFrom reshape2 melt acast
#' @export
FitModel_gllvm <- function(
  phyloseq,
  site_col,
  sample_col = "Name",
  replicate_col = NULL,
  otu_col = "OTU",
  count_col = "y",
  abundance_rhs,
  capture_formula = NULL,
  occupancy_covars = NULL,
  min_species_sum = 50,
  min_detection_replicates = 1,
  abundance_threshold = 1,
  n_iter = 50,
  burn_in = 10,
  abundance_family = c("poisson", "nbinom", "zip", "zinb"),
  num_lv_c = 2,
  verbose = TRUE
) {

  # ============================================================
  # 0. Basic argument checks
  # ============================================================

  abundance_family <- match.arg(abundance_family)

  if (!is.numeric(n_iter) || length(n_iter) != 1L ||
      n_iter < 1 || n_iter != as.integer(n_iter)) {
    stop("n_iter must be a positive integer.")
  }

  if (!is.numeric(burn_in) || length(burn_in) != 1L ||
      burn_in < 0 || burn_in != as.integer(burn_in)) {
    stop("burn_in must be a non-negative integer.")
  }

  if (burn_in >= n_iter) {
    stop("burn_in must be smaller than n_iter.")
  }

  if (!is.numeric(num_lv_c) || length(num_lv_c) != 1L ||
      num_lv_c < 0 || num_lv_c != as.integer(num_lv_c)) {
    stop("num_lv_c must be a non-negative integer.")
  }

  # ============================================================
  # Helper functions
  # ============================================================

  bt <- function(x) {
    paste0("`", x, "`")
  }

  # Force structural identifiers to factors
  to_factor_cols <- function(df, cols) {

    cols <- unique(stats::na.omit(cols))
    cols <- intersect(cols, names(df))

    for (col in cols) {
      df[[col]] <- as.factor(as.character(df[[col]]))
    }

    df
  }

  # Convert character predictors to factors, while preserving
  # genuinely numeric covariates as numeric
  character_to_factor <- function(df, cols) {

    cols <- unique(stats::na.omit(cols))
    cols <- intersect(cols, names(df))

    for (col in cols) {

      if (is.character(df[[col]])) {
        df[[col]] <- factor(df[[col]])
      }
    }

    df
  }

  get_formula_vars <- function(formula, response) {
    setdiff(all.vars(formula), response)
  }

  # ============================================================
  # 1. SAFE abundance formula handling
  #
  # Supports BOTH:
  #
  # abundance_rhs =
  #   (1 | OTU) + (1 | Samplingmonth / OTU)
  #
  # AND:
  #
  # abundance_rhs =
  #   y ~ (1 | OTU) + (1 | Samplingmonth / OTU)
  #
  # IMPORTANT:
  # substitute() MUST be called before abundance_rhs is evaluated.
  # ============================================================

  abundance_rhs_expr <- substitute(abundance_rhs)

  abundance_rhs_txt <- paste(
    deparse(abundance_rhs_expr),
    collapse = " "
  )

  is_full_formula <-
    is.call(abundance_rhs_expr) &&
    identical(abundance_rhs_expr[[1L]], as.name("~"))

  if (is_full_formula) {

    abundance_formula <- stats::as.formula(
      abundance_rhs_txt,
      env = parent.frame()
    )

  } else {

    abundance_formula <- stats::as.formula(
      paste(
        bt(count_col),
        "~",
        abundance_rhs_txt
      ),
      env = parent.frame()
    )
  }

  # ============================================================
  # 2. Default capture formula
  # ============================================================

  if (is.null(capture_formula)) {

    capture_formula <- stats::as.formula(
      paste(
        "a_sim ~ 1 + (1 |",
        bt(otu_col),
        ")"
      ),
      env = parent.frame()
    )
  }

  if (!inherits(capture_formula, "formula")) {
    stop("capture_formula must be a formula.")
  }

  # ============================================================
  # 3. Prepare long data
  # ============================================================

  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col,
    nested_cols = unique(
      stats::na.omit(
        c(sample_col, replicate_col)
      )
    )
  )

  long_df <- prep$long_df

  # ============================================================
  # 4. Required-column checks
  # ============================================================

  needed_cols <- unique(
    c(
      site_col,
      sample_col,
      otu_col,
      count_col
    )
  )

  missing_cols <- setdiff(
    needed_cols,
    names(long_df)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in long_df: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.null(replicate_col) &&
      !(replicate_col %in% names(long_df))) {

    stop(
      "replicate_col '",
      replicate_col,
      "' not found in long_df."
    )
  }

  if (!is.null(occupancy_covars) &&
      length(occupancy_covars) > 0) {

    missing_occ <- setdiff(
      occupancy_covars,
      names(long_df)
    )

    if (length(missing_occ) > 0) {
      stop(
        "Missing occupancy covariates in long_df: ",
        paste(missing_occ, collapse = ", ")
      )
    }
  }

  # ============================================================
  # 5. Type handling
  # ============================================================

  long_df[[site_col]]   <- as.character(long_df[[site_col]])
  long_df[[sample_col]] <- as.character(long_df[[sample_col]])
  long_df[[otu_col]]    <- as.character(long_df[[otu_col]])
  long_df[[count_col]]  <- as.numeric(long_df[[count_col]])

  if (!is.null(replicate_col)) {
    long_df[[replicate_col]] <-
      as.character(long_df[[replicate_col]])
  }

  long_df <- long_df[
    !is.na(long_df[[count_col]]),
    ,
    drop = FALSE
  ]

  if (nrow(long_df) == 0) {
    stop("No non-missing count observations remain.")
  }

  # ============================================================
  # 6. Compute sample-level sequencing depth
  #
  # Use site + sample rather than sample alone so the function
  # remains valid even when sample labels repeat among sites.
  # ============================================================

  sample_id_keys <- unique(
    c(site_col, sample_col)
  )

  total_reads_df <- long_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(sample_id_keys)
      )
    ) |>
    dplyr::summarise(
      total_reads = sum(
        .data[[count_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  long_df <- dplyr::left_join(
    long_df,
    total_reads_df,
    by = sample_id_keys
  )

  # ============================================================
  # 7. Offset handling
  # ============================================================

  uses_offset <- any(
    grepl(
      "offset\\s*\\(",
      deparse(abundance_formula)
    )
  )

  if (uses_offset) {

    zero_samples <- total_reads_df |>
      dplyr::filter(
        is.na(.data$total_reads) |
        .data$total_reads <= 0
      )

    if (nrow(zero_samples) > 0) {

      if (verbose) {
        message(
          "Offset detected -> removing ",
          nrow(zero_samples),
          " zero-read sample(s)."
        )
      }

      zero_samples$.drop_sample <- TRUE

      long_df <- long_df |>
        dplyr::left_join(
          zero_samples |>
            dplyr::select(
              dplyr::all_of(sample_id_keys),
              .data$.drop_sample
            ),
          by = sample_id_keys
        ) |>
        dplyr::filter(
          is.na(.data$.drop_sample)
        ) |>
        dplyr::select(
          -dplyr::any_of(".drop_sample")
        )
    }

    if (any(
      long_df$total_reads <= 0,
      na.rm = TRUE
    )) {

      stop(
        "total_reads must be positive when using ",
        "offset(log(total_reads))."
      )
    }

  } else {

    if (verbose) {
      message("No abundance offset detected.")
    }
  }

  # ============================================================
  # 8. OTU filtering
  # ============================================================

  otu_stats <- long_df |>
    dplyr::group_by(
      .data[[otu_col]]
    ) |>
    dplyr::summarise(
      total_count = sum(
        .data[[count_col]],
        na.rm = TRUE
      ),
      detected_replicates = sum(
        .data[[count_col]] > abundance_threshold,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  keep_otus <- otu_stats |>
    dplyr::filter(
      .data$total_count >= min_species_sum,
      .data$detected_replicates >=
        min_detection_replicates
    ) |>
    dplyr::pull(
      .data[[otu_col]]
    )

  long_df <- long_df |>
    dplyr::filter(
      .data[[otu_col]] %in% keep_otus
    )

  if (nrow(long_df) == 0) {
    stop("No OTUs remain after filtering.")
  }

  if (length(keep_otus) < 2) {
    stop(
      "At least two OTUs are required for the GLLVM occupancy model."
    )
  }

  if (verbose) {

    message(
      "OTUs before filtering: ",
      dplyr::n_distinct(
        otu_stats[[otu_col]]
      )
    )

    message(
      "OTUs after filtering: ",
      dplyr::n_distinct(
        long_df[[otu_col]]
      )
    )
  }

  # ============================================================
  # 9. Factor / covariate handling
  # ============================================================

  id_factor_cols <- unique(
    stats::na.omit(
      c(
        site_col,
        sample_col,
        replicate_col,
        otu_col
      )
    )
  )

  formula_predictors <- unique(
    c(
      occupancy_covars,
      all.vars(capture_formula),
      all.vars(abundance_formula)
    )
  )

  formula_predictors <- setdiff(
    formula_predictors,
    c(
      count_col,
      "total_reads",
      "z_sim",
      "a_sim"
    )
  )

  # IDs always factors
  long_df <- to_factor_cols(
    long_df,
    id_factor_cols
  )

  # Character covariates -> factors,
  # numeric covariates stay numeric
  long_df <- character_to_factor(
    long_df,
    formula_predictors
  )

  # ============================================================
  # 10. Set reference OTU
  # ============================================================

  otu_abundances <- long_df |>
    dplyr::group_by(
      .data[[otu_col]]
    ) |>
    dplyr::summarise(
      total_count = sum(
        .data[[count_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  top_otu <- otu_abundances |>
    dplyr::arrange(
      dplyr::desc(.data$total_count)
    ) |>
    dplyr::slice(1) |>
    dplyr::pull(
      .data[[otu_col]]
    )

  long_df[[otu_col]] <-
    stats::relevel(
      factor(long_df[[otu_col]]),
      ref = as.character(top_otu)
    )

  # ============================================================
  # 11. Define hierarchical keys
  # ============================================================

  site_keys <- unique(
    c(site_col, otu_col)
  )

  sample_keys <- unique(
    c(site_col, sample_col, otu_col)
  )

  pcr_keys <- unique(
    c(
      site_col,
      sample_col,
      replicate_col,
      otu_col
    )
  )

  pcr_keys <- pcr_keys[
    !is.na(pcr_keys)
  ]

  # ============================================================
  # 12. Build site-level latent occupancy data Z
  # ============================================================

  site_keep_vars <- if (
    !is.null(occupancy_covars)
  ) {

    occupancy_covars

  } else {

    character(0)
  }

  reduced_data <- long_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(site_keys)
      )
    ) |>
    dplyr::summarise(
      z_obs = as.integer(
        any(
          .data[[count_col]] >
            abundance_threshold
        )
      ),
      dplyr::across(
        dplyr::all_of(site_keep_vars),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      z_sim = .data$z_obs
    )

  reduced_data <- to_factor_cols(
    reduced_data,
    intersect(
      id_factor_cols,
      names(reduced_data)
    )
  )

  reduced_data <- character_to_factor(
    reduced_data,
    occupancy_covars
  )

  # ============================================================
  # 13. Build biological-sample latent capture data A
  # ============================================================

  cap_vars <- get_formula_vars(
    capture_formula,
    "a_sim"
  )

  sample_keep_vars <- setdiff(
    cap_vars,
    c(
      sample_keys,
      "a_sim"
    )
  )

  sample_keep_vars <- intersect(
    sample_keep_vars,
    names(long_df)
  )

  sample_data <- long_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(sample_keys)
      )
    ) |>
    dplyr::summarise(
      a_obs = as.integer(
        any(
          .data[[count_col]] >
            abundance_threshold
        )
      ),
      dplyr::across(
        dplyr::all_of(sample_keep_vars),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      a_sim = .data$a_obs
    )

  sample_data <- to_factor_cols(
    sample_data,
    intersect(
      id_factor_cols,
      names(sample_data)
    )
  )

  sample_data <- character_to_factor(
    sample_data,
    cap_vars
  )

  # ============================================================
  # 14. Formula validation
  # ============================================================

  validate_formula <- function(
    formula,
    data,
    expected_response,
    model_name
  ) {

    response <- all.vars(formula)[1]

    if (!identical(
      response,
      expected_response
    )) {

      stop(
        model_name,
        " formula response must be '",
        expected_response,
        "', but got '",
        response,
        "'."
      )
    }

    missing_vars <- setdiff(
      all.vars(formula),
      names(data)
    )

    if (length(missing_vars) > 0) {

      stop(
        "Missing variables in ",
        model_name,
        " formula: ",
        paste(
          missing_vars,
          collapse = ", "
        )
      )
    }

    invisible(TRUE)
  }

  validate_formula(
    capture_formula,
    sample_data,
    "a_sim",
    "capture"
  )

  validate_formula(
    abundance_formula,
    long_df,
    count_col,
    "abundance"
  )

  # ============================================================
  # 15. Abundance family setup
  # ============================================================

  abundance_glmm_family <- switch(
    abundance_family,

    poisson =
      stats::poisson(),

    nbinom =
      glmmTMB::nbinom2(),

    zip =
      stats::poisson(),

    zinb =
      glmmTMB::nbinom2()
  )

  zi_formula <- if (
    abundance_family %in% c(
      "zip",
      "zinb"
    )
  ) {

    ~1

  } else {

    ~0
  }

  # ============================================================
  # 16. Zero probability from abundance model
  # ============================================================

  get_p0_abundance <- function(
    fit,
    newdata,
    family_name
  ) {

    # For zero-inflated models we need the conditional count mean,
    # not the zero-inflation-adjusted response expectation.
    pred_type <- if (
      family_name %in% c(
        "zip",
        "zinb"
      )
    ) {

      "conditional"

    } else {

      "response"
    }

    mu <- as.numeric(
      stats::predict(
        fit,
        type = pred_type,
        newdata = newdata,
        allow.new.levels = TRUE
      )
    )

    mu <- pmax(
      mu,
      1e-12
    )

    zi <- rep(
      0,
      length(mu)
    )

    if (
      family_name %in% c(
        "zip",
        "zinb"
      )
    ) {

      zi <- tryCatch(
        as.numeric(
          stats::predict(
            fit,
            type = "zprob",
            newdata = newdata,
            allow.new.levels = TRUE
          )
        ),
        error = function(e) {
          rep(
            0,
            length(mu)
          )
        }
      )

      zi <- pmin(
        pmax(
          zi,
          0
        ),
        1
      )
    }

    p0_cond <- switch(
      family_name,

      poisson =
        stats::dpois(
          0,
          lambda = mu
        ),

      zip =
        stats::dpois(
          0,
          lambda = mu
        ),

      nbinom = {

        theta <- tryCatch(
          as.numeric(
            stats::sigma(fit)
          ),
          error = function(e) NA_real_
        )

        if (
          !is.finite(theta) ||
          theta <= 0
        ) {

          stop(
            "Invalid negative-binomial dispersion parameter."
          )
        }

        stats::dnbinom(
          0,
          mu = mu,
          size = theta
        )
      },

      zinb = {

        theta <- tryCatch(
          as.numeric(
            stats::sigma(fit)
          ),
          error = function(e) NA_real_
        )

        if (
          !is.finite(theta) ||
          theta <= 0
        ) {

          stop(
            "Invalid negative-binomial dispersion parameter."
          )
        }

        stats::dnbinom(
          0,
          mu = mu,
          size = theta
        )
      }
    )

    p0 <- zi +
      (1 - zi) *
      p0_cond

    pmin(
      pmax(
        p0,
        1e-12
      ),
      1
    )
  }

  # ============================================================
  # 17. Precompute site x OTU template
  # ============================================================

  z_template <- reshape2::acast(
    reduced_data,
    stats::as.formula(
      paste(
        bt(site_col),
        "~",
        bt(otu_col)
      )
    ),
    value.var = "z_sim",
    fill = 0
  )

  z_sites <- rownames(z_template)
  z_otus  <- colnames(z_template)

  if (length(z_sites) < 2) {
    stop(
      "At least two sites are required for the GLLVM occupancy model."
    )
  }

  # ============================================================
  # 18. Occupancy design matrix
  # ============================================================

  X_cov <- NULL

  if (
    !is.null(occupancy_covars) &&
    length(occupancy_covars) > 0
  ) {

    cov_cols <- unique(
      c(
        site_col,
        occupancy_covars
      )
    )

    cov_df <- reduced_data |>
      dplyr::select(
        dplyr::all_of(cov_cols)
      ) |>
      dplyr::distinct()

    cov_df[[site_col]] <-
      as.character(
        cov_df[[site_col]]
      )

    # Covariates must be constant within site
    cov_count <- cov_df |>
      dplyr::count(
        .data[[site_col]],
        name = ".n_cov_rows"
      )

    bad_sites <- cov_count |>
      dplyr::filter(
        .data$.n_cov_rows > 1
      )

    if (nrow(bad_sites) > 0) {

      stop(
        "occupancy_covars must have one unique combination per site. ",
        "Covariates vary within ",
        nrow(bad_sites),
        " site(s)."
      )
    }

    cov_df <- character_to_factor(
      cov_df,
      occupancy_covars
    )

    cov_df <- cov_df[
      match(
        z_sites,
        cov_df[[site_col]]
      ),
      ,
      drop = FALSE
    ]

    if (anyNA(cov_df[[site_col]])) {

      stop(
        "Could not align occupancy covariates with all GLLVM sites."
      )
    }

    X_full <- stats::model.matrix(
      ~ .,
      data = cov_df[
        ,
        occupancy_covars,
        drop = FALSE
      ]
    )

    # Remove intercept column only if present
    if (
      ncol(X_full) > 0 &&
      "(Intercept)" %in% colnames(X_full)
    ) {

      X_cov <- X_full[
        ,
        colnames(X_full) != "(Intercept)",
        drop = FALSE
      ]

    } else {

      X_cov <- X_full
    }

    if (ncol(X_cov) == 0) {
      X_cov <- NULL
    }
  }

  # ============================================================
  # 19. Storage
  # ============================================================

  psi_list <- vector(
    "list",
    n_iter
  )

  capture_list <- vector(
    "list",
    n_iter
  )

  lambda_list <- vector(
    "list",
    n_iter
  )

  p_detect_list <- vector(
    "list",
    n_iter
  )

  occupancy_models <- vector(
    "list",
    n_iter
  )

  capture_models <- vector(
    "list",
    n_iter
  )

  abundance_models <- vector(
    "list",
    n_iter
  )

  lv_sites_list <- vector(
    "list",
    n_iter
  )

  lv_species_list <- vector(
    "list",
    n_iter
  )

  iteration_success <- rep(
    FALSE,
    n_iter
  )

  diagnostic_AIC <- data.frame(
    iteration = integer(),
    occupancy_AIC = numeric(),
    capture_AIC = numeric(),
    abundance_AIC = numeric(),
    abundance_logLik = numeric(),
    stringsAsFactors = FALSE
  )

  # ============================================================
  # 20. EM-like stochastic iterations
  # ============================================================

  for (i in seq_len(n_iter)) {

    if (verbose) {
      message(
        "Iteration ",
        i
      )
    }

    t_iter <- Sys.time()

    # ----------------------------------------------------------
    # 20.1 Update occupancy matrix
    # ----------------------------------------------------------

    z_matrix <- z_template
    z_matrix[,] <- 0

    idx <- cbind(
      match(
        as.character(
          reduced_data[[site_col]]
        ),
        z_sites
      ),
      match(
        as.character(
          reduced_data[[otu_col]]
        ),
        z_otus
      )
    )

    if (anyNA(idx)) {

      warning(
        "Could not align site/OTU indices at iteration ",
        i,
        ". Iteration skipped."
      )

      next
    }

    z_matrix[idx] <-
      reduced_data$z_sim

    # ----------------------------------------------------------
    # 20.2 GLLVM occupancy model
    # ----------------------------------------------------------

    model_occupancy <- tryCatch(

      gllvm::gllvm(
        y = z_matrix,
        X = X_cov,
        family = "binomial",
        num.lv = num_lv_c
      ),

      error = function(e) {

        message(
          "GLLVM failed at iteration ",
          i,
          ": ",
          e$message
        )

        NULL
      }
    )

    if (is.null(model_occupancy)) {
      next
    }

    occupancy_models[[i]] <-
      model_occupancy

    # ----------------------------------------------------------
    # 20.3 GLLVM latent site variables
    # ----------------------------------------------------------

    if (
      !is.null(model_occupancy$lvs) &&
      length(model_occupancy$lvs) > 0 &&
      ncol(
        as.data.frame(
          model_occupancy$lvs
        )
      ) > 0
    ) {

      lv_sites <- as.data.frame(
        model_occupancy$lvs
      )

      colnames(lv_sites) <-
        paste0(
          "LV",
          seq_len(
            ncol(lv_sites)
          )
        )

      site_names <- rownames(
        model_occupancy$lvs
      )

      if (
        is.null(site_names) ||
        length(site_names) != nrow(lv_sites)
      ) {

        site_names <- z_sites
      }

      lv_sites[[site_col]] <-
        as.character(site_names)

      lv_sites$Iteration <- i

    } else {

      lv_sites <- data.frame()
    }

    # ----------------------------------------------------------
    # 20.4 Taxon latent loadings
    # ----------------------------------------------------------

    theta_obj <- tryCatch(
      model_occupancy$params$theta,
      error = function(e) NULL
    )

    if (
      !is.null(theta_obj) &&
      length(theta_obj) > 0
    ) {

      theta_mat <- as.matrix(
        theta_obj
      )

      if (
        nrow(theta_mat) > 0 &&
        ncol(theta_mat) > 0
      ) {

        lv_species <- as.data.frame(
          theta_mat
        )

        colnames(lv_species) <-
          paste0(
            "LV",
            seq_len(
              ncol(lv_species)
            )
          )

        otu_names <- rownames(
          theta_mat
        )

        if (
          is.null(otu_names) ||
          length(otu_names) != nrow(lv_species)
        ) {

          otu_names <- z_otus[
            seq_len(
              min(
                length(z_otus),
                nrow(lv_species)
              )
            )
          ]
        }

        if (
          length(otu_names) ==
          nrow(lv_species)
        ) {

          lv_species[[otu_col]] <-
            as.character(otu_names)

        } else {

          lv_species[[otu_col]] <-
            as.character(
              seq_len(
                nrow(lv_species)
              )
            )
        }

        lv_species$Iteration <- i

      } else {

        lv_species <- data.frame()
      }

    } else {

      lv_species <- data.frame()
    }

    lv_sites_list[[i]] <-
      lv_sites

    lv_species_list[[i]] <-
      lv_species

    # ----------------------------------------------------------
    # 20.5 Occupancy probability
    # ----------------------------------------------------------

    psi_prob <- tryCatch(

      stats::predict(
        model_occupancy,
        type = "response"
      ),

      error = function(e) {

        message(
          "GLLVM prediction failed at iteration ",
          i,
          ": ",
          e$message
        )

        NULL
      }
    )

    if (is.null(psi_prob)) {
      next
    }

    psi_prob <- as.matrix(
      psi_prob
    )

    if (!identical(
      dim(psi_prob),
      dim(z_matrix)
    )) {

      warning(
        "GLLVM prediction dimensions do not match occupancy matrix ",
        "at iteration ",
        i,
        ". Iteration skipped."
      )

      next
    }

    rownames(psi_prob) <-
      rownames(z_matrix)

    colnames(psi_prob) <-
      colnames(z_matrix)

    psi_long <- reshape2::melt(
      psi_prob,
      varnames = c(
        site_col,
        otu_col
      ),
      value.name = "psi_prob"
    )

    psi_long[[site_col]] <-
      as.character(
        psi_long[[site_col]]
      )

    psi_long[[otu_col]] <-
      as.character(
        psi_long[[otu_col]]
      )

    psi_long$psi_prob <-
      pmin(
        pmax(
          as.numeric(
            psi_long$psi_prob
          ),
          1e-6
        ),
        1 - 1e-6
      )

    psi_long$eta <-
      stats::qlogis(
        psi_long$psi_prob
      )

    psi_list[[i]] <-
      psi_long[
        ,
        c(
          site_col,
          otu_col,
          "eta"
        ),
        drop = FALSE
      ]

    # ----------------------------------------------------------
    # 20.6 Add current Z state to sample data
    # ----------------------------------------------------------

    sample_data <- sample_data |>
      dplyr::select(
        -dplyr::any_of(
          c(
            "z_sim",
            "p0_sample",
            "capture_prob"
          )
        )
      ) |>
      dplyr::left_join(
        reduced_data |>
          dplyr::select(
            dplyr::all_of(
              c(
                site_keys,
                "z_sim"
              )
            )
          ),
        by = site_keys
      )

    capture_fit_data <- sample_data |>
      dplyr::filter(
        .data$z_sim == 1
      )

    if (nrow(capture_fit_data) == 0) {

      message(
        "No rows available for capture model at iteration ",
        i,
        ". Iteration skipped."
      )

      next
    }

    # ----------------------------------------------------------
    # 20.7 Capture GLMM
    # ----------------------------------------------------------

    cap_fit <- tryCatch(

      glmmTMB::glmmTMB(
        formula = capture_formula,
        data = capture_fit_data,
        family = stats::binomial()
      ),

      error = function(e) {

        message(
          "Capture model failed at iteration ",
          i,
          ": ",
          e$message
        )

        NULL
      }
    )

    if (is.null(cap_fit)) {
      next
    }

    capture_models[[i]] <-
      cap_fit

    eta_capture <- tryCatch(

      as.numeric(
        stats::predict(
          cap_fit,
          type = "link",
          newdata = sample_data,
          allow.new.levels = TRUE
        )
      ),

      error = function(e) {

        message(
          "Capture prediction failed at iteration ",
          i,
          ": ",
          e$message
        )

        NULL
      }
    )

    if (
      is.null(eta_capture) ||
      length(eta_capture) != nrow(sample_data)
    ) {

      next
    }

    capture_prob <-
      stats::plogis(
        eta_capture
      )

    capture_list[[i]] <-
      data.frame(
        sample_data[
          sample_keys
        ],
        eta = eta_capture,
        check.names = FALSE
      )

    # ----------------------------------------------------------
    # 20.8 Abundance data conditional on A = 1
    # ----------------------------------------------------------

    abundance_data <- long_df |>
      dplyr::left_join(
        sample_data |>
          dplyr::select(
            dplyr::all_of(
              c(
                sample_keys,
                "a_sim"
              )
            )
          ),
        by = sample_keys
      ) |>
      dplyr::filter(
        .data$a_sim == 1
      )

    if (nrow(abundance_data) == 0) {

      message(
        "No rows available for abundance model at iteration ",
        i,
        ". Iteration skipped."
      )

      next
    }

    abundance_data <- to_factor_cols(
      abundance_data,
      intersect(
        id_factor_cols,
        names(abundance_data)
      )
    )

    abundance_data <- character_to_factor(
      abundance_data,
      formula_predictors
    )

    # ----------------------------------------------------------
    # 20.9 Abundance GLMM
    # ----------------------------------------------------------

    model_abundance <- tryCatch(

      glmmTMB::glmmTMB(
        formula = abundance_formula,
        data = abundance_data,
        family = abundance_glmm_family,
        ziformula = zi_formula
      ),

      error = function(e) {

        message(
          "Abundance model failed at iteration ",
          i,
          ": ",
          e$message
        )

        NULL
      }
    )

    if (is.null(model_abundance)) {
      next
    }

    abundance_models[[i]] <-
      model_abundance

    # ----------------------------------------------------------
    # 20.10 Abundance predictions
    # ----------------------------------------------------------

    lambda_eta <- tryCatch(

      as.numeric(
        stats::predict(
          model_abundance,
          type = "link",
          newdata = long_df,
          allow.new.levels = TRUE
        )
      ),

      error = function(e) {

        message(
          "Abundance prediction failed at iteration ",
          i,
          ": ",
          e$message
        )

        NULL
      }
    )

    if (
      is.null(lambda_eta) ||
      length(lambda_eta) != nrow(long_df)
    ) {

      next
    }

    lambda_list[[i]] <-
      data.frame(
        long_df[
          site_keys
        ],
        eta = lambda_eta,
        check.names = FALSE
      ) |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            site_keys
          )
        )
      ) |>
      dplyr::summarise(
        eta = mean(
          .data$eta,
          na.rm = TRUE
        ),
        .groups = "drop"
      )

    # ----------------------------------------------------------
    # 20.11 PCR/read-level zero probabilities
    # ----------------------------------------------------------

    pred_data_all <- long_df |>
      dplyr::left_join(
        sample_data |>
          dplyr::select(
            dplyr::all_of(
              c(
                sample_keys,
                "a_sim"
              )
            )
          ),
        by = sample_keys
      )

    pred_data_all <- to_factor_cols(
      pred_data_all,
      intersect(
        id_factor_cols,
        names(pred_data_all)
      )
    )

    pred_data_all <- character_to_factor(
      pred_data_all,
      formula_predictors
    )

    p0_pcr <- tryCatch(

      get_p0_abundance(
        fit = model_abundance,
        newdata = pred_data_all,
        family_name = abundance_family
      ),

      error = function(e) {

        message(
          "Zero-probability calculation failed at iteration ",
          i,
          ": ",
          e$message
        )

        NULL
      }
    )

    if (
      is.null(p0_pcr) ||
      length(p0_pcr) != nrow(pred_data_all)
    ) {

      next
    }

    # ----------------------------------------------------------
    # 20.12 Collapse PCR zeros to sample level
    # ----------------------------------------------------------

    sample_p0 <- pred_data_all |>
      dplyr::mutate(
        p0_pcr = p0_pcr
      ) |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            sample_keys
          )
        )
      ) |>
      dplyr::summarise(
        p0_sample = prod(
          .data$p0_pcr,
          na.rm = TRUE
        ),
        .groups = "drop"
      )

    sample_data <- sample_data |>
      dplyr::select(
        -dplyr::any_of(
          c(
            "p0_sample",
            "capture_prob"
          )
        )
      ) |>
      dplyr::mutate(
        capture_prob =
          capture_prob
      ) |>
      dplyr::left_join(
        sample_p0,
        by = sample_keys
      )

    sample_data$p0_sample[
      is.na(sample_data$p0_sample)
    ] <- 1

    # ----------------------------------------------------------
    # 20.13 Update latent A
    # ----------------------------------------------------------

    zero_sample <- which(
      sample_data$a_obs == 0 &
      sample_data$z_sim == 1
    )

    if (length(zero_sample) > 0) {

      numerator_a <-
        sample_data$capture_prob[
          zero_sample
        ] *
        sample_data$p0_sample[
          zero_sample
        ]

      denominator_a <-
        (
          1 -
          sample_data$capture_prob[
            zero_sample
          ]
        ) +
        numerator_a

      posterior_a <-
        numerator_a /
        pmax(
          denominator_a,
          1e-12
        )

      posterior_a <-
        pmin(
          pmax(
            posterior_a,
            1e-4
          ),
          1 - 1e-4
        )

      sample_data$a_sim[
        zero_sample
      ] <- stats::rbinom(
        n = length(
          zero_sample
        ),
        size = 1,
        prob = posterior_a
      )
    }

    # Positive observations remain positive
    sample_data$a_sim[
      sample_data$a_obs == 1
    ] <- 1

    # Cannot capture if site absent
    sample_data$a_sim[
      sample_data$z_sim == 0
    ] <- 0

    # ----------------------------------------------------------
    # 20.14 Collapse sample probabilities to site level
    # ----------------------------------------------------------

    site_p0 <- sample_data |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            site_keys
          )
        )
      ) |>
      dplyr::summarise(
        p0_site = prod(
          (
            1 -
            .data$capture_prob
          ) +
          .data$capture_prob *
          .data$p0_sample,
          na.rm = TRUE
        ),
        .groups = "drop"
      )

    site_p0$p0_site <-
      pmin(
        pmax(
          site_p0$p0_site,
          0
        ),
        1
      )

    # ----------------------------------------------------------
    # 20.15 Detection probability on cloglog scale
    # ----------------------------------------------------------

    p_detect_list[[i]] <-
      site_p0 |>
      dplyr::mutate(
        .p0_safe = pmin(
          pmax(
            .data$p0_site,
            1e-12
          ),
          1 - 1e-12
        ),
        eta = log(
          -log(
            .data$.p0_safe
          )
        )
      ) |>
      dplyr::select(
        dplyr::all_of(
          site_keys
        ),
        .data$eta
      )

    # ----------------------------------------------------------
    # 20.16 Update site-level p0
    # ----------------------------------------------------------

    reduced_data <- reduced_data |>
      dplyr::select(
        -dplyr::any_of(
          "p0_site"
        )
      ) |>
      dplyr::left_join(
        site_p0,
        by = site_keys
      ) |>
      dplyr::mutate(
        p0_site =
          tidyr::replace_na(
            .data$p0_site,
            1
          )
      )

    # ----------------------------------------------------------
    # 20.17 Merge GLLVM psi
    # ----------------------------------------------------------

    z_merge <- reduced_data |>
      dplyr::left_join(
        psi_long[
          ,
          c(
            site_col,
            otu_col,
            "psi_prob"
          ),
          drop = FALSE
        ],
        by = site_keys
      )

    if (anyNA(z_merge$psi_prob)) {

      message(
        "Missing GLLVM occupancy probabilities at iteration ",
        i,
        ". Iteration skipped."
      )

      next
    }

    # ----------------------------------------------------------
    # 20.18 Update latent Z
    # ----------------------------------------------------------

    zero_indices <- which(
      z_merge$z_obs == 0
    )

    if (length(zero_indices) > 0) {

      numerator_z <-
        z_merge$psi_prob[
          zero_indices
        ] *
        z_merge$p0_site[
          zero_indices
        ]

      denominator_z <-
        (
          1 -
          z_merge$psi_prob[
            zero_indices
          ]
        ) +
        numerator_z

      posterior_z <-
        numerator_z /
        pmax(
          denominator_z,
          1e-12
        )

      posterior_z <-
        pmin(
          pmax(
            posterior_z,
            1e-4
          ),
          1 - 1e-4
        )

      reduced_data$z_sim[
        zero_indices
      ] <- stats::rbinom(
        n = length(
          zero_indices
        ),
        size = 1,
        prob = posterior_z
      )
    }

    # Positive site detections remain occupied
    reduced_data$z_sim[
      reduced_data$z_obs == 1
    ] <- 1

    # ----------------------------------------------------------
    # 20.19 Diagnostics
    # ----------------------------------------------------------

    diagnostic_AIC <- rbind(
      diagnostic_AIC,
      data.frame(
        iteration = i,

        occupancy_AIC =
          tryCatch(
            stats::AIC(
              model_occupancy
            ),
            error = function(e) {
              NA_real_
            }
          ),

        capture_AIC =
          tryCatch(
            stats::AIC(
              cap_fit
            ),
            error = function(e) {
              NA_real_
            }
          ),

        abundance_AIC =
          tryCatch(
            stats::AIC(
              model_abundance
            ),
            error = function(e) {
              NA_real_
            }
          ),

        abundance_logLik =
          tryCatch(
            as.numeric(
              stats::logLik(
                model_abundance
              )
            ),
            error = function(e) {
              NA_real_
            }
          )
      )
    )

    iteration_success[i] <- TRUE

    if (verbose) {

      message(
        "Iteration ",
        i,
        " finished in ",
        round(
          difftime(
            Sys.time(),
            t_iter,
            units = "secs"
          ),
          2
        ),
        " seconds"
      )
    }
  }

  # ============================================================
  # 21. Retained successful iterations
  # ============================================================

  keep_requested <- seq.int(
    burn_in + 1,
    n_iter
  )

  keep <- keep_requested[
    iteration_success[
      keep_requested
    ]
  ]

  if (length(keep) == 0) {

    stop(
      "No successful iterations remained after burn-in."
    )
  }

  if (
    verbose &&
    length(keep) <
    length(keep_requested)
  ) {

    message(
      length(keep),
      " of ",
      length(keep_requested),
      " post-burn-in iterations were fully successful."
    )
  }

  # ============================================================
  # 22. Generic link-scale summarizer
  # ============================================================

  summarise_link <- function(
    lst,
    link_type = c(
      "logit",
      "log",
      "cloglog"
    )
  ) {

    link_type <- match.arg(
      link_type
    )

    df <- dplyr::bind_rows(
      lst
    )

    if (nrow(df) == 0) {
      return(
        data.frame()
      )
    }

    inv_link <- switch(
      link_type,

      logit =
        stats::plogis,

      log =
        exp,

      cloglog =
        function(x) {
          1 - exp(-exp(x))
        }
    )

    df |>
      dplyr::mutate(
        value =
          inv_link(
            .data$eta
          )
      ) |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            site_keys
          )
        )
      ) |>
      dplyr::summarise(
        mean = mean(
          .data$value,
          na.rm = TRUE
        ),
        median = stats::median(
          .data$value,
          na.rm = TRUE
        ),
        lwr = stats::quantile(
          .data$value,
          0.025,
          na.rm = TRUE,
          names = FALSE
        ),
        upr = stats::quantile(
          .data$value,
          0.975,
          na.rm = TRUE,
          names = FALSE
        ),
        .groups = "drop"
      )
  }

  # ============================================================
  # 23. Capture summaries
  # ============================================================

  summarise_capture <- function(
    lst
  ) {

    df <- dplyr::bind_rows(
      lst
    )

    if (nrow(df) == 0) {
      return(
        data.frame()
      )
    }

    df |>
      dplyr::mutate(
        capture =
          stats::plogis(
            .data$eta
          )
      ) |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            sample_keys
          )
        )
      ) |>
      dplyr::summarise(
        capture_mean =
          mean(
            .data$capture,
            na.rm = TRUE
          ),

        capture_median =
          stats::median(
            .data$capture,
            na.rm = TRUE
          ),

        capture_lwr =
          stats::quantile(
            .data$capture,
            0.025,
            na.rm = TRUE,
            names = FALSE
          ),

        capture_upr =
          stats::quantile(
            .data$capture,
            0.975,
            na.rm = TRUE,
            names = FALSE
          ),

        .groups = "drop"
      )
  }

  # ============================================================
  # 24. Site-level capture summaries
  # ============================================================

  summarise_capture_site <- function(
    capture_summary
  ) {

    if (nrow(capture_summary) == 0) {
      return(
        data.frame()
      )
    }

    capture_summary |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            site_keys
          )
        )
      ) |>
      dplyr::summarise(
        capture_mean =
          mean(
            .data$capture_mean,
            na.rm = TRUE
          ),

        capture_median =
          stats::median(
            .data$capture_median,
            na.rm = TRUE
          ),

        capture_lwr =
          stats::quantile(
            .data$capture_lwr,
            0.025,
            na.rm = TRUE,
            names = FALSE
          ),

        capture_upr =
          stats::quantile(
            .data$capture_upr,
            0.975,
            na.rm = TRUE,
            names = FALSE
          ),

        .groups = "drop"
      )
  }

  # ============================================================
  # 25. Occupancy summaries
  # ============================================================

  psi_summary <-
    summarise_link(
      psi_list[keep],
      "logit"
    ) |>
    dplyr::rename_with(
      ~ paste0(
        "psi_",
        .x
      ),
      -dplyr::all_of(
        site_keys
      )
    )

  # ============================================================
  # 26. Abundance summaries
  # ============================================================

  lambda_summary <-
    summarise_link(
      lambda_list[keep],
      "log"
    ) |>
    dplyr::rename_with(
      ~ paste0(
        "lambda_",
        .x
      ),
      -dplyr::all_of(
        site_keys
      )
    )

  # ============================================================
  # 27. Detection summaries
  # ============================================================

  p_detect_summary <-
    summarise_link(
      p_detect_list[keep],
      "cloglog"
    ) |>
    dplyr::rename_with(
      ~ paste0(
        "p_detect_",
        .x
      ),
      -dplyr::all_of(
        site_keys
      )
    )

  # ============================================================
  # 28. Capture summaries
  # ============================================================

  capture_summary <-
    summarise_capture(
      capture_list[keep]
    )

  capture_site_summary <-
    summarise_capture_site(
      capture_summary
    )

  # ============================================================
  # 29. Final site x OTU summary
  # ============================================================

  final_summary <- psi_summary

  if (nrow(capture_site_summary) > 0) {

    final_summary <-
      dplyr::left_join(
        final_summary,
        capture_site_summary,
        by = site_keys
      )
  }

  if (nrow(lambda_summary) > 0) {

    final_summary <-
      dplyr::left_join(
        final_summary,
        lambda_summary,
        by = site_keys
      )
  }

  if (nrow(p_detect_summary) > 0) {

    final_summary <-
      dplyr::left_join(
        final_summary,
        p_detect_summary,
        by = site_keys
      )
  }

  # ============================================================
  # 30. Combine GLLVM latent variables
  # ============================================================

  lv_sites_combined <-
    dplyr::bind_rows(
      lv_sites_list[keep]
    )

  lv_species_combined <-
    dplyr::bind_rows(
      lv_species_list[keep]
    )

  # ============================================================
  # 31. Mean site latent scores
  # ============================================================

  mean_lv_sites <- if (
    nrow(lv_sites_combined) > 0 &&
    site_col %in% names(lv_sites_combined)
  ) {

    lv_sites_combined |>
      dplyr::group_by(
        .data[[site_col]]
      ) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with(
            "LV"
          ),
          ~ mean(
            .x,
            na.rm = TRUE
          )
        ),
        .groups = "drop"
      )

  } else {

    data.frame()
  }

  # ============================================================
  # 32. Mean OTU latent loadings
  # ============================================================

  mean_lv_species <- if (
    nrow(lv_species_combined) > 0 &&
    otu_col %in% names(lv_species_combined)
  ) {

    lv_species_combined |>
      dplyr::group_by(
        .data[[otu_col]]
      ) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with(
            "LV"
          ),
          ~ mean(
            .x,
            na.rm = TRUE
          )
        ),
        .groups = "drop"
      )

  } else {

    data.frame()
  }

  # ============================================================
  # 33. Return results
  # ============================================================

  list(

    summary = final_summary,

    capture = capture_summary,

    capture_site =
      capture_site_summary,

    psi_list =
      psi_list[keep],

    capture_list =
      capture_list[keep],

    lambda_list =
      lambda_list[keep],

    p_detect_list =
      p_detect_list[keep],

    occupancy_models =
      occupancy_models[keep],

    capture_models =
      capture_models[keep],

    abundance_models =
      abundance_models[keep],

    reduced_data =
      reduced_data,

    sample_data =
      sample_data,

    long_df =
      long_df,

    lv_sites =
      lv_sites_combined,

    lv_species =
      lv_species_combined,

    mean_lv_sites =
      mean_lv_sites,

    mean_lv_species =
      mean_lv_species,

    filter_summary = list(
      otu_stats =
        otu_stats,

      kept_otus =
        keep_otus,

      min_species_sum =
        min_species_sum,

      min_detection_replicates =
        min_detection_replicates,

      abundance_threshold =
        abundance_threshold
    ),

    diagnostic_AIC =
      diagnostic_AIC,

    iteration_success =
      iteration_success,

    retained_iterations =
      keep,

    abundance_formula =
      abundance_formula,

    capture_formula =
      capture_formula,

    note = paste(
      "Three-level GLLVM-GLMM model.",
      "Hierarchy: site occupancy Z estimated jointly across taxa using gllvm,",
      "biological-sample capture A estimated using a binomial GLMM,",
      "and replicate-level abundance Y estimated using a count GLMM."
    )
  )
}
