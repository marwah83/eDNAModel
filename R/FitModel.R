#' Fit Hierarchical Occupancy–Capture–Abundance Models for eDNA Data
#'
#' Fits a hierarchical model for environmental DNA (eDNA) and metabarcoding
#' count data using an iterative stochastic EM-like procedure. The model
#' represents three linked processes: species occupancy at the site level
#' (\eqn{\psi}), capture of DNA in biological samples (\eqn{p}), and sequencing
#' abundance in PCR replicates (\eqn{\lambda}). Each model component is fitted
#' using generalized linear mixed models implemented in \code{glmmTMB}, allowing
#' separate fixed effects, random effects, offsets, and distributional
#' assumptions for the occupancy, capture, and abundance processes. Latent
#' occupancy and capture states are updated iteratively using probabilities
#' derived from the fitted component models. Fixed-effect uncertainty is
#' combined across retained iterations using Rubin-style pooling, incorporating
#' both the average within-fit covariance and the between-iteration covariance.
#' Prediction uncertainty for occupancy, capture, and abundance is propagated
#' from the fitted linear predictors and their standard errors and transformed
#' to the natural scale. Supported abundance distributions include Poisson,
#' Negative Binomial, Zero-Inflated Poisson (ZIP), and Zero-Inflated Negative
#' Binomial (ZINB). Random effects and sequencing-depth offsets can be specified
#' through the underlying \code{glmmTMB} model formulas.
#'
#' @param phyloseq A \code{phyloseq} object containing OTU counts and sample
#'   metadata.
#' @param site_col Character. Column identifying the site-level occupancy unit.
#' @param sample_col Character. Column identifying biological samples.
#' @param replicate_col Character or \code{NULL}. Column identifying PCR
#'   replicates.
#' @param occupancy_formula Formula describing occupancy probability
#'   (\eqn{\psi}); the response must be \code{z_sim}.
#' @param capture_formula Formula describing capture probability
#'   (\eqn{p}); the response must be \code{a_sim}.
#' @param abundance_formula Formula describing sequencing abundance
#'   (\eqn{\lambda}); the response must match \code{count_col}.
#' @param otu_col Character. OTU identifier column.
#' @param count_col Character. Count column containing sequencing read counts.
#' @param min_species_sum Minimum total count required to retain an OTU.
#' @param min_detection_replicates Minimum number of positive observations
#'   required to retain an OTU.
#' @param abundance_threshold Threshold above which a count is treated as a
#'   detection.
#' @param n_iter Number of stochastic fitting iterations.
#' @param burn_in Number of initial iterations discarded before pooling.
#' @param abundance_family One of \code{"poisson"}, \code{"nbinom"},
#'   \code{"zip"}, or \code{"zinb"}.
#' @param verbose Logical. If \code{TRUE}, print fitting progress.
#' @param n_sim_ci Number of simulations used for prediction-scale uncertainty
#'   propagation.
#'
#' @return A list containing
#' \describe{
#'   \item{psi}{Summaries of estimated occupancy probabilities.}
#'   \item{capture}{Summaries of estimated capture probabilities.}
#'   \item{lambda}{Summaries of estimated sequencing abundance.}
#'   \item{p_detect}{Detection probabilities implied by the abundance model.}
#'   \item{beta}{Pooled fixed-effect estimates, standard errors, confidence
#'     intervals, and variance components on the link scale.}
#'   \item{beta_natural}{Intercept estimates transformed to the natural scale
#'     where appropriate.}
#'   \item{beta_covariance}{Within-iteration, between-iteration, and total
#'     pooled fixed-effect covariance matrices, when returned by the fitted
#'     implementation.}
#'   \item{beta_matrices}{Fixed-effect estimates from retained stochastic
#'     iterations, when returned by the fitted implementation.}
#'   \item{occ_fit, cap_fit, abund_fit}{Final fitted component GLMMs.}
#'   \item{occupancy_models, capture_models, abundance_models}{Fitted component
#'     models from retained iterations.}
#'   \item{site_data, sample_data, long_df}{Processed datasets used internally
#'     for model fitting.}
#'   \item{filter_summary}{Information on OTU filtering.}
#'   \item{diagnostic_AIC}{Component-wise AIC values across iterations.}
#' }
#'
#' @details
#' The hierarchical model is
#'
#' \deqn{
#' Z_i \sim \mathrm{Bernoulli}(\psi_i)
#' }
#'
#' \deqn{
#' A_{ij} \mid Z_i \sim \mathrm{Bernoulli}(p_{ij} Z_i)
#' }
#'
#' \deqn{
#' Y_{ijk} \mid A_{ij} \sim
#' \mathrm{Count}(\lambda_{ijk} A_{ij})
#' }
#'
#' where \eqn{i} indexes sites, \eqn{j} indexes biological samples, and
#' \eqn{k} indexes PCR replicates. The abundance distribution may be Poisson,
#' Negative Binomial, ZIP, or ZINB.
#'
#' The abundance component also determines the probability of obtaining a
#' positive sequencing observation conditional on capture:
#'
#' \deqn{
#' p_{\mathrm{detect}}
#' =
#' 1 - P(Y=0 \mid A=1).
#' }
#'
#' For the Poisson model,
#'
#' \deqn{
#' p_{\mathrm{detect}}
#' =
#' 1-\exp(-\lambda).
#' }
#'
#' For NB, ZIP, and ZINB models, the corresponding model-specific zero
#' probability is used. Thus, \code{p_detect} is a derived quantity rather
#' than an independently estimated model parameter.
#'
#' @section Stochastic EM-like algorithm:
#' Each iteration:
#' \enumerate{
#'   \item Fits the occupancy model using the current latent occupancy states
#'     \eqn{Z}.
#'   \item Fits the capture model conditional on \eqn{Z=1}.
#'   \item Fits the abundance model conditional on \eqn{A=1}.
#'   \item Computes zero and detection probabilities from the abundance model.
#'   \item Updates latent capture states \eqn{A}.
#'   \item Updates latent occupancy states \eqn{Z}.
#' }
#'
#' @section Fixed-effect uncertainty:
#' For each retained iteration \eqn{m}, let
#' \eqn{\hat{\boldsymbol{\beta}}_m} denote the estimated fixed-effect vector and
#' \eqn{\mathbf{U}_m} its estimated covariance matrix. The pooled fixed-effect
#' estimate is
#'
#' \deqn{
#' \bar{\boldsymbol{\beta}}
#' =
#' \frac{1}{M}
#' \sum_{m=1}^{M}
#' \hat{\boldsymbol{\beta}}_m.
#' }
#'
#' The average within-iteration covariance is
#'
#' \deqn{
#' \bar{\mathbf{U}}
#' =
#' \frac{1}{M}
#' \sum_{m=1}^{M}
#' \mathbf{U}_m,
#' }
#'
#' and the between-iteration covariance is
#'
#' \deqn{
#' \mathbf{B}
#' =
#' \frac{1}{M-1}
#' \sum_{m=1}^{M}
#' (\hat{\boldsymbol{\beta}}_m-\bar{\boldsymbol{\beta}})
#' (\hat{\boldsymbol{\beta}}_m-\bar{\boldsymbol{\beta}})^\top.
#' }
#'
#' The Rubin-style total covariance is
#'
#' \deqn{
#' \mathbf{T}
#' =
#' \bar{\mathbf{U}}
#' +
#' \left(1+\frac{1}{M}\right)\mathbf{B}.
#' }
#'
#' This incorporates uncertainty within the individual GLMM fits together
#' with variability caused by the stochastic latent-state updates.
#'
#' @section Prediction uncertainty:
#' Prediction uncertainty for occupancy, capture, and abundance is propagated
#' on the linear-predictor scale using an approximate normal distribution,
#'
#' \deqn{
#' \eta^{*}
#' \sim
#' \mathcal{N}(\hat{\eta}, \widehat{\mathrm{SE}}_{\eta}^{\,2}),
#' }
#'
#' followed by transformation through the appropriate inverse-link function.
#'
#' @section Model features:
#' \itemize{
#'   \item Explicit three-level eDNA hierarchy.
#'   \item Separate occupancy, capture, and abundance regression models.
#'   \item Fixed and random effects through \code{glmmTMB}.
#'   \item Poisson, NB, ZIP, and ZINB abundance distributions.
#'   \item Sequencing-depth and other offsets.
#'   \item Pooled fixed-effect covariance across stochastic iterations.
#' }
#'
#' @section Caveats:
#' \itemize{
#'   \item The method is an approximate stochastic latent-state procedure and
#'     does not maximize a single observed-data joint likelihood.
#'   \item Rubin-style pooling is used as an approximation for uncertainty
#'     propagation across stochastic iterations.
#'   \item Component-wise AIC values are diagnostic and are not equivalent to
#'     AIC from a single joint likelihood.
#'   \item Detection probabilities depend on the assumed abundance
#'     distribution.
#' }
#'
#' @examples
#' \dontrun{
#' fit <- FitModel(
#'   phyloseq = ps,
#'   site_col = "Site",
#'   sample_col = "Sample",
#'   replicate_col = "Replicate",
#'   otu_col = "OTU",
#'   count_col = "y",
#'   occupancy_formula = z_sim ~ 1 + (1 | OTU),
#'   capture_formula = a_sim ~ 1 + (1 | OTU),
#'   abundance_formula = y ~ offset(log(total_reads)) + (1 | OTU),
#'   abundance_family = "nbinom",
#'   n_iter = 20,
#'   burn_in = 5
#' )
#'
#' fit$beta
#' head(fit$psi)
#' head(fit$capture)
#' head(fit$lambda)
#' }
#'
#' @importFrom glmmTMB glmmTMB nbinom2
#' @importFrom stats plogis predict rbinom quantile median poisson
#' @import dplyr
#' @import rlang
#'
#' @export
FitModel <- function(
    phyloseq,
    site_col,
    sample_col,
    replicate_col = NULL,
    occupancy_formula,
    capture_formula,
    abundance_formula,
    otu_col = NULL,
    count_col = NULL,
    min_species_sum = 10,
    min_detection_replicates = 1,
    abundance_threshold = 0,
    n_iter = 20,
    burn_in = 5,
    abundance_family = c(
      "poisson",
      "nbinom",
      "zip",
      "zinb"
    ),
    verbose = TRUE,
    n_sim_ci = 100
) {

  # ============================================================
  # Basic checks
  # ============================================================

  abundance_family <- match.arg(
    abundance_family
  )

  if (is.null(otu_col)) {
    stop("Please specify otu_col.")
  }

  if (is.null(count_col)) {
    stop("Please specify count_col.")
  }

  if (burn_in >= n_iter) {
    stop("burn_in must be < n_iter.")
  }


  # ============================================================
  # Helper: formula variables
  # ============================================================

  get_formula_vars <- function(
      formula,
      response
  ) {

    setdiff(
      all.vars(formula),
      response
    )
  }


  # ============================================================
  # NEW HELPER:
  # Extract fixed-effect estimates and model-based variances
  # ============================================================

  extract_beta_glmmTMB <- function(
      fit,
      component = c(
        "cond",
        "zi"
      )
  ) {

    component <- match.arg(
      component
    )

    if (component == "cond") {

      sm <- tryCatch(
        summary(fit)$coefficients$cond,
        error = function(e) NULL
      )

    } else {

      sm <- tryCatch(
        summary(fit)$coefficients$zi,
        error = function(e) NULL
      )
    }

    if (
      is.null(sm) ||
      length(sm) == 0 ||
      nrow(sm) == 0
    ) {
      return(
        data.frame()
      )
    }

    sm <- as.data.frame(
      sm
    )

    if (
      !all(
        c(
          "Estimate",
          "Std. Error"
        ) %in% names(sm)
      )
    ) {
      return(
        data.frame()
      )
    }

    out <- data.frame(
      term = rownames(sm),

      estimate =
        as.numeric(
          sm$Estimate
        ),

      se =
        as.numeric(
          sm$`Std. Error`
        ),

      variance =
        as.numeric(
          sm$`Std. Error`
        )^2,

      stringsAsFactors = FALSE
    )

    rownames(out) <- NULL

    out
  }


  # ============================================================
  # NEW HELPER:
  # Rubin-style pooling of beta estimates
  #
  # beta_bar = mean(beta_r)
  #
  # U_bar = mean(within-iteration variance)
  #
  # B = between-iteration variance
  #
  # T = U_bar + (1 + 1/M) B
  #
  # SE = sqrt(T)
  # ============================================================

  pool_beta_MI <- function(
      beta_df
  ) {

    required_names <- c(
      "component",
      "term",
      "estimate",
      "variance"
    )

    if (
      is.null(beta_df) ||
      nrow(beta_df) == 0 ||
      !all(
        required_names %in%
        names(beta_df)
      )
    ) {

      return(
        data.frame()
      )
    }

    beta_df |>
      dplyr::filter(
        is.finite(.data$estimate),
        is.finite(.data$variance),
        .data$variance >= 0
      ) |>
      dplyr::group_by(
        .data$component,
        .data$term
      ) |>
      dplyr::group_modify(
        function(.x, .y) {

          q <- .x$estimate
          u <- .x$variance

          M <- length(q)

          if (M == 0) {
            return(
              data.frame()
            )
          }

          # ----------------------------------------------------
          # Pooled point estimate
          # ----------------------------------------------------

          beta_bar <- mean(
            q,
            na.rm = TRUE
          )

          # ----------------------------------------------------
          # Average within-imputation variance
          # ----------------------------------------------------

          U_bar <- mean(
            u,
            na.rm = TRUE
          )

          # ----------------------------------------------------
          # Between-imputation variance
          # ----------------------------------------------------

          B <- if (M > 1) {

            stats::var(
              q,
              na.rm = TRUE
            )

          } else {

            0
          }

          if (!is.finite(B)) {
            B <- 0
          }

          # ----------------------------------------------------
          # Total pooled variance
          # ----------------------------------------------------

          total_variance <-
            U_bar +
            (1 + 1 / M) * B

          pooled_se <- sqrt(
            total_variance
          )


          # ----------------------------------------------------
          # Rubin-style relative increase in variance
          # ----------------------------------------------------

          if (
            is.finite(U_bar) &&
            U_bar > 0
          ) {

            r <- (
              (1 + 1 / M) * B
            ) / U_bar

          } else {

            r <- Inf
          }


          # ----------------------------------------------------
          # Approximate Rubin degrees of freedom
          # ----------------------------------------------------

          if (
            M <= 1 ||
            B <= .Machine$double.eps
          ) {

            df <- Inf

          } else if (
            !is.finite(r)
          ) {

            df <- M - 1

          } else {

            df <-
              (M - 1) *
              (1 + 1 / r)^2
          }


          # ----------------------------------------------------
          # Critical value
          # ----------------------------------------------------

          crit <- if (
            is.finite(df)
          ) {

            stats::qt(
              0.975,
              df = df
            )

          } else {

            stats::qnorm(
              0.975
            )
          }


          # ----------------------------------------------------
          # 95% pooled CI
          # ----------------------------------------------------

          lower <-
            beta_bar -
            crit * pooled_se

          upper <-
            beta_bar +
            crit * pooled_se


          # ----------------------------------------------------
          # Fraction of total variance due to latent-state
          # / between-imputation variation
          # ----------------------------------------------------

          missing_information <- if (
            is.finite(total_variance) &&
            total_variance > 0
          ) {

            (
              (1 + 1 / M) * B
            ) /
              total_variance

          } else {

            NA_real_
          }


          data.frame(
            M = M,

            estimate =
              beta_bar,

            within_variance =
              U_bar,

            between_variance =
              B,

            total_variance =
              total_variance,

            std_error =
              pooled_se,

            df =
              df,

            lower =
              lower,

            upper =
              upper,

            fraction_missing_information =
              missing_information
          )
        }
      ) |>
      dplyr::ungroup()
  }


  # ============================================================
  # Prepare long data
  # ============================================================

  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col,
    nested_cols = unique(
      c(
        sample_col,
        replicate_col
      )
    )
  )

  long_df <- prep$long_df


  # ============================================================
  # Validate columns
  # ============================================================

  needed_cols <- c(
    site_col,
    sample_col,
    otu_col,
    count_col
  )

  missing_cols <- setdiff(
    needed_cols,
    names(long_df)
  )

  if (
    length(missing_cols) > 0
  ) {

    stop(
      "Missing columns in long_df: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }


  # ============================================================
  # Add sample metadata
  # ============================================================

  sample_meta <- data.frame(
    phyloseq::sample_data(
      phyloseq
    ),
    check.names = FALSE
  )

  sample_meta$.__sample_id__ <-
    rownames(
      sample_meta
    )

  if (
    !(sample_col %in%
      names(sample_meta))
  ) {

    sample_meta[[sample_col]] <-
      sample_meta$.__sample_id__
  }

  sample_meta[[sample_col]] <-
    as.character(
      sample_meta[[sample_col]]
    )

  long_df[[sample_col]] <-
    as.character(
      long_df[[sample_col]]
    )

  meta_cols_to_add <- setdiff(
    names(sample_meta),
    names(long_df)
  )

  meta_cols_to_add <- setdiff(
    meta_cols_to_add,
    ".__sample_id__"
  )

  if (
    length(meta_cols_to_add) > 0
  ) {

    long_df <- dplyr::left_join(
      long_df,

      sample_meta[
        ,
        c(
          sample_col,
          meta_cols_to_add
        ),
        drop = FALSE
      ],

      by = sample_col
    )
  }


  # ============================================================
  # Type handling
  # ============================================================

  long_df[[site_col]] <-
    as.character(
      long_df[[site_col]]
    )

  long_df[[otu_col]] <-
    as.character(
      long_df[[otu_col]]
    )

  long_df[[sample_col]] <-
    as.character(
      long_df[[sample_col]]
    )

  long_df[[count_col]] <-
    as.numeric(
      long_df[[count_col]]
    )

  if (
    !is.null(replicate_col) &&
    replicate_col %in%
      names(long_df)
  ) {

    long_df[[replicate_col]] <-
      as.character(
        long_df[[replicate_col]]
      )
  }

  long_df <- long_df[
    !is.na(
      long_df[[count_col]]
    ),
    ,
    drop = FALSE
  ]


  # ============================================================
  # Compute total_reads
  # ============================================================

  total_reads_df <- long_df |>
    dplyr::group_by(
      .data[[sample_col]]
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
    by = sample_col
  )

  uses_offset <- any(
    grepl(
      "offset",
      deparse(
        abundance_formula
      )
    )
  )

  if (uses_offset) {

    if (verbose) {
      message(
        "Offset detected -> cleaning zero-read samples..."
      )
    }

    zero_samples <- unique(
      long_df[[sample_col]][
        long_df$total_reads <= 0
      ]
    )

    if (
      length(zero_samples) > 0
    ) {

      long_df <- long_df |>
        dplyr::filter(
          !(
            .data[[sample_col]] %in%
              zero_samples
          )
        )

      if (verbose) {

        message(
          length(zero_samples),
          " zero-read samples removed."
        )
      }
    }

    if (
      any(
        long_df$total_reads <= 0,
        na.rm = TRUE
      )
    ) {

      stop(
        "total_reads must be positive when using offset(log(total_reads))."
      )
    }

  } else {

    if (verbose) {
      message(
        "No offset used."
      )
    }
  }


  # ============================================================
  # Formula variables
  # ============================================================

  occ_vars <- get_formula_vars(
    occupancy_formula,
    "z_sim"
  )

  cap_vars <- get_formula_vars(
    capture_formula,
    "a_sim"
  )

  abund_vars <- get_formula_vars(
    abundance_formula,
    count_col
  )

  all_formula_vars <- unique(
    c(
      occ_vars,
      cap_vars,
      abund_vars
    )
  )

  all_formula_vars <- setdiff(
    all_formula_vars,
    c(
      count_col,
      "total_reads"
    )
  )

  missing_formula_vars <- setdiff(
    all_formula_vars,
    names(long_df)
  )

  if (
    length(
      missing_formula_vars
    ) > 0
  ) {

    stop(
      "The following formula variables are missing from long_df: ",
      paste(
        missing_formula_vars,
        collapse = ", "
      )
    )
  }

  for (
    v in all_formula_vars
  ) {

    if (
      v %in% names(long_df) &&
      is.character(
        long_df[[v]]
      )
    ) {

      long_df[[v]] <-
        factor(
          long_df[[v]]
        )
    }
  }


  # ============================================================
  # OTU filtering
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
        .data[[count_col]] >
          abundance_threshold,
        na.rm = TRUE
      ),

      .groups = "drop"
    )

  keep_otus <- otu_stats |>
    dplyr::filter(
      .data$total_count >=
        min_species_sum,

      .data$detected_replicates >=
        min_detection_replicates
    ) |>
    dplyr::pull(
      .data[[otu_col]]
    )

  long_df <- long_df |>
    dplyr::filter(
      .data[[otu_col]] %in%
        keep_otus
    )

  if (
    nrow(long_df) == 0
  ) {

    stop(
      "No OTUs remain after filtering."
    )
  }

  long_df[[site_col]] <-
    factor(
      long_df[[site_col]]
    )

  long_df[[otu_col]] <-
    factor(
      long_df[[otu_col]]
    )

  long_df[[sample_col]] <-
    factor(
      long_df[[sample_col]]
    )

  if (
    !is.null(replicate_col) &&
    replicate_col %in%
      names(long_df)
  ) {

    long_df[[replicate_col]] <-
      factor(
        long_df[[replicate_col]]
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
  # Hierarchy keys
  # ============================================================

  site_keys <- c(
    site_col,
    otu_col
  )

  sample_keys <- c(
    site_col,
    sample_col,
    otu_col
  )

  pcr_keys <- c(
    site_col,
    sample_col,
    otu_col
  )

  if (
    !is.null(replicate_col) &&
    replicate_col %in%
      names(long_df)
  ) {

    pcr_keys <- c(
      site_col,
      sample_col,
      replicate_col,
      otu_col
    )
  }


  # ============================================================
  # Build site-level Z data
  # ============================================================

  site_keep_vars <- setdiff(
    occ_vars,
    c(
      site_keys,
      "z_sim"
    )
  )

  site_keep_vars <- intersect(
    site_keep_vars,
    names(long_df)
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

  site_data <- long_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          site_keys
        )
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
        dplyr::all_of(
          site_keep_vars
        ),
        ~ dplyr::first(.x)
      ),

      .groups = "drop"
    ) |>
    dplyr::mutate(
      z_sim =
        .data$z_obs
    )


  # ============================================================
  # Build sample-level A data
  # ============================================================

  sample_data <- long_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          sample_keys
        )
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
        dplyr::all_of(
          sample_keep_vars
        ),
        ~ dplyr::first(.x)
      ),

      .groups = "drop"
    ) |>
    dplyr::mutate(
      a_sim =
        .data$a_obs
    )


  # ============================================================
  # Validate formulas
  # ============================================================

  validate_formula <- function(
      formula,
      data,
      expected_response,
      model_name
  ) {

    if (
      !inherits(
        formula,
        "formula"
      )
    ) {

      stop(
        model_name,
        " formula must be a valid formula."
      )
    }

    response <- all.vars(
      formula
    )[1]

    if (
      response !=
        expected_response
    ) {

      stop(
        model_name,
        " formula response must be '",
        expected_response,
        "', but got '",
        response,
        "'."
      )
    }

    vars <- all.vars(
      formula
    )

    missing_vars <- setdiff(
      vars,
      names(data)
    )

    if (
      length(missing_vars) > 0
    ) {

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
    occupancy_formula,
    site_data,
    "z_sim",
    "occupancy"
  )

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
  # Family setup
  # ============================================================

  fam <- switch(
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
    abundance_family %in%
      c(
        "zip",
        "zinb"
      )
  ) {

    ~ 1

  } else {

    ~ 0
  }


  # ============================================================
  # Helper: PCR-level zero probability
  # ============================================================

  get_p0_pcr <- function(
      fit,
      newdata,
      family_name
  ) {

    mu <- as.numeric(
      predict(
        fit,
        type = "response",
        newdata = newdata,
        allow.new.levels = TRUE
      )
    )

    zi <- rep(
      0,
      length(mu)
    )

    if (
      family_name %in%
        c(
          "zip",
          "zinb"
        )
    ) {

      zi <- tryCatch(
        as.numeric(
          predict(
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
          stats::sigma(
            fit
          ),
          error = function(e)
            NA_real_
        )

        if (
          !is.finite(theta) ||
          theta <= 0
        ) {

          stats::dpois(
            0,
            lambda = mu
          )

        } else {

          stats::dnbinom(
            0,
            mu = mu,
            size = theta
          )
        }
      },

      zinb = {

        theta <- tryCatch(
          stats::sigma(
            fit
          ),
          error = function(e)
            NA_real_
        )

        if (
          !is.finite(theta) ||
          theta <= 0
        ) {

          stats::dpois(
            0,
            lambda = mu
          )

        } else {

          stats::dnbinom(
            0,
            mu = mu,
            size = theta
          )
        }
      }
    )


    p0 <-
      zi +
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
  # Storage
  # ============================================================

  psi_list <- list()
  capture_list <- list()
  lambda_list <- list()
  p_detect_list <- list()

  occupancy_models <- list()
  capture_models <- list()
  abundance_models <- list()


  # ------------------------------------------------------------
  # NEW:
  # Fixed-effect coefficients from each completed-data GLMM
  # ------------------------------------------------------------

  beta_occ_list <- list()
  beta_capture_list <- list()
  beta_abundance_list <- list()
  beta_zi_list <- list()


  diagnostic_AIC <- data.frame(
    iteration = integer(),
    occupancy_AIC = numeric(),
    capture_AIC = numeric(),
    abundance_AIC = numeric()
  )


  # ============================================================
  # EM-like / stochastic-imputation iterations
  # ============================================================

  for (
    i in seq_len(
      n_iter
    )
  ) {

    if (verbose) {
      message(
        "Iteration ",
        i
      )
    }


    # ==========================================================
    # 1. Occupancy model
    # ==========================================================

    occ_fit <- glmmTMB::glmmTMB(
      formula =
        occupancy_formula,

      data =
        site_data,

      family =
        stats::binomial()
    )


    pred_psi <- predict(
      occ_fit,

      type =
        "link",

      newdata =
        site_data,

      allow.new.levels =
        TRUE,

      se.fit =
        TRUE
    )


    eta_psi <-
      as.numeric(
        pred_psi$fit
      )

    psi <-
      stats::plogis(
        eta_psi
      )


    psi_list[[i]] <-
      data.frame(
        site_data[
          site_keys
        ],

        eta =
          eta_psi,

        se =
          as.numeric(
            pred_psi$se.fit
          )
      )


    occupancy_models[[i]] <-
      occ_fit


    # ----------------------------------------------------------
    # NEW:
    # Store occupancy beta + within-fit variance
    # ----------------------------------------------------------

    beta_occ_i <-
      extract_beta_glmmTMB(
        occ_fit,
        component = "cond"
      )

    if (
      nrow(beta_occ_i) > 0
    ) {

      beta_occ_i$iteration <- i
      beta_occ_i$component <-
        "occupancy"

      beta_occ_list[[i]] <-
        beta_occ_i
    }


    # ==========================================================
    # 2. Capture model
    # ==========================================================

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
        site_data |>
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


    if (
      nrow(
        capture_fit_data
      ) == 0
    ) {

      stop(
        "No rows available for capture model at iteration ",
        i
      )
    }


    cap_fit <- glmmTMB::glmmTMB(
      formula =
        capture_formula,

      data =
        capture_fit_data,

      family =
        stats::binomial()
    )


    pred_capture <- predict(
      cap_fit,

      type =
        "link",

      newdata =
        sample_data,

      allow.new.levels =
        TRUE,

      se.fit =
        TRUE
    )


    eta_capture <-
      as.numeric(
        pred_capture$fit
      )

    capture_prob <-
      stats::plogis(
        eta_capture
      )


    capture_list[[i]] <-
      data.frame(
        sample_data[
          sample_keys
        ],

        eta =
          eta_capture,

        se =
          as.numeric(
            pred_capture$se.fit
          )
      )


    capture_models[[i]] <-
      cap_fit


    # ----------------------------------------------------------
    # NEW:
    # Store capture beta + within-fit variance
    # ----------------------------------------------------------

    beta_capture_i <-
      extract_beta_glmmTMB(
        cap_fit,
        component = "cond"
      )

    if (
      nrow(beta_capture_i) > 0
    ) {

      beta_capture_i$iteration <- i

      beta_capture_i$component <-
        "capture"

      beta_capture_list[[i]] <-
        beta_capture_i
    }


    # ==========================================================
    # 3. Abundance model
    # ==========================================================

    abund_fit_data <- long_df |>
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


    if (
      nrow(
        abund_fit_data
      ) == 0
    ) {

      stop(
        "No rows available for abundance model at iteration ",
        i
      )
    }


    abund_fit <- glmmTMB::glmmTMB(
      formula =
        abundance_formula,

      data =
        abund_fit_data,

      family =
        fam,

      ziformula =
        zi_formula
    )


    pred_lambda <- predict(
      abund_fit,

      type =
        "link",

      newdata =
        long_df,

      allow.new.levels =
        TRUE,

      se.fit =
        TRUE
    )


    eta_lambda_all <-
      as.numeric(
        pred_lambda$fit
      )

    lambda_all <-
      exp(
        eta_lambda_all
      )


    lambda_list[[i]] <-
      data.frame(
        long_df[
          pcr_keys
        ],

        eta =
          eta_lambda_all,

        se =
          as.numeric(
            pred_lambda$se.fit
          )
      )


    p0_pcr <- get_p0_pcr(
      fit =
        abund_fit,

      newdata =
        long_df,

      family_name =
        abundance_family
    )


    p_detect_eta <-
      log(
        -log(
          pmax(
            p0_pcr,
            1e-12
          )
        )
      )


    p_detect_list[[i]] <-
      data.frame(
        long_df[
          pcr_keys
        ],

        eta =
          p_detect_eta,

        se =
          NA_real_
      )


    abundance_models[[i]] <-
      abund_fit


    # ----------------------------------------------------------
    # NEW:
    # Store abundance beta + within-fit variance
    # ----------------------------------------------------------

    beta_abundance_i <-
      extract_beta_glmmTMB(
        abund_fit,
        component = "cond"
      )

    if (
      nrow(
        beta_abundance_i
      ) > 0
    ) {

      beta_abundance_i$iteration <- i

      beta_abundance_i$component <-
        "abundance"

      beta_abundance_list[[i]] <-
        beta_abundance_i
    }


    # ----------------------------------------------------------
    # NEW:
    # Store zero-inflation beta if ZIP/ZINB
    # ----------------------------------------------------------

    if (
      abundance_family %in%
        c(
          "zip",
          "zinb"
        )
    ) {

      beta_zi_i <-
        extract_beta_glmmTMB(
          abund_fit,
          component = "zi"
        )

      if (
        nrow(beta_zi_i) > 0
      ) {

        beta_zi_i$iteration <- i

        beta_zi_i$component <-
          "zero_inflation"

        beta_zi_list[[i]] <-
          beta_zi_i
      }
    }


    # ==========================================================
    # 4. Collapse PCR probabilities to sample level
    # ==========================================================

    sample_p0 <- long_df |>
      dplyr::mutate(
        p0_pcr =
          p0_pcr
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
      is.na(
        sample_data$p0_sample
      )
    ] <- 1


    # ==========================================================
    # 5. Update latent capture state A
    # ==========================================================

    zero_sample <- which(
      sample_data$a_obs == 0 &
      sample_data$z_sim == 1
    )


    if (
      length(zero_sample) > 0
    ) {

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
        sample_data$capture_prob[
          zero_sample
        ] *
        sample_data$p0_sample[
          zero_sample
        ]


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
            0.001
          ),
          0.999
        )


      sample_data$a_sim[
        zero_sample
      ] <-
        stats::rbinom(
          n =
            length(
              zero_sample
            ),

          size =
            1,

          prob =
            posterior_a
        )
    }


    sample_data$a_sim[
      sample_data$a_obs == 1
    ] <- 1


    sample_data$a_sim[
      sample_data$z_sim == 0
    ] <- 0


    # ==========================================================
    # 6. Update latent occupancy Z
    # ==========================================================

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


    site_data <- site_data |>
      dplyr::select(
        -dplyr::any_of(
          "p0_site"
        )
      ) |>
      dplyr::left_join(
        site_p0,
        by = site_keys
      )


    site_data$p0_site[
      is.na(
        site_data$p0_site
      )
    ] <- 1


    zero_site <- which(
      site_data$z_obs == 0
    )


    if (
      length(zero_site) > 0
    ) {

      numerator_z <-
        psi[
          zero_site
        ] *
        site_data$p0_site[
          zero_site
        ]


      denominator_z <-
        (
          1 -
          psi[
            zero_site
          ]
        ) +
        psi[
          zero_site
        ] *
        site_data$p0_site[
          zero_site
        ]


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
            0.001
          ),
          0.999
        )


      site_data$z_sim[
        zero_site
      ] <-
        stats::rbinom(
          n =
            length(
              zero_site
            ),

          size =
            1,

          prob =
            posterior_z
        )
    }


    site_data$z_sim[
      site_data$z_obs == 1
    ] <- 1


    # ==========================================================
    # Diagnostic AIC
    # ==========================================================

    diagnostic_AIC <- rbind(
      diagnostic_AIC,

      data.frame(
        iteration =
          i,

        occupancy_AIC =
          stats::AIC(
            occ_fit
          ),

        capture_AIC =
          stats::AIC(
            cap_fit
          ),

        abundance_AIC =
          stats::AIC(
            abund_fit
          )
      )
    )
  }


  # ============================================================
  # Existing prediction summaries
  # ============================================================

  summarise_link <- function(
      lst,
      link_name,
      prefix,
      n_sim = n_sim_ci
  ) {

    inv_link <- switch(
      link_name,

      logit =
        stats::plogis,

      log =
        exp,

      cloglog =
        function(x) {
          1 -
            exp(
              -exp(x)
            )
        },

      stop(
        "Unknown link_name: ",
        link_name
      )
    )


    df <- dplyr::bind_rows(
      lst
    )


    if (
      nrow(df) == 0
    ) {

      return(
        data.frame()
      )
    }


    if (
      !(
        "eta" %in%
          names(df)
      )
    ) {

      stop(
        "Each list element must contain an eta column."
      )
    }


    if (
      !(
        "se" %in%
          names(df)
      )
    ) {

      df$se <- NA_real_
    }


    keys <- setdiff(
      names(df),
      c(
        "eta",
        "se"
      )
    )


    out <- df |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            keys
          )
        )
      ) |>
      dplyr::group_modify(
        function(
            .x,
            .y
        ) {

          eta_hat <-
            .x$eta

          eta_se <-
            .x$se


          if (
            all(
              is.na(
                eta_se
              )
            )
          ) {

            eta_sim <- rep(
              eta_hat,
              each = n_sim
            )

          } else {

            eta_se[
              is.na(eta_se) |
              eta_se <= 0
            ] <- 1e-8


            eta_sim <- unlist(
              Map(
                function(
                    mu,
                    sig
                ) {

                  stats::rnorm(
                    n_sim,
                    mean = mu,
                    sd = sig
                  )
                },

                eta_hat,
                eta_se
              )
            )
          }


          value_sim <-
            inv_link(
              eta_sim
            )


          data.frame(
            mean =
              mean(
                value_sim,
                na.rm = TRUE
              ),

            median =
              stats::median(
                value_sim,
                na.rm = TRUE
              ),

            lwr =
              stats::quantile(
                value_sim,
                0.025,
                na.rm = TRUE
              ),

            upr =
              stats::quantile(
                value_sim,
                0.975,
                na.rm = TRUE
              )
          )
        }
      ) |>
      dplyr::ungroup()


    names(out)[
      names(out) == "mean"
    ] <-
      paste0(
        prefix,
        "_mean"
      )


    names(out)[
      names(out) == "median"
    ] <-
      paste0(
        prefix,
        "_median"
      )


    names(out)[
      names(out) == "lwr"
    ] <-
      paste0(
        prefix,
        "_lwr"
      )


    names(out)[
      names(out) == "upr"
    ] <-
      paste0(
        prefix,
        "_upr"
      )


    out
  }


  # ============================================================
  # Retained post-burn-in iterations
  # ============================================================

  keep <- seq.int(
    burn_in + 1,
    n_iter
  )


  # ============================================================
  # NEW:
  # Combine fixed-effect beta estimates after burn-in
  # ============================================================

  beta_occ_draws <-
    dplyr::bind_rows(
      beta_occ_list[
        keep
      ]
    )


  beta_capture_draws <-
    dplyr::bind_rows(
      beta_capture_list[
        keep
      ]
    )


  beta_abundance_draws <-
    dplyr::bind_rows(
      beta_abundance_list[
        keep
      ]
    )


  beta_zi_draws <-
    dplyr::bind_rows(
      beta_zi_list[
        keep
      ]
    )


  beta_all_draws <-
    dplyr::bind_rows(
      beta_occ_draws,
      beta_capture_draws,
      beta_abundance_draws,
      beta_zi_draws
    )


  # ============================================================
  # NEW:
  # Pool beta estimates and their variances
  # ============================================================

  beta_pooled <-
    pool_beta_MI(
      beta_all_draws
    )


  # ============================================================
  # NEW:
  # Natural-scale transformation for INTERCEPTS
  #
  # Occupancy:
  #   beta -> psi
  #
  # Capture:
  #   beta -> p
  #
  # Abundance:
  #   beta -> lambda
  #
  # Zero inflation:
  #   beta -> pi
  #
  # Slopes remain on coefficient scale because their natural-scale
  # interpretation depends on the covariate configuration.
  # ============================================================

  if (
    nrow(
      beta_pooled
    ) > 0
  ) {

    beta_pooled_natural <-
      beta_pooled |>
      dplyr::mutate(

        natural_estimate =
          dplyr::case_when(

            .data$component %in%
              c(
                "occupancy",
                "capture",
                "zero_inflation"
              ) &
              .data$term ==
                "(Intercept)" ~

              stats::plogis(
                .data$estimate
              ),


            .data$component ==
              "abundance" &
              .data$term ==
                "(Intercept)" ~

              exp(
                .data$estimate
              ),


            TRUE ~
              NA_real_
          ),


        natural_lower =
          dplyr::case_when(

            .data$component %in%
              c(
                "occupancy",
                "capture",
                "zero_inflation"
              ) &
              .data$term ==
                "(Intercept)" ~

              stats::plogis(
                .data$lower
              ),


            .data$component ==
              "abundance" &
              .data$term ==
                "(Intercept)" ~

              exp(
                .data$lower
              ),


            TRUE ~
              NA_real_
          ),


        natural_upper =
          dplyr::case_when(

            .data$component %in%
              c(
                "occupancy",
                "capture",
                "zero_inflation"
              ) &
              .data$term ==
                "(Intercept)" ~

              stats::plogis(
                .data$upper
              ),


            .data$component ==
              "abundance" &
              .data$term ==
                "(Intercept)" ~

              exp(
                .data$upper
              ),


            TRUE ~
              NA_real_
          )
      )

  } else {

    beta_pooled_natural <-
      data.frame()
  }


  # ============================================================
  # Final models
  # ============================================================

  # ============================================================
# Final models
# ============================================================

final_iter <- n_iter

occ_fit_final <- occupancy_models[[final_iter]]

cap_fit_final <- capture_models[[final_iter]]

abund_fit_final <- abundance_models[[final_iter]]


  # ============================================================
  # Return
  # ============================================================

  list(

    # ----------------------------------------------------------
    # Prediction-level summaries
    # ----------------------------------------------------------

    psi =
      summarise_link(
        psi_list[
          keep
        ],
        "logit",
        "psi"
      ),


    capture =
      summarise_link(
        capture_list[
          keep
        ],
        "logit",
        "capture"
      ),


    lambda =
      summarise_link(
        lambda_list[
          keep
        ],
        "log",
        "lambda"
      ),


    p_detect =
      summarise_link(
        p_detect_list[
          keep
        ],
        "cloglog",
        "p_detect"
      ),


    # ----------------------------------------------------------
    # NEW:
    # Properly pooled fixed-effect beta estimates
    # ----------------------------------------------------------

    beta =
      beta_pooled,


    beta_natural =
      beta_pooled_natural,


    beta_iterations =
      list(

        occupancy =
          beta_occ_draws,

        capture =
          beta_capture_draws,

        abundance =
          beta_abundance_draws,

        zero_inflation =
          beta_zi_draws
      ),


    # ----------------------------------------------------------
    # Existing prediction-level retained iterations
    # ----------------------------------------------------------

    psi_list =
      psi_list[
        keep
      ],


    capture_list =
      capture_list[
        keep
      ],


    lambda_list =
      lambda_list[
        keep
      ],


    p_detect_list =
      p_detect_list[
        keep
      ],


    # ----------------------------------------------------------
    # Final completed-data GLMM fits
    # ----------------------------------------------------------

    occ_fit =
      occ_fit_final,


    cap_fit =
      cap_fit_final,


    abund_fit =
      abund_fit_final,


    # ----------------------------------------------------------
    # Family
    # ----------------------------------------------------------

    abundance_family =
      abundance_family,


    # ----------------------------------------------------------
    # All retained GLMM fits
    # ----------------------------------------------------------

    occupancy_models =
      occupancy_models[
        keep
      ],


    capture_models =
      capture_models[
        keep
      ],


    abundance_models =
      abundance_models[
        keep
      ],


    # ----------------------------------------------------------
    # Data
    # ----------------------------------------------------------

    site_data =
      site_data,


    sample_data =
      sample_data,


    long_df =
      long_df,


    # ----------------------------------------------------------
    # Filtering
    # ----------------------------------------------------------

    filter_summary =
      list(

        otu_stats =
          otu_stats,

        kept_otus =
          keep_otus
      ),


    # ----------------------------------------------------------
    # Diagnostic component AIC
    # ----------------------------------------------------------

    diagnostic_AIC =
      diagnostic_AIC,


    # ----------------------------------------------------------
    # Method note
    # ----------------------------------------------------------

    note =
      paste(
        "Approximate three-level stochastic-imputation GLMM.",
        "Hierarchy: site occupancy Z -> biological-sample capture A -> PCR replicate counts Y.",
        "Fixed-effect coefficients are pooled across retained iterations.",
        "Their uncertainty combines average within-GLMM variance and",
        "between-iteration latent-state variability using Rubin-style pooling.",
        "Prediction summaries continue to integrate eta uncertainty by simulation.",
        "AIC values are exploratory component-model diagnostics, not a formal joint-likelihood criterion."
      )
  )
}
