#' Fit a Hierarchical Occupancy–Detection–Abundance Model for Microbiome Data
#'
#' Fits a three-level hierarchical model for eDNA / microbiome count data using an
#' EM-like iterative procedure based on repeated GLMM estimation via \code{glmmTMB}.
#'
#' The model explicitly represents the ecological observation process:
#'
#' \itemize{
#'   \item \strong{Site level}: true presence/absence (\eqn{Z})
#'   \item \strong{Biological sample level}: capture/detection (\eqn{A})
#'   \item \strong{PCR replicate level}: observed read counts (\eqn{Y})
#' }
#'
#' Covariates specified in the model formulas are automatically preserved and
#' propagated into the internally constructed datasets (site-level and sample-level),
#' allowing flexible inclusion of environmental variables, treatments, interactions,
#' and random effects.
#'
#' @param phyloseq A \code{phyloseq} object containing OTU counts and sample metadata.
#' @param site_col Character. Column identifying sites (occupancy level).
#' @param sample_col Character. Column identifying biological samples (capture level).
#' @param replicate_col Character or \code{NULL}. Column identifying PCR replicates.
#' @param occupancy_formula Formula for occupancy model (response must be \code{z_sim}).
#' @param capture_formula Formula for capture model (response must be \code{a_sim}).
#' @param abundance_formula Formula for abundance model (response must match \code{count_col}).
#' @param otu_col Character. OTU identifier column.
#' @param count_col Character. Count column (e.g. read counts).
#' @param min_species_sum Minimum total counts required to retain an OTU.
#' @param min_detection_replicates Minimum number of detections required per OTU.
#' @param abundance_threshold Threshold defining detection (default = 0).
#' @param n_iter Number of EM-like iterations.
#' @param burn_in Number of initial iterations discarded.
#' @param abundance_family One of \code{"poisson"}, \code{"nbinom"}, \code{"zip"}, \code{"zinb"}.
#' @param verbose Logical; print progress during fitting.
#' @param n_sim_ci Number of simulations used for uncertainty propagation (ψ, capture, λ only).
#'
#' @return A list containing:
#' \describe{
#'   \item{psi}{Occupancy probability summaries (\code{psi_mean}, \code{psi_lwr}, etc.).}
#'   \item{capture}{Capture probability summaries.}
#'   \item{lambda}{Abundance (expected counts) summaries.}
#'   \item{p_detect}{Detection probability derived from abundance (natural scale).}
#'   \item{psi_list, capture_list, lambda_list}{Per-iteration linear predictors (eta) and SE.}
#'   \item{p_detect_list}{Per-iteration detection probabilities (natural scale, no SE).}
#'   \item{occupancy_models, capture_models, abundance_models}{Fitted GLMM objects.}
#'   \item{site_data, sample_data, long_df}{Processed datasets.}
#'   \item{filter_summary}{OTU filtering information.}
#'   \item{diagnostic_AIC}{Component-wise AIC per iteration.}
#' }
#'
#' @details
#' The hierarchical model is defined as:
#'
#' \deqn{
#' Z_{i} \sim \text{Bernoulli}(\psi_i)
#' }
#' \deqn{
#' A_{ij} \mid Z_i \sim \text{Bernoulli}(p_{ij} \cdot Z_i)
#' }
#' \deqn{
#' Y_{ijk} \mid A_{ij} \sim \text{Count}(\lambda_{ijk} \cdot A_{ij})
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{i} indexes sites
#'   \item \eqn{j} indexes biological samples
#'   \item \eqn{k} indexes PCR replicates
#' }
#'
#' The abundance model defines \eqn{\lambda}, which induces a detection probability:
#'
#' \deqn{
#' p_{\text{detect}} = 1 - P(Y = 0 \mid A = 1)
#' }
#'
#' For a Poisson model:
#'
#' \deqn{
#' p_{\text{detect}} = 1 - \exp(-\lambda)
#' }
#'
#' For NB / ZIP / ZINB models, the zero probability is computed using the
#' corresponding distribution (including zero-inflation when applicable).
#'
#' \strong{Important:}
#' \itemize{
#'   \item \code{p_detect} is computed directly from the abundance model.
#'   \item It is not an independently estimated parameter.
#'   \item It is stored on the natural scale (not as a linear predictor).
#' }
#'
#' @section EM-like algorithm:
#' Each iteration performs:
#' \enumerate{
#'   \item Fit occupancy model using current \eqn{Z}
#'   \item Fit capture model conditional on \eqn{Z = 1}
#'   \item Fit abundance model conditional on \eqn{A = 1}
#'   \item Compute detection probability from abundance model
#'   \item Update \eqn{A} using capture + abundance probabilities
#'   \item Update \eqn{Z} using capture + abundance probabilities
#' }
#'
#' @section Uncertainty estimation:
#' Uncertainty is propagated via simulation for:
#'
#' \itemize{
#'   \item Occupancy (\eqn{\psi})
#'   \item Capture probability
#'   \item Abundance (\eqn{\lambda})
#' }
#'
#' using:
#'
#' \deqn{\eta \sim \mathcal{N}(\hat{\eta}, \text{SE}^2)}
#'
#' followed by transformation to the natural scale.
#'
#' Detection probability (\code{p_detect}) is summarized empirically across iterations,
#' as it is a deterministic transformation of the abundance model.
#'
#' @section Model features:
#' \itemize{
#'   \item Explicit 3-level eDNA hierarchy
#'   \item Supports random effects via \code{glmmTMB}
#'   \item Handles Poisson, NB, ZIP, ZINB families
#'   \item Supports offsets (e.g. \code{offset(log(total_reads))})
#'   \item Automatically propagates covariates across hierarchy levels
#' }
#'
#' @section Caveats:
#' \itemize{
#'   \item Approximate EM-like method (not a full joint likelihood)
#'   \item AIC values are component-wise diagnostics only
#'   \item Detection depends on abundance distribution assumptions
#'   \item Large \eqn{\lambda} values imply detection saturation (\eqn{p_{\text{detect}} \approx 1})
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
#'
#'   occupancy_formula = z_sim ~ 1 + (1 | OTU),
#'   capture_formula   = a_sim ~ 1 + (1 | OTU),
#'   abundance_formula = y ~ offset(log(total_reads)) + (1 | OTU),
#'
#'   abundance_family = "nbinom",
#'   n_iter = 20,
#'   burn_in = 5
#' )
#'
#' head(fit$psi)
#' head(fit$capture)
#' head(fit$lambda)
#' head(fit$p_detect)
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
    n_iter = 50,
    burn_in = 10,
    abundance_family = c("poisson", "nbinom", "zip", "zinb"),
    verbose = TRUE,
    n_sim_ci = 1000
) {

  abundance_family <- match.arg(abundance_family)

  if (is.null(otu_col)) stop("Please specify otu_col.")
  if (is.null(count_col)) stop("Please specify count_col.")
  if (burn_in >= n_iter) stop("burn_in must be < n_iter.")

  get_formula_vars <- function(formula, response) {
    setdiff(all.vars(formula), response)
  }

  # ------------------------------------------------------------
  # Prepare long data
  # ------------------------------------------------------------

  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col,
    nested_cols = unique(c(sample_col, replicate_col))
  )

  long_df <- prep$long_df

  needed_cols <- c(site_col, sample_col, otu_col, count_col)
  missing_cols <- setdiff(needed_cols, names(long_df))

  if (length(missing_cols) > 0) {
    stop("Missing columns in long_df: ", paste(missing_cols, collapse = ", "))
  }

  # ------------------------------------------------------------
  # Add sample metadata
  # ------------------------------------------------------------

  sample_meta <- data.frame(
    phyloseq::sample_data(phyloseq),
    check.names = FALSE
  )

  sample_meta$.__sample_id__ <- rownames(sample_meta)

  if (!(sample_col %in% names(sample_meta))) {
    sample_meta[[sample_col]] <- sample_meta$.__sample_id__
  }

  sample_meta[[sample_col]] <- as.character(sample_meta[[sample_col]])
  long_df[[sample_col]] <- as.character(long_df[[sample_col]])

  meta_cols_to_add <- setdiff(names(sample_meta), names(long_df))
  meta_cols_to_add <- setdiff(meta_cols_to_add, ".__sample_id__")

  if (length(meta_cols_to_add) > 0) {
    long_df <- dplyr::left_join(
      long_df,
      sample_meta[, c(sample_col, meta_cols_to_add), drop = FALSE],
      by = sample_col
    )
  }

  # ------------------------------------------------------------
  # Type handling
  # ------------------------------------------------------------

  long_df[[site_col]]   <- as.character(long_df[[site_col]])
  long_df[[otu_col]]    <- as.character(long_df[[otu_col]])
  long_df[[sample_col]] <- as.character(long_df[[sample_col]])
  long_df[[count_col]]  <- as.numeric(long_df[[count_col]])

  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    long_df[[replicate_col]] <- as.character(long_df[[replicate_col]])
  }

  long_df <- long_df[!is.na(long_df[[count_col]]), ]

  # ------------------------------------------------------------
  # Compute total_reads
  # ------------------------------------------------------------

  total_reads_df <- long_df |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      total_reads = sum(.data[[count_col]], na.rm = TRUE),
      .groups = "drop"
    )

  long_df <- dplyr::left_join(long_df, total_reads_df, by = sample_col)

  uses_offset <- any(grepl("offset", deparse(abundance_formula)))

  if (uses_offset) {
    if (verbose) message("Offset detected -> cleaning zero-read samples...")

    zero_samples <- unique(long_df[[sample_col]][long_df$total_reads <= 0])

    if (length(zero_samples) > 0) {
      long_df <- long_df |>
        dplyr::filter(!(.data[[sample_col]] %in% zero_samples))

      if (verbose) {
        message(length(zero_samples), " zero-read samples removed.")
      }
    }

    if (any(long_df$total_reads <= 0, na.rm = TRUE)) {
      stop("total_reads must be positive when using offset(log(total_reads)).")
    }
  } else {
    if (verbose) message("No offset used.")
  }

  # ------------------------------------------------------------
  # Check formula variables
  # ------------------------------------------------------------

  occ_vars <- get_formula_vars(occupancy_formula, "z_sim")
  cap_vars <- get_formula_vars(capture_formula, "a_sim")
  abund_vars <- get_formula_vars(abundance_formula, count_col)

  all_formula_vars <- unique(c(occ_vars, cap_vars, abund_vars))
  all_formula_vars <- setdiff(all_formula_vars, c(count_col, "total_reads"))

  missing_formula_vars <- setdiff(all_formula_vars, names(long_df))

  if (length(missing_formula_vars) > 0) {
    stop(
      "The following formula variables are missing from long_df: ",
      paste(missing_formula_vars, collapse = ", ")
    )
  }

  for (v in all_formula_vars) {
    if (v %in% names(long_df) && is.character(long_df[[v]])) {
      long_df[[v]] <- factor(long_df[[v]])
    }
  }

  # ------------------------------------------------------------
  # OTU filtering
  # ------------------------------------------------------------

  otu_stats <- long_df |>
    dplyr::group_by(.data[[otu_col]]) |>
    dplyr::summarise(
      total_count = sum(.data[[count_col]], na.rm = TRUE),
      detected_replicates = sum(.data[[count_col]] > abundance_threshold, na.rm = TRUE),
      .groups = "drop"
    )

  keep_otus <- otu_stats |>
    dplyr::filter(
      .data$total_count >= min_species_sum,
      .data$detected_replicates >= min_detection_replicates
    ) |>
    dplyr::pull(.data[[otu_col]])

  long_df <- long_df |>
    dplyr::filter(.data[[otu_col]] %in% keep_otus)

  if (nrow(long_df) == 0) {
    stop("No OTUs remain after filtering.")
  }

  long_df[[site_col]]   <- factor(long_df[[site_col]])
  long_df[[otu_col]]    <- factor(long_df[[otu_col]])
  long_df[[sample_col]] <- factor(long_df[[sample_col]])

  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    long_df[[replicate_col]] <- factor(long_df[[replicate_col]])
  }

  if (verbose) {
    message("OTUs before filtering: ", dplyr::n_distinct(otu_stats[[otu_col]]))
    message("OTUs after filtering: ", dplyr::n_distinct(long_df[[otu_col]]))
  }

  # ------------------------------------------------------------
  # Define hierarchy keys
  # ------------------------------------------------------------

  site_keys <- c(site_col, otu_col)
  sample_keys <- c(site_col, sample_col, otu_col)

  pcr_keys <- c(site_col, sample_col, otu_col)

  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    pcr_keys <- c(site_col, sample_col, replicate_col, otu_col)
  }

  # ------------------------------------------------------------
  # Build Z and A data
  # ------------------------------------------------------------

  site_keep_vars <- setdiff(occ_vars, c(site_keys, "z_sim"))
  site_keep_vars <- intersect(site_keep_vars, names(long_df))

  sample_keep_vars <- setdiff(cap_vars, c(sample_keys, "a_sim"))
  sample_keep_vars <- intersect(sample_keep_vars, names(long_df))

  site_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
    dplyr::summarise(
      z_obs = as.integer(any(.data[[count_col]] > abundance_threshold)),
      dplyr::across(
        dplyr::all_of(site_keep_vars),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(z_sim = .data$z_obs)

  sample_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
    dplyr::summarise(
      a_obs = as.integer(any(.data[[count_col]] > abundance_threshold)),
      dplyr::across(
        dplyr::all_of(sample_keep_vars),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(a_sim = .data$a_obs)

  # ------------------------------------------------------------
  # Validate formulas
  # ------------------------------------------------------------

  validate_formula <- function(formula, data, expected_response, model_name) {

    if (!inherits(formula, "formula")) {
      stop(model_name, " formula must be a valid formula.")
    }

    response <- all.vars(formula)[1]

    if (response != expected_response) {
      stop(
        model_name, " formula response must be '", expected_response,
        "', but got '", response, "'."
      )
    }

    vars <- all.vars(formula)
    missing_vars <- setdiff(vars, names(data))

    if (length(missing_vars) > 0) {
      stop(
        "Missing variables in ", model_name, " formula: ",
        paste(missing_vars, collapse = ", ")
      )
    }

    invisible(TRUE)
  }

  validate_formula(occupancy_formula, site_data, "z_sim", "occupancy")
  validate_formula(capture_formula, sample_data, "a_sim", "capture")
  validate_formula(abundance_formula, long_df, count_col, "abundance")

  # ------------------------------------------------------------
  # Family setup
  # ------------------------------------------------------------

  fam <- switch(
    abundance_family,
    poisson = stats::poisson(),
    nbinom  = glmmTMB::nbinom2(),
    zip     = stats::poisson(),
    zinb    = glmmTMB::nbinom2()
  )

  zi_formula <- if (abundance_family %in% c("zip", "zinb")) ~1 else ~0

  # ------------------------------------------------------------
  # Helper: PCR-level zero probability
  # ------------------------------------------------------------

  get_p0_pcr <- function(fit, newdata, family_name) {

    mu <- as.numeric(predict(
      fit,
      type = "response",
      newdata = newdata,
      allow.new.levels = TRUE
    ))

    zi <- rep(0, length(mu))

    if (family_name %in% c("zip", "zinb")) {
      zi <- tryCatch(
        as.numeric(predict(
          fit,
          type = "zprob",
          newdata = newdata,
          allow.new.levels = TRUE
        )),
        error = function(e) rep(0, length(mu))
      )
    }

    p0_cond <- switch(
      family_name,
      poisson = stats::dpois(0, lambda = mu),
      zip     = stats::dpois(0, lambda = mu),
      nbinom  = {
        theta <- tryCatch(stats::sigma(fit), error = function(e) NA_real_)
        if (!is.finite(theta) || theta <= 0) {
          stats::dpois(0, lambda = mu)
        } else {
          stats::dnbinom(0, mu = mu, size = theta)
        }
      },
      zinb = {
        theta <- tryCatch(stats::sigma(fit), error = function(e) NA_real_)
        if (!is.finite(theta) || theta <= 0) {
          stats::dpois(0, lambda = mu)
        } else {
          stats::dnbinom(0, mu = mu, size = theta)
        }
      }
    )

    p0 <- zi + (1 - zi) * p0_cond
    pmin(pmax(p0, 1e-12), 1)
  }

  # ------------------------------------------------------------
  # Storage
  # ------------------------------------------------------------

  psi_list <- list()
  capture_list <- list()
  lambda_list <- list()
  p_detect_list <- list()

  occupancy_models <- list()
  capture_models <- list()
  abundance_models <- list()

  diagnostic_AIC <- data.frame(
    iteration = integer(),
    occupancy_AIC = numeric(),
    capture_AIC = numeric(),
    abundance_AIC = numeric()
  )

  # ------------------------------------------------------------
  # EM-like iterations
  # ------------------------------------------------------------

  for (i in seq_len(n_iter)) {

    if (verbose) message("Iteration ", i)

    # ----------------------------------------------------------
    # 1. Occupancy model
    # ----------------------------------------------------------

    occ_fit <- glmmTMB::glmmTMB(
      formula = occupancy_formula,
      data = site_data,
      family = stats::binomial()
    )

    pred_psi <- predict(
      occ_fit,
      type = "link",
      newdata = site_data,
      allow.new.levels = TRUE,
      se.fit = TRUE
    )

    eta_psi <- as.numeric(pred_psi$fit)
    psi <- stats::plogis(eta_psi)

    psi_list[[i]] <- data.frame(
      site_data[site_keys],
      eta = eta_psi,
      se = as.numeric(pred_psi$se.fit)
    )

    occupancy_models[[i]] <- occ_fit

    # ----------------------------------------------------------
    # 2. Capture model
    # ----------------------------------------------------------

    sample_data <- sample_data |>
      dplyr::select(-dplyr::any_of(c("z_sim", "p0_sample", "capture_prob"))) |>
      dplyr::left_join(
        site_data |>
          dplyr::select(dplyr::all_of(c(site_keys, "z_sim"))),
        by = site_keys
      )

    capture_fit_data <- sample_data |>
      dplyr::filter(.data$z_sim == 1)

    if (nrow(capture_fit_data) == 0) {
      stop("No rows available for capture model at iteration ", i)
    }

    cap_fit <- glmmTMB::glmmTMB(
      formula = capture_formula,
      data = capture_fit_data,
      family = stats::binomial()
    )

    pred_capture <- predict(
      cap_fit,
      type = "link",
      newdata = sample_data,
      allow.new.levels = TRUE,
      se.fit = TRUE
    )

    eta_capture <- as.numeric(pred_capture$fit)
    capture_prob <- stats::plogis(eta_capture)

    capture_list[[i]] <- data.frame(
      sample_data[sample_keys],
      eta = eta_capture,
      se = as.numeric(pred_capture$se.fit)
    )

    capture_models[[i]] <- cap_fit

    # ----------------------------------------------------------
    # 3. Abundance model
    # ----------------------------------------------------------

    abund_fit_data <- long_df |>
      dplyr::left_join(
        sample_data |>
          dplyr::select(dplyr::all_of(c(sample_keys, "a_sim"))),
        by = sample_keys
      ) |>
      dplyr::filter(.data$a_sim == 1)

    if (nrow(abund_fit_data) == 0) {
      stop("No rows available for abundance model at iteration ", i)
    }

    abund_fit <- glmmTMB::glmmTMB(
      formula = abundance_formula,
      data = abund_fit_data,
      family = fam,
      ziformula = zi_formula
    )

    pred_lambda <- predict(
      abund_fit,
      type = "link",
      newdata = long_df,
      allow.new.levels = TRUE,
      se.fit = TRUE
    )

    eta_lambda_all <- as.numeric(pred_lambda$fit)
    lambda_all <- exp(eta_lambda_all)

    lambda_list[[i]] <- data.frame(
      long_df[pcr_keys],
      eta = eta_lambda_all,
      se = as.numeric(pred_lambda$se.fit)
    )

    p0_pcr <- get_p0_pcr(
      fit = abund_fit,
      newdata = long_df,
      family_name = abundance_family
    )

    p_detect_eta <- log(-log(pmax(p0_pcr, 1e-12)))

    p_detect_list[[i]] <- data.frame(
      long_df[pcr_keys],
      eta = p_detect_eta,
      se = NA_real_
    )

    abundance_models[[i]] <- abund_fit

    # ----------------------------------------------------------
    # 4. Collapse PCR probabilities to sample level
    # ----------------------------------------------------------

    sample_p0 <- long_df |>
      dplyr::mutate(p0_pcr = p0_pcr) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
      dplyr::summarise(
        p0_sample = prod(.data$p0_pcr, na.rm = TRUE),
        .groups = "drop"
      )

    sample_data <- sample_data |>
      dplyr::select(-dplyr::any_of(c("p0_sample", "capture_prob"))) |>
      dplyr::mutate(capture_prob = capture_prob) |>
      dplyr::left_join(sample_p0, by = sample_keys)

    sample_data$p0_sample[is.na(sample_data$p0_sample)] <- 1

    # ----------------------------------------------------------
    # 5. Update A
    # ----------------------------------------------------------

    zero_sample <- which(sample_data$a_obs == 0 & sample_data$z_sim == 1)

    if (length(zero_sample) > 0) {

      numerator_a <- sample_data$capture_prob[zero_sample] *
        sample_data$p0_sample[zero_sample]

      denominator_a <- (1 - sample_data$capture_prob[zero_sample]) +
        sample_data$capture_prob[zero_sample] *
        sample_data$p0_sample[zero_sample]

      posterior_a <- numerator_a / pmax(denominator_a, 1e-12)
      posterior_a <- pmin(pmax(posterior_a, 0.001), 0.999)

      sample_data$a_sim[zero_sample] <- stats::rbinom(
        n = length(zero_sample),
        size = 1,
        prob = posterior_a
      )
    }

    sample_data$a_sim[sample_data$a_obs == 1] <- 1
    sample_data$a_sim[sample_data$z_sim == 0] <- 0

    # ----------------------------------------------------------
    # 6. Update Z
    # ----------------------------------------------------------

    site_p0 <- sample_data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(
        p0_site = prod(
          (1 - .data$capture_prob) +
            .data$capture_prob * .data$p0_sample,
          na.rm = TRUE
        ),
        .groups = "drop"
      )

    site_data <- site_data |>
      dplyr::select(-dplyr::any_of("p0_site")) |>
      dplyr::left_join(site_p0, by = site_keys)

    site_data$p0_site[is.na(site_data$p0_site)] <- 1

    zero_site <- which(site_data$z_obs == 0)

    if (length(zero_site) > 0) {

      numerator_z <- psi[zero_site] *
        site_data$p0_site[zero_site]

      denominator_z <- (1 - psi[zero_site]) +
        psi[zero_site] * site_data$p0_site[zero_site]

      posterior_z <- numerator_z / pmax(denominator_z, 1e-12)
      posterior_z <- pmin(pmax(posterior_z, 0.001), 0.999)

      site_data$z_sim[zero_site] <- stats::rbinom(
        n = length(zero_site),
        size = 1,
        prob = posterior_z
      )
    }

    site_data$z_sim[site_data$z_obs == 1] <- 1

    diagnostic_AIC <- rbind(
      diagnostic_AIC,
      data.frame(
        iteration = i,
        occupancy_AIC = stats::AIC(occ_fit),
        capture_AIC = stats::AIC(cap_fit),
        abundance_AIC = stats::AIC(abund_fit)
      )
    )
  }

  # ------------------------------------------------------------
  # Summaries on natural scale with eta uncertainty
  # ------------------------------------------------------------

  summarise_link <- function(lst, link_name, prefix, n_sim = n_sim_ci) {

    inv_link <- switch(
      link_name,
      logit   = stats::plogis,
      log     = exp,
      cloglog = function(x) 1 - exp(-exp(x)),
      stop("Unknown link_name: ", link_name)
    )

    df <- dplyr::bind_rows(lst)

    if (nrow(df) == 0) return(data.frame())

    if (!("eta" %in% names(df))) {
      stop("Each list element must contain an eta column.")
    }

    if (!("se" %in% names(df))) {
      df$se <- NA_real_
    }

    keys <- setdiff(names(df), c("eta", "se"))

    out <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::group_modify(function(.x, .y) {

        eta_hat <- .x$eta
        eta_se  <- .x$se

        if (all(is.na(eta_se))) {
          eta_sim <- rep(eta_hat, each = n_sim)
        } else {
          eta_se[is.na(eta_se) | eta_se <= 0] <- 1e-8

          eta_sim <- unlist(
            Map(
              function(mu, sig) {
                stats::rnorm(n_sim, mean = mu, sd = sig)
              },
              eta_hat,
              eta_se
            )
          )
        }

        value_sim <- inv_link(eta_sim)

        data.frame(
          mean   = mean(value_sim, na.rm = TRUE),
          median = stats::median(value_sim, na.rm = TRUE),
          lwr    = stats::quantile(value_sim, 0.025, na.rm = TRUE),
          upr    = stats::quantile(value_sim, 0.975, na.rm = TRUE)
        )
      }) |>
      dplyr::ungroup()

    names(out)[names(out) == "mean"]   <- paste0(prefix, "_mean")
    names(out)[names(out) == "median"] <- paste0(prefix, "_median")
    names(out)[names(out) == "lwr"]    <- paste0(prefix, "_lwr")
    names(out)[names(out) == "upr"]    <- paste0(prefix, "_upr")

    out
  }

  keep <- seq.int(burn_in + 1, n_iter)

  final_iter <- n_iter

  occ_fit_final   <- occupancy_models[[final_iter]]
  cap_fit_final   <- capture_models[[final_iter]]
  abund_fit_final <- abundance_models[[final_iter]]

  list(
    psi = summarise_link(psi_list[keep], "logit", "psi"),
    capture = summarise_link(capture_list[keep], "logit", "capture"),
    lambda = summarise_link(lambda_list[keep], "log", "lambda"),
    p_detect = summarise_link(p_detect_list[keep], "cloglog", "p_detect"),

    psi_list = psi_list[keep],
    capture_list = capture_list[keep],
    lambda_list = lambda_list[keep],
    p_detect_list = p_detect_list[keep],

    occ_fit = occ_fit_final,
    cap_fit = cap_fit_final,
    abund_fit = abund_fit_final,

    abundance_family = abundance_family,

    occupancy_models = occupancy_models[keep],
    capture_models = capture_models[keep],
    abundance_models = abundance_models[keep],

    site_data = site_data,
    sample_data = sample_data,
    long_df = long_df,

    filter_summary = list(
      otu_stats = otu_stats,
      kept_otus = keep_otus
    ),

    diagnostic_AIC = diagnostic_AIC,

    note = paste(
      "Approximate three-level EM-like GLMM.",
      "Hierarchy: site occupancy Z -> biological-sample capture A -> PCR replicate counts Y.",
      "Intervals are based on simulation from eta ~ Normal(eta_hat, se_eta^2), then transformed to the natural scale.",
      "AIC values are exploratory component-model diagnostics, not a formal joint-likelihood criterion."
    )
  )
}
