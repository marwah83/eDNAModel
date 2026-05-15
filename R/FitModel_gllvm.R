#' Fit a Multispecies 3-Level Occupancy–Detection–Abundance Model using GLLVM and GLMM
#'
#' This function fits a hierarchical three-level model for eDNA / microbiome data
#' stored as a `phyloseq` object. The model explicitly separates:
#'
#' \itemize{
#'   \item \strong{Site level (Z)}: true OTU presence/absence
#'   \item \strong{Sample level (A)}: capture / detection process
#'   \item \strong{PCR replicate level (Y)}: observed read counts
#' }
#'
#' Occupancy is modeled at the site-by-OTU level using a binomial Generalized
#' Linear Latent Variable Model (GLLVM), while capture and abundance are modeled
#' using GLMMs via `glmmTMB`.
#'
#' @param phyloseq A `phyloseq` object containing OTU count data and sample metadata.
#' @param site_col Character string identifying sites (occupancy level).
#' @param sample_col Character string identifying biological samples (capture level).
#' @param replicate_col Character or `NULL` identifying PCR replicates (abundance level).
#' @param abundance_rhs Right-hand side of the abundance model, e.g.
#'   `(1 | OTU)` or `offset(log(total_reads)) + (1 | OTU) + (1 | Name / OTU)`.
#' @param capture_formula Formula for capture model (response must be `a_sim`).
#' @param occupancy_covars Optional character vector of site-level covariates used in the GLLVM occupancy model.
#' @param min_species_sum Minimum total read count required to retain an OTU.
#' @param min_detection_replicates Minimum number of detections required per OTU.
#' @param abundance_threshold Threshold defining detection (default = 1).
#' @param n_iter Number of EM-like iterations.
#' @param burn_in Number of initial iterations discarded.
#' @param abundance_family One of `"poisson"`, `"nbinom"`, `"zip"`, `"zinb"`.
#' @param num_lv_c Number of latent variables in the GLLVM occupancy model.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{Posterior summaries for occupancy (`psi_*`), capture (`capture_*`),
#'   abundance (`lambda_*`), and detection (`p_detect_*`).}
#'   \item{psi_list}{Per-iteration occupancy linear predictors.}
#'   \item{capture_list}{Per-iteration capture linear predictors.}
#'   \item{lambda_list}{Per-iteration abundance linear predictors.}
#'   \item{p_detect_list}{Per-iteration detection linear predictors (no SE; deterministic from abundance).}
#'   \item{occupancy_models}{List of fitted GLLVM occupancy models.}
#'   \item{capture_models}{List of fitted GLMM capture models.}
#'   \item{abundance_models}{List of fitted GLMM abundance models.}
#'   \item{site_data, sample_data, long_df}{Processed hierarchical datasets.}
#'   \item{lv_sites, lv_species}{Latent variable coordinates.}
#'   \item{mean_lv_sites, mean_lv_species}{Posterior mean latent variables.}
#'   \item{filter_summary}{OTU filtering information.}
#' }
#'
#' @details
#' The model defines a three-level hierarchical structure:
#'
#' \deqn{
#' Z_{i,m} \sim \text{Bernoulli}(\psi_{i,m})
#' }
#'
#' \deqn{
#' A_{i,j,m} \mid Z_{i,m} \sim \text{Bernoulli}(p_{i,j,m} \cdot Z_{i,m})
#' }
#'
#' \deqn{
#' Y_{i,j,k,m} \mid A_{i,j,m} \sim \text{Count}(\lambda_{i,j,k,m} \cdot A_{i,j,m})
#' }
#'
#' where:
#' \itemize{
#'   \item \(i\): site
#'   \item \(j\): biological sample
#'   \item \(k\): PCR replicate
#'   \item \(m\): OTU
#' }
#'
#' Observed detection is defined as:
#'
#' \deqn{
#' z^{obs}_{i,m} = I\left(\max_{j,k} Y_{i,j,k,m} > c\right)
#' }
#'
#' where \(c\) is the `abundance_threshold`.
#'
#' The occupancy model uses a GLLVM:
#'
#' \deqn{
#' \text{logit}(\psi_{i,m}) = X_i \beta_m + \text{latent variables}
#' }
#'
#' The capture model is a binomial GLMM:
#'
#' \deqn{
#' \text{logit}(p_{i,j,m}) = W_{i,j} \gamma_m
#' }
#'
#' The abundance model defines \(\lambda\), which induces a detection probability:
#'
#' \deqn{
#' p_{\text{detect}} = 1 - P(Y = 0 \mid A = 1)
#' }
#'
#' This is computed using the specified count distribution (Poisson, NB, ZIP, ZINB).
#'
#' \strong{Important:}
#' \itemize{
#'   \item `p_detect` is a deterministic function of the abundance model
#'   \item No standard error is estimated for `p_detect`
#' }
#'
#' @section EM-like algorithm:
#' Each iteration performs:
#' \enumerate{
#'   \item Fit GLLVM occupancy model for \(Z\)
#'   \item Fit GLMM capture model for \(A\) conditional on \(Z=1\)
#'   \item Fit GLMM abundance model for \(Y\) conditional on \(A=1\)
#'   \item Update \(A\) using capture and abundance probabilities
#'   \item Update \(Z\) using capture and abundance probabilities
#' }
#'
#' @section Uncertainty estimation:
#' Uncertainty is propagated using simulation on the link scale:
#'
#' \deqn{
#' \eta \sim \mathcal{N}(\hat{\eta}, \text{SE}^2)
#' }
#'
#' For detection (`p_detect`), uncertainty is propagated without SE because it is
#' a deterministic transformation of the abundance model.
#'
#' @section Model features:
#' \itemize{
#'   \item Full 3-level eDNA hierarchy (Z → A → Y)
#'   \item Latent ecological structure via GLLVM
#'   \item Flexible GLMMs for capture and abundance
#'   \item Supports Poisson, NB, ZIP, ZINB
#'   \item Handles offsets (e.g. `offset(log(total_reads))`)
#' }
#'
#' @section Caveats:
#' \itemize{
#'   \item Approximate EM-like inference, not full joint likelihood
#'   \item Detection depends on abundance model assumptions
#'   \item Large \(\lambda\) implies detection saturation
#' }
#'
#' @examples
#' \dontrun{
#' fit <- FitModel_gllvm(
#'   phyloseq = ps,
#'   site_col = "site_month",
#'   sample_col = "Name",
#'   replicate_col = "Replicate",
#'
#'   abundance_rhs = (1 | OTU),
#'   capture_formula = a_sim ~ 1 + (1 | OTU),
#'
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
  abundance_rhs,
  capture_formula = a_sim ~ 1 + (1 | OTU),
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

  abundance_family <- match.arg(abundance_family)

  if (burn_in >= n_iter) {
    stop("burn_in must be smaller than n_iter.")
  }

  abundance_rhs_expr <- substitute(abundance_rhs)
  abundance_rhs_txt <- paste(deparse(abundance_rhs_expr), collapse = " ")

  abundance_formula <- stats::as.formula(
    paste("y ~", abundance_rhs_txt),
    env = parent.frame()
  )

  to_factor_cols <- function(df, cols) {
    for (col in cols) {
      if (!is.null(col) && col %in% names(df)) {
        df[[col]] <- as.factor(as.character(df[[col]]))
      }
    }
    df
  }

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

  needed_cols <- c(site_col, sample_col, "OTU", "y")
  missing_cols <- setdiff(needed_cols, names(long_df))

  if (length(missing_cols) > 0) {
    stop("Missing columns in long_df: ", paste(missing_cols, collapse = ", "))
  }

  # ------------------------------------------------------------
  # Type handling
  # ------------------------------------------------------------

  long_df[[site_col]]   <- as.character(long_df[[site_col]])
  long_df[[sample_col]] <- as.character(long_df[[sample_col]])
  long_df$OTU           <- as.character(long_df$OTU)
  long_df$y             <- as.numeric(long_df$y)

  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    long_df[[replicate_col]] <- as.character(long_df[[replicate_col]])
  }

  long_df <- long_df[!is.na(long_df$y), ]

  # ------------------------------------------------------------
  # Compute total_reads for optional offset(log(total_reads))
  # ------------------------------------------------------------

  total_reads_df <- long_df |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      total_reads = sum(.data$y, na.rm = TRUE),
      .groups = "drop"
    )

  long_df <- dplyr::left_join(long_df, total_reads_df, by = sample_col)

  uses_offset <- any(grepl("offset", deparse(abundance_formula)))

  if (uses_offset) {
    zero_samples <- unique(long_df[[sample_col]][long_df$total_reads <= 0])

    if (length(zero_samples) > 0) {
      long_df <- long_df |>
        dplyr::filter(!(.data[[sample_col]] %in% zero_samples))
    }

    if (any(long_df$total_reads <= 0, na.rm = TRUE)) {
      stop("total_reads must be positive when using offset(log(total_reads)).")
    }
  }

  # ------------------------------------------------------------
  # OTU filtering
  # ------------------------------------------------------------

  otu_stats <- long_df |>
    dplyr::group_by(.data$OTU) |>
    dplyr::summarise(
      total_count = sum(.data$y, na.rm = TRUE),
      detected_replicates = sum(.data$y > abundance_threshold, na.rm = TRUE),
      .groups = "drop"
    )

  keep_otus <- otu_stats |>
    dplyr::filter(
      .data$total_count >= min_species_sum,
      .data$detected_replicates >= min_detection_replicates
    ) |>
    dplyr::pull(.data$OTU)

  long_df <- long_df |>
    dplyr::filter(.data$OTU %in% keep_otus)

  if (nrow(long_df) == 0) {
    stop("No OTUs remain after filtering.")
  }

  if (verbose) {
    message("OTUs before filtering: ", dplyr::n_distinct(otu_stats$OTU))
    message("OTUs after filtering: ", dplyr::n_distinct(long_df$OTU))
  }

  # ------------------------------------------------------------
  # Factor handling
  # ------------------------------------------------------------

  factor_cols <- unique(c(
    site_col,
    sample_col,
    replicate_col,
    "OTU",
    occupancy_covars,
    all.vars(capture_formula),
    all.vars(abundance_formula)
  ))

  factor_cols <- setdiff(factor_cols, c("y", "total_reads", "z_sim", "a_sim"))

  long_df <- to_factor_cols(long_df, factor_cols)

  otu_abundances <- long_df |>
    dplyr::group_by(.data$OTU) |>
    dplyr::summarise(
      total_count = sum(.data$y, na.rm = TRUE),
      .groups = "drop"
    )

  top_otu <- otu_abundances |>
    dplyr::arrange(dplyr::desc(.data$total_count)) |>
    dplyr::slice(1) |>
    dplyr::pull(.data$OTU)

  long_df$OTU <- stats::relevel(factor(long_df$OTU), ref = as.character(top_otu))

  # ------------------------------------------------------------
  # Define keys
  # ------------------------------------------------------------

  site_keys <- c(site_col, "OTU")
  sample_keys <- c(site_col, sample_col, "OTU")

  pcr_keys <- c(site_col, sample_col, "OTU")

  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    pcr_keys <- c(site_col, sample_col, replicate_col, "OTU")
  }

  # ------------------------------------------------------------
  # Build site-level Z
  # ------------------------------------------------------------

  reduced_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
    dplyr::summarise(
      z_obs = as.integer(any(.data$y > abundance_threshold)),
      dplyr::across(
        .cols = -dplyr::all_of("y"),
        .fns = dplyr::first
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(z_sim = .data$z_obs)

  reduced_data <- to_factor_cols(reduced_data, factor_cols)

  # ------------------------------------------------------------
  # Build biological-sample-level A
  # ------------------------------------------------------------

  cap_vars <- get_formula_vars(capture_formula, "a_sim")
  sample_keep_vars <- setdiff(cap_vars, c(sample_keys, "a_sim"))
  sample_keep_vars <- intersect(sample_keep_vars, names(long_df))

  sample_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
    dplyr::summarise(
      a_obs = as.integer(any(.data$y > abundance_threshold)),
      dplyr::across(
        dplyr::all_of(sample_keep_vars),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(a_sim = .data$a_obs)

  sample_data <- to_factor_cols(sample_data, factor_cols)

  # ------------------------------------------------------------
  # Validate capture and abundance formulas
  # ------------------------------------------------------------

  validate_formula <- function(formula, data, expected_response, model_name) {
    response <- all.vars(formula)[1]

    if (response != expected_response) {
      stop(
        model_name, " formula response must be '", expected_response,
        "', but got '", response, "'."
      )
    }

    missing_vars <- setdiff(all.vars(formula), names(data))

    if (length(missing_vars) > 0) {
      stop(
        "Missing variables in ", model_name, " formula: ",
        paste(missing_vars, collapse = ", ")
      )
    }

    invisible(TRUE)
  }

  validate_formula(capture_formula, sample_data, "a_sim", "capture")
  validate_formula(abundance_formula, long_df, "y", "abundance")

  # ------------------------------------------------------------
  # Family setup
  # ------------------------------------------------------------

  abundance_glmm_family <- switch(
    abundance_family,
    poisson = stats::poisson(),
    nbinom  = glmmTMB::nbinom2(),
    zip     = stats::poisson(),
    zinb    = glmmTMB::nbinom2()
  )

  zi_formula <- if (abundance_family %in% c("zip", "zinb")) ~1 else ~0

  # ------------------------------------------------------------
  # Helper: zero probability from abundance model
  # ------------------------------------------------------------

  get_p0_abundance <- function(fit, newdata, family_name) {

    mu <- as.numeric(stats::predict(
      fit,
      type = "response",
      newdata = newdata,
      allow.new.levels = TRUE
    ))

    zi <- rep(0, length(mu))

    if (family_name %in% c("zip", "zinb")) {
      zi <- tryCatch(
        as.numeric(stats::predict(
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

      zip = stats::dpois(0, lambda = mu),

      nbinom = {
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

  lv_sites_list <- list()
  lv_species_list <- list()

  # ------------------------------------------------------------
  # Main EM-like iterations
  # ------------------------------------------------------------

  for (i in seq_len(n_iter)) {

    if (verbose) message("Iteration ", i)

    reduced_data[[site_col]] <- as.character(reduced_data[[site_col]])
    long_df[[site_col]] <- as.character(long_df[[site_col]])

    # ----------------------------------------------------------
    # 1. GLLVM occupancy model: Z
    # ----------------------------------------------------------

    z_matrix <- reshape2::acast(
      reduced_data,
      stats::as.formula(paste(site_col, "~ OTU")),
      value.var = "z_sim",
      fill = 0
    )

    z_sites <- rownames(z_matrix)

    cov_cols <- unique(c(site_col, occupancy_covars))

    cov_df <- reduced_data |>
      dplyr::select(dplyr::all_of(cov_cols)) |>
      dplyr::distinct()

    cov_df <- to_factor_cols(cov_df, occupancy_covars)
    cov_df <- cov_df[match(z_sites, cov_df[[site_col]]), , drop = FALSE]

    X_cov <- if (!is.null(occupancy_covars) && length(occupancy_covars) > 0) {
      stats::model.matrix(
        ~ .,
        data = cov_df[, occupancy_covars, drop = FALSE]
      )[, -1, drop = FALSE]
    } else {
      NULL
    }

    model_occupancy <- gllvm::gllvm(
      y = z_matrix,
      X = X_cov,
      family = "binomial",
      num.lv = num_lv_c
    )

    occupancy_models[[i]] <- model_occupancy

    # ----------------------------------------------------------
    # Latent variables
    # ----------------------------------------------------------

    lv_sites <- as.data.frame(model_occupancy$lvs)
    colnames(lv_sites) <- paste0("LV", seq_len(ncol(lv_sites)))
    lv_sites$Iteration <- i
    lv_sites$Site <- rownames(model_occupancy$lvs)

    lv_species <- as.data.frame(model_occupancy$params$theta)
    colnames(lv_species) <- paste0("LV", seq_len(ncol(lv_species)))
    lv_species$Iteration <- i
    lv_species$OTU <- rownames(model_occupancy$params$theta)

    lv_sites_list[[i]] <- lv_sites
    lv_species_list[[i]] <- lv_species

    # ----------------------------------------------------------
    # Predict occupancy probability ψ
    # ----------------------------------------------------------

    psi_prob <- stats::predict(model_occupancy, type = "response")
    rownames(psi_prob) <- rownames(z_matrix)

    psi_long <- reshape2::melt(
      psi_prob,
      varnames = c(site_col, "OTU"),
      value.name = "psi_prob"
    )

    psi_long$psi_prob <- pmin(pmax(psi_long$psi_prob, 1e-6), 1 - 1e-6)
    psi_long$eta <- stats::qlogis(psi_long$psi_prob)

    psi_list[[i]] <- psi_long[, c(site_col, "OTU", "eta")]

    # ----------------------------------------------------------
    # 2. Capture model: A | Z = 1
    # ----------------------------------------------------------

    sample_data <- sample_data |>
      dplyr::select(-dplyr::any_of(c("z_sim", "p0_sample", "capture_prob"))) |>
      dplyr::left_join(
        reduced_data |>
          dplyr::select(dplyr::all_of(c(site_col, "OTU", "z_sim"))),
        by = c(site_col, "OTU")
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

    capture_models[[i]] <- cap_fit

    eta_capture <- as.numeric(stats::predict(
      cap_fit,
      type = "link",
      newdata = sample_data,
      allow.new.levels = TRUE
    ))

    capture_prob <- stats::plogis(eta_capture)

    capture_list[[i]] <- data.frame(
      sample_data[sample_keys],
      eta = eta_capture
    )

    # ----------------------------------------------------------
    # 3. Abundance model: Y | A = 1
    # ----------------------------------------------------------

    abundance_data <- long_df |>
      dplyr::left_join(
        sample_data |>
          dplyr::select(dplyr::all_of(c(sample_keys, "a_sim"))),
        by = sample_keys
      ) |>
      dplyr::filter(.data$a_sim == 1)

    if (nrow(abundance_data) == 0) {
      stop("No rows available for abundance model at iteration ", i)
    }

    abundance_data <- to_factor_cols(abundance_data, factor_cols)

    model_abundance <- glmmTMB::glmmTMB(
      formula = abundance_formula,
      data = abundance_data,
      family = abundance_glmm_family,
      ziformula = zi_formula
    )

    abundance_models[[i]] <- model_abundance

    lambda_pred <- stats::predict(
      model_abundance,
      type = "link",
      newdata = long_df,
      allow.new.levels = TRUE
    )

    lambda_eta <- as.numeric(lambda_pred)

    lambda_list[[i]] <- data.frame(
      long_df[site_keys],
      eta = lambda_eta
    ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(
        eta = mean(.data$eta, na.rm = TRUE),
        .groups = "drop"
      )

    # ----------------------------------------------------------
    # 4. p0 at PCR/read level for all rows
    # ----------------------------------------------------------

    pred_data_all <- long_df |>
      dplyr::left_join(
        sample_data |>
          dplyr::select(dplyr::all_of(c(sample_keys, "a_sim"))),
        by = sample_keys
      )

    pred_data_all <- to_factor_cols(pred_data_all, factor_cols)

    p0_pcr <- get_p0_abundance(
      fit = model_abundance,
      newdata = pred_data_all,
      family_name = abundance_family
    )

    # ----------------------------------------------------------
    # 5. Collapse PCR/read zeros to sample level
    # ----------------------------------------------------------

    sample_p0 <- pred_data_all |>
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
    # 6. Update latent A
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
    # 7. Collapse sample probabilities to site level
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

    # Site-level detection probability
    p_detect_list[[i]] <- site_p0 |>
      dplyr::mutate(
        p_detect = 1 - .data$p0_site
      ) |>
      dplyr::select(dplyr::all_of(site_keys), p_detect)

    # ----------------------------------------------------------
    # 8. Update latent Z
    # ----------------------------------------------------------

    reduced_data <- reduced_data |>
      dplyr::select(-dplyr::any_of("p0_site")) |>
      dplyr::left_join(site_p0, by = site_keys) |>
      dplyr::mutate(
        p0_site = tidyr::replace_na(.data$p0_site, 1)
      )

    z_merge <- reduced_data |>
      dplyr::left_join(
        psi_long[, c(site_col, "OTU", "psi_prob")],
        by = site_keys
      )

    zero_indices <- which(z_merge$z_obs == 0)

    if (length(zero_indices) > 0) {

      numerator_z <- z_merge$psi_prob[zero_indices] *
        z_merge$p0_site[zero_indices]

      denominator_z <- (1 - z_merge$psi_prob[zero_indices]) +
        z_merge$psi_prob[zero_indices] *
        z_merge$p0_site[zero_indices]

      posterior_z <- numerator_z / pmax(denominator_z, 1e-12)
      posterior_z <- pmin(pmax(posterior_z, 0.001), 0.999)

      reduced_data$z_sim[zero_indices] <- stats::rbinom(
        n = length(zero_indices),
        size = 1,
        prob = posterior_z
      )
    }

    reduced_data$z_sim[reduced_data$z_obs == 1] <- 1
  }

  # ------------------------------------------------------------
  # Summaries
  # ------------------------------------------------------------

  summarise_link <- function(lst, link_type = c("logit", "log")) {

    link_type <- match.arg(link_type)

    df <- dplyr::bind_rows(lst)

    if (nrow(df) == 0) {
      return(data.frame())
    }

    inv_link <- switch(
      link_type,
      logit = stats::plogis,
      log = exp
    )

    df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(
        mean = mean(inv_link(.data$eta), na.rm = TRUE),
        median = stats::median(inv_link(.data$eta), na.rm = TRUE),
        lwr = stats::quantile(inv_link(.data$eta), 0.025, na.rm = TRUE),
        upr = stats::quantile(inv_link(.data$eta), 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }

  summarise_capture <- function(lst) {

    df <- dplyr::bind_rows(lst)

    if (nrow(df) == 0) {
      return(data.frame())
    }

    df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
      dplyr::summarise(
        capture_mean = mean(stats::plogis(.data$eta), na.rm = TRUE),
        capture_median = stats::median(stats::plogis(.data$eta), na.rm = TRUE),
        capture_lwr = stats::quantile(stats::plogis(.data$eta), 0.025, na.rm = TRUE),
        capture_upr = stats::quantile(stats::plogis(.data$eta), 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }

  summarise_pdetect <- function(lst) {

    df <- dplyr::bind_rows(lst)

    if (nrow(df) == 0) {
      return(data.frame())
    }

    df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(
        p_detect_mean = mean(.data$p_detect, na.rm = TRUE),
        p_detect_median = stats::median(.data$p_detect, na.rm = TRUE),
        p_detect_lwr = stats::quantile(.data$p_detect, 0.025, na.rm = TRUE),
        p_detect_upr = stats::quantile(.data$p_detect, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }

  keep <- seq.int(burn_in + 1, n_iter)

  psi_summary <- summarise_link(psi_list[keep], "logit") |>
    dplyr::rename_with(~ paste0("psi_", .x), -dplyr::all_of(site_keys))

  lambda_summary <- summarise_link(lambda_list[keep], "log") |>
    dplyr::rename_with(~ paste0("lambda_", .x), -dplyr::all_of(site_keys))

  p_detect_summary <- summarise_pdetect(p_detect_list[keep])

  capture_summary <- summarise_capture(capture_list[keep])

  final_summary <- psi_summary |>
    dplyr::left_join(lambda_summary, by = site_keys) |>
    dplyr::left_join(p_detect_summary, by = site_keys)

  # ------------------------------------------------------------
  # Latent variables
  # ------------------------------------------------------------

  lv_sites_combined <- dplyr::bind_rows(lv_sites_list[keep])
  lv_species_combined <- dplyr::bind_rows(lv_species_list[keep])

  mean_lv_sites <- lv_sites_combined |>
    dplyr::group_by(.data$Site) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::starts_with("LV"),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  mean_lv_species <- lv_species_combined |>
    dplyr::group_by(.data$OTU) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::starts_with("LV"),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  list(
    summary = final_summary,
    capture = capture_summary,

    psi_list = psi_list[keep],
    capture_list = capture_list[keep],
    lambda_list = lambda_list[keep],
    p_detect_list = p_detect_list[keep],

    occupancy_models = occupancy_models[keep],
    capture_models = capture_models[keep],
    abundance_models = abundance_models[keep],

    reduced_data = reduced_data,
    sample_data = sample_data,
    long_df = long_df,

    lv_sites = lv_sites_combined,
    lv_species = lv_species_combined,
    mean_lv_sites = mean_lv_sites,
    mean_lv_species = mean_lv_species,

    filter_summary = list(
      otu_stats = otu_stats,
      kept_otus = keep_otus,
      min_species_sum = min_species_sum,
      min_detection_replicates = min_detection_replicates,
      abundance_threshold = abundance_threshold
    ),

    note = paste(
      "Three-level GLLVM-GLMM model.",
      "Hierarchy: site occupancy Z estimated by gllvm, biological-sample capture A estimated by GLMM,",
      "and replicate-level abundance Y estimated by GLMM."
    )
  )
}
                          
