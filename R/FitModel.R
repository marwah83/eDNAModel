#' Fit a Hierarchical Occupancy–Detection–Abundance Model for Microbiome Data
#'
#' Implements an EM-like iterative framework using repeated GLMM fitting
#' to jointly estimate microbial occupancy, detection (capture), and abundance
#' for OTU-level microbiome data stored in a `phyloseq` object.
#'
#' The model explicitly follows a **three-level eDNA hierarchy**:
#' \itemize{
#'   \item Site level: true presence/absence (\eqn{z})
#'   \item Biological sample level: capture/detection (\eqn{a})
#'   \item PCR replicate level: observed counts (\eqn{y})
#' }
#'
#' Formula covariates are automatically preserved and propagated into the
#' internally constructed site-level and sample-level datasets. This allows
#' covariates such as `in_out`, `Samplingmonth`, treatments, interactions,
#' and nested random effects to be used directly in model formulas.
#'
#' @param phyloseq A `phyloseq` object containing OTU table and sample data.
#' @param site_col Character. Column representing sites (occupancy level).
#' @param sample_col Character. Column representing biological samples (capture level).
#' @param replicate_col Character or `NULL`. Column representing PCR replicates.
#' @param occupancy_formula Formula for occupancy model (response must be `z_sim`).
#' @param capture_formula Formula for capture model (response must be `a_sim`).
#' @param abundance_formula Formula for abundance model (response must match `count_col`).
#' @param otu_col Character. OTU identifier column.
#' @param count_col Character. Count column (e.g., read counts).
#' @param min_species_sum Minimum total counts required to retain an OTU.
#' @param min_detection_replicates Minimum detections required per OTU.
#' @param abundance_threshold Threshold defining detection (default = 0).
#' @param n_iter Number of EM-like iterations.
#' @param burn_in Number of initial iterations discarded.
#' @param abundance_family One of `"poisson"`, `"nbinom"`, `"zip"`, `"zinb"`.
#' @param verbose Logical; print progress.
#'
#' @return A list with:
#' \describe{
#'   \item{psi}{Occupancy probabilities (`psi_mean`, `psi_median`, `psi_lwr`, `psi_upr`).}
#'   \item{capture}{Capture probabilities (`capture_mean`, etc.).}
#'   \item{lambda}{Abundance (`lambda_mean`, etc.).}
#'   \item{p_detect}{Detection probability derived from abundance.}
#'   \item{psi_list, capture_list, lambda_list, p_detect_list}{Per-iteration draws.}
#'   \item{occupancy_models, capture_models, abundance_models}{Fitted GLMMs.}
#'   \item{site_data, sample_data, long_df}{Processed data.}
#'   \item{filter_summary}{OTU filtering summary.}
#'   \item{diagnostic_AIC}{AIC values per iteration.}
#'   \item{note}{Model description.}
#' }
#'
#' @details
#' The model decomposes microbial occurrence into:
#'
#' \enumerate{
#'   \item \strong{Occupancy (\eqn{\psi})}:
#'   Probability an OTU is present at a site.
#'
#'   \item \strong{Capture (\eqn{p})}:
#'   Probability a biological sample captures DNA given presence.
#'
#'   \item \strong{Abundance (\eqn{\lambda})}:
#'   Expected read counts at PCR replicate level given capture.
#' }
#'
#' The hierarchical structure is:
#'
#' \deqn{
#' z_{i} \sim \text{Bernoulli}(\psi_i), \quad
#' a_{ij} \sim \text{Bernoulli}(p_{ij} \cdot z_i), \quad
#' y_{ijk} \sim \text{Count}(\lambda_{ijk} \cdot a_{ij})
#' }
#'
#' Detection probability implied by abundance is:
#'
#' \deqn{p_{\text{detect}} = 1 - \exp(-\lambda)}
#'
#'
#' \strong{Important interpretation:}
#' \itemize{
#'   \item For small \eqn{\lambda}: \eqn{p_{\text{detect}} \approx \lambda}
#'   \item For large \eqn{\lambda}: \eqn{p_{\text{detect}} \to 1}
#' }
#'
#' Thus, `p_detect` is not independent of `lambda` and should not be interpreted
#' as a separate parameter.
#'
#' @section EM-like algorithm:
#' \enumerate{
#'   \item Fit occupancy model using current \eqn{z_sim}
#'   \item Fit capture model conditional on \eqn{z_sim = 1}
#'   \item Fit abundance model conditional on \eqn{a_sim = 1}
#'   \item Update \eqn{a_sim} using capture + abundance
#'   \item Update \eqn{z_sim} using capture + abundance
#' }
#'
#' Posterior summaries are computed from empirical distributions across
#' iterations (after burn-in), i.e., **true empirical quantiles**, not Wald intervals.
#'
#' @section Model features:
#' \itemize{
#'   \item Explicit 3-level eDNA hierarchy (site → sample → PCR replicate)
#'   \item Supports complex random-effects structures
#'   \item Automatically handles covariates from `sample_data`
#'   \item Supports offsets via `offset(log(total_reads))`
#'   \item Works with Poisson, NB, ZIP, ZINB
#' }
#'
#' @section Caveats:
#' \itemize{
#'   \item This is an approximate EM-like method, not a full joint likelihood.
#'   \item AIC values are component-wise diagnostics only.
#'   \item Detection from abundance assumes Poisson-type logic (approximate for NB/ZIP).
#'   \item Large \eqn{\lambda} leads to saturated detection probabilities (~1).
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
    verbose = TRUE
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
  # Compute total_reads from processed data
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
  # Build site-level Z and biological-sample-level A
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
    # 1. Occupancy model: Z_site,OTU
    # ----------------------------------------------------------
    
    occ_fit <- glmmTMB::glmmTMB(
      formula = occupancy_formula,
      data = site_data,
      family = stats::binomial()
    )
    
    eta_psi <- as.numeric(predict(
      occ_fit,
      type = "link",
      newdata = site_data,
      allow.new.levels = TRUE
    ))
    
    psi <- stats::plogis(eta_psi)
    
    psi_list[[i]] <- data.frame(
      site_data[site_keys],
      eta = eta_psi
    )
    
    occupancy_models[[i]] <- occ_fit
    
    # ----------------------------------------------------------
    # 2. Capture model: A_sample,OTU | Z = 1
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
    
    eta_capture <- as.numeric(predict(
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
    
    capture_models[[i]] <- cap_fit
    
    # ----------------------------------------------------------
    # 3. Abundance model: PCR reads Y | A = 1
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
    
    # Predict conditional abundance for all PCR rows
    eta_lambda_all <- as.numeric(predict(
      abund_fit,
      type = "link",
      newdata = long_df,
      allow.new.levels = TRUE
    ))
    
    lambda_all <- exp(eta_lambda_all)
    
    lambda_list[[i]] <- data.frame(
      long_df[pcr_keys],
      eta = eta_lambda_all
    )
    
    # PCR-level P(Y = 0 | A = 1)
    p0_pcr <- get_p0_pcr(
      fit = abund_fit,
      newdata = long_df,
      family_name = abundance_family
    )
    
    p_detect_pcr <- 1 - p0_pcr
    p_detect_eta <- log(-log(pmax(p0_pcr, 1e-12)))
    
    p_detect_list[[i]] <- data.frame(
      long_df[pcr_keys],
      eta = p_detect_eta
    )
    
    abundance_models[[i]] <- abund_fit
    
    # ----------------------------------------------------------
    # 4. Collapse PCR probabilities to biological sample level
    #    P(all PCR replicates zero | A = 1)
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
    # 5. Update A_sample,OTU
    #    If any PCR is positive, A = 1.
    #    If all PCRs are zero, sample A using posterior probability.
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
    # 6. Update Z_site,OTU
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
  # Summaries on natural scale
  # ------------------------------------------------------------
  
  summarise_link <- function(lst, link_name, prefix) {
    
    inv_link <- switch(
      link_name,
      logit   = stats::plogis,
      log     = exp,
      cloglog = function(x) 1 - exp(-exp(x)),
      stop("Unknown link_name: ", link_name)
    )
    
    df <- dplyr::bind_rows(lst)
    
    if (nrow(df) == 0) return(data.frame())
    
    keys <- setdiff(names(df), "eta")
    
    out <- df |>
      dplyr::mutate(value = inv_link(.data$eta)) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::summarise(
        mean = mean(.data$value, na.rm = TRUE),
        median = stats::median(.data$value, na.rm = TRUE),
        lwr = stats::quantile(.data$value, 0.025, na.rm = TRUE),
        upr = stats::quantile(.data$value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    
    names(out)[names(out) == "mean"]   <- paste0(prefix, "_mean")
    names(out)[names(out) == "median"] <- paste0(prefix, "_median")
    names(out)[names(out) == "lwr"]    <- paste0(prefix, "_lwr")
    names(out)[names(out) == "upr"]    <- paste0(prefix, "_upr")
    
    out
  }
  
  keep <- seq.int(burn_in + 1, n_iter)
  
  list(
    psi = summarise_link(psi_list[keep], "logit", "psi"),
    capture = summarise_link(capture_list[keep], "logit", "capture"),
    lambda = summarise_link(lambda_list[keep], "log", "lambda"),
    p_detect = summarise_link(p_detect_list[keep], "cloglog", "p_detect"),
    
    psi_list = psi_list[keep],
    capture_list = capture_list[keep],
    lambda_list = lambda_list[keep],
    p_detect_list = p_detect_list[keep],
    
    occupancy_models = occupancy_models,
    capture_models = capture_models,
    abundance_models = abundance_models,
    
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
      "Formula covariates are automatically preserved in site/sample data.",
      "AIC values are exploratory component-model diagnostics, not a formal joint-likelihood criterion."
    )
  )
}
                          
  
