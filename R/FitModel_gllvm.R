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
  
  abundance_family <- match.arg(abundance_family)
  
  if (burn_in >= n_iter) {
    stop("burn_in must be smaller than n_iter.")
  }
  
  bt <- function(x) paste0("`", x, "`")
  
  to_factor_cols <- function(df, cols) {
    cols <- unique(stats::na.omit(cols))
    for (col in cols) {
      if (col %in% names(df)) {
        df[[col]] <- as.factor(as.character(df[[col]]))
      }
    }
    df
  }
  
  get_formula_vars <- function(formula, response) {
    setdiff(all.vars(formula), response)
  }
  
  # ------------------------------------------------------------
  # Handle abundance formula (RHS OR full formula)
  # ------------------------------------------------------------
  
  if (inherits(abundance_rhs, "formula")) {
    
    # user provided full formula: y ~ ...
    abundance_formula <- abundance_rhs
    
  } else {
    
    # user provided RHS: 1 + (1 | OTU)
    abundance_rhs_expr <- substitute(abundance_rhs)
    abundance_rhs_txt <- paste(deparse(abundance_rhs_expr), collapse = " ")
    
    abundance_formula <- stats::as.formula(
      paste(bt(count_col), "~", abundance_rhs_txt),
      env = parent.frame()
    )
  }
  
  if (is.null(capture_formula)) {
    capture_formula <- stats::as.formula(
      paste("a_sim ~ 1 + (1 |", bt(otu_col), ")"),
      env = parent.frame()
    )
  }
  
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
  
  long_df[[site_col]]   <- as.character(long_df[[site_col]])
  long_df[[sample_col]] <- as.character(long_df[[sample_col]])
  long_df[[otu_col]]    <- as.character(long_df[[otu_col]])
  long_df[[count_col]]  <- as.numeric(long_df[[count_col]])
  
  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    long_df[[replicate_col]] <- as.character(long_df[[replicate_col]])
  }
  
  long_df <- long_df[!is.na(long_df[[count_col]]), ]
  
  total_reads_df <- long_df |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      total_reads = sum(.data[[count_col]], na.rm = TRUE),
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
  
  if (verbose) {
    message("OTUs before filtering: ", dplyr::n_distinct(otu_stats[[otu_col]]))
    message("OTUs after filtering: ", dplyr::n_distinct(long_df[[otu_col]]))
  }
  
  factor_cols <- unique(c(
    site_col,
    sample_col,
    replicate_col,
    otu_col,
    occupancy_covars,
    all.vars(capture_formula),
    all.vars(abundance_formula)
  ))
  
  factor_cols <- setdiff(
    factor_cols,
    c(count_col, "total_reads", "z_sim", "a_sim")
  )
  
  long_df <- to_factor_cols(long_df, factor_cols)
  
  otu_abundances <- long_df |>
    dplyr::group_by(.data[[otu_col]]) |>
    dplyr::summarise(
      total_count = sum(.data[[count_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  top_otu <- otu_abundances |>
    dplyr::arrange(dplyr::desc(.data$total_count)) |>
    dplyr::slice(1) |>
    dplyr::pull(.data[[otu_col]])
  
  long_df[[otu_col]] <- stats::relevel(
    factor(long_df[[otu_col]]),
    ref = as.character(top_otu)
  )
  
  site_keys <- c(site_col, otu_col)
  sample_keys <- c(site_col, sample_col, otu_col)
  
  site_keep_vars <- intersect(occupancy_covars, names(long_df))
  
  reduced_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
    dplyr::summarise(
      z_obs = as.integer(any(.data[[count_col]] > abundance_threshold)),
      dplyr::across(dplyr::all_of(site_keep_vars), ~ dplyr::first(.x)),
      .groups = "drop"
    ) |>
    dplyr::mutate(z_sim = .data$z_obs)
  
  reduced_data <- to_factor_cols(reduced_data, factor_cols)
  
  cap_vars <- get_formula_vars(capture_formula, "a_sim")
  sample_keep_vars <- setdiff(cap_vars, c(sample_keys, "a_sim"))
  sample_keep_vars <- intersect(sample_keep_vars, names(long_df))
  
  sample_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
    dplyr::summarise(
      a_obs = as.integer(any(.data[[count_col]] > abundance_threshold)),
      dplyr::across(dplyr::all_of(sample_keep_vars), ~ dplyr::first(.x)),
      .groups = "drop"
    ) |>
    dplyr::mutate(a_sim = .data$a_obs)
  
  sample_data <- to_factor_cols(sample_data, factor_cols)
  
  validate_formula <- function(formula, data, expected_response, model_name) {
    response <- all.vars(formula)[1]
    
    if (response != expected_response) {
      stop(model_name, " formula response must be '", expected_response,
           "', but got '", response, "'.")
    }
    
    missing_vars <- setdiff(all.vars(formula), names(data))
    
    if (length(missing_vars) > 0) {
      stop("Missing variables in ", model_name, " formula: ",
           paste(missing_vars, collapse = ", "))
    }
    
    invisible(TRUE)
  }
  
  validate_formula(capture_formula, sample_data, "a_sim", "capture")
  validate_formula(abundance_formula, long_df, count_col, "abundance")
  
  abundance_glmm_family <- switch(
    abundance_family,
    poisson = stats::poisson(),
    nbinom  = glmmTMB::nbinom2(),
    zip     = stats::poisson(),
    zinb    = glmmTMB::nbinom2()
  )
  
  zi_formula <- if (abundance_family %in% c("zip", "zinb")) ~1 else ~0
  
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
  
  ## ------------------------------------------------------------
  ## Precompute site x OTU template and occupancy covariates once
  ## ------------------------------------------------------------
  
  z_template <- reshape2::acast(
    reduced_data,
    stats::as.formula(paste(bt(site_col), "~", bt(otu_col))),
    value.var = "z_sim",
    fill = 0
  )
  
  z_sites <- rownames(z_template)
  z_otus  <- colnames(z_template)
  
  cov_cols <- unique(c(site_col, occupancy_covars))
  
  cov_df <- reduced_data |>
    dplyr::select(dplyr::all_of(cov_cols)) |>
    dplyr::distinct()
  
  cov_df[[site_col]] <- as.character(cov_df[[site_col]])
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
  
  psi_list <- vector("list", n_iter)
  capture_list <- vector("list", n_iter)
  lambda_list <- vector("list", n_iter)
  p_detect_list <- vector("list", n_iter)
  
  occupancy_models <- vector("list", n_iter)
  capture_models <- vector("list", n_iter)
  abundance_models <- vector("list", n_iter)
  
  lv_sites_list <- vector("list", n_iter)
  lv_species_list <- vector("list", n_iter)
  
  diagnostic_AIC <- data.frame(
    iteration = integer(),
    occupancy_AIC = numeric(),
    capture_AIC = numeric(),
    abundance_AIC = numeric(),
    abundance_logLik = numeric()
  )
  
  for (i in seq_len(n_iter)) {
    
    if (verbose) message("Iteration ", i)
    t_iter <- Sys.time()
    
    ## Update z_matrix without rebuilding acast
    z_matrix <- z_template
    z_matrix[,] <- 0
    
    idx <- cbind(
      match(as.character(reduced_data[[site_col]]), z_sites),
      match(as.character(reduced_data[[otu_col]]), z_otus)
    )
    
    z_matrix[idx] <- reduced_data$z_sim
    
    model_occupancy <- tryCatch(
      gllvm::gllvm(
        y = z_matrix,
        X = X_cov,
        family = "binomial",
        num.lv = num_lv_c
      ),
      error = function(e) {
        message("GLLVM failed at iteration ", i, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(model_occupancy)) next
    
    occupancy_models[[i]] <- model_occupancy
    
    ## Safe latent variables: works also when num_lv_c = 0
    if (!is.null(model_occupancy$lvs) &&
        length(model_occupancy$lvs) > 0 &&
        ncol(as.data.frame(model_occupancy$lvs)) > 0) {
      
      lv_sites <- as.data.frame(model_occupancy$lvs)
      colnames(lv_sites) <- paste0("LV", seq_len(ncol(lv_sites)))
      lv_sites[[site_col]] <- rownames(model_occupancy$lvs)
      lv_sites$Iteration <- i
      
    } else {
      lv_sites <- data.frame()
    }
    
    if (!is.null(model_occupancy$params$theta) &&
        length(model_occupancy$params$theta) > 0 &&
        ncol(as.matrix(model_occupancy$params$theta)) > 0) {
      
      theta_mat <- as.matrix(model_occupancy$params$theta)
      
      lv_species <- as.data.frame(theta_mat)
      colnames(lv_species) <- paste0("LV", seq_len(ncol(lv_species)))
      lv_species[[otu_col]] <- rownames(theta_mat)
      lv_species$Iteration <- i
      
    } else {
      lv_species <- data.frame()
    }
    
    lv_sites_list[[i]] <- lv_sites
    lv_species_list[[i]] <- lv_species
    
    psi_prob <- stats::predict(model_occupancy, type = "response")
    rownames(psi_prob) <- rownames(z_matrix)
    colnames(psi_prob) <- colnames(z_matrix)
    
    psi_long <- reshape2::melt(
      psi_prob,
      varnames = c(site_col, otu_col),
      value.name = "psi_prob"
    )
    
    psi_long$psi_prob <- pmin(pmax(psi_long$psi_prob, 1e-6), 1 - 1e-6)
    psi_long$eta <- stats::qlogis(psi_long$psi_prob)
    
    psi_list[[i]] <- psi_long[, c(site_col, otu_col, "eta")]
    
    sample_data <- sample_data |>
      dplyr::select(-dplyr::any_of(c("z_sim", "p0_sample", "capture_prob"))) |>
      dplyr::left_join(
        reduced_data |>
          dplyr::select(dplyr::all_of(c(site_keys, "z_sim"))),
        by = site_keys
      )
    
    capture_fit_data <- sample_data |>
      dplyr::filter(.data$z_sim == 1)
    
    if (nrow(capture_fit_data) == 0) {
      stop("No rows available for capture model at iteration ", i)
    }
    
    cap_fit <- tryCatch(
      glmmTMB::glmmTMB(
        formula = capture_formula,
        data = capture_fit_data,
        family = stats::binomial()
      ),
      error = function(e) {
        message("Capture model failed at iteration ", i, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(cap_fit)) next
    
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
    
    model_abundance <- tryCatch(
      glmmTMB::glmmTMB(
        formula = abundance_formula,
        data = abundance_data,
        family = abundance_glmm_family,
        ziformula = zi_formula
      ),
      error = function(e) {
        message("Abundance model failed at iteration ", i, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(model_abundance)) next
    
    abundance_models[[i]] <- model_abundance
    
    lambda_eta <- as.numeric(stats::predict(
      model_abundance,
      type = "link",
      newdata = long_df,
      allow.new.levels = TRUE
    ))
    
    lambda_list[[i]] <- data.frame(
      long_df[site_keys],
      eta = lambda_eta
    ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(
        eta = mean(.data$eta, na.rm = TRUE),
        .groups = "drop"
      )
    
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
    
    zero_sample <- which(sample_data$a_obs == 0 & sample_data$z_sim == 1)
    
    if (length(zero_sample) > 0) {
      numerator_a <- sample_data$capture_prob[zero_sample] *
        sample_data$p0_sample[zero_sample]
      
      denominator_a <- (1 - sample_data$capture_prob[zero_sample]) +
        sample_data$capture_prob[zero_sample] *
        sample_data$p0_sample[zero_sample]
      
      posterior_a <- numerator_a / pmax(denominator_a, 1e-12)
      posterior_a <- pmin(pmax(posterior_a, 1e-4), 1 - 1e-4)
      
      sample_data$a_sim[zero_sample] <- stats::rbinom(
        n = length(zero_sample),
        size = 1,
        prob = posterior_a
      )
    }
    
    sample_data$a_sim[sample_data$a_obs == 1] <- 1
    sample_data$a_sim[sample_data$z_sim == 0] <- 0
    
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
    
    p_detect_list[[i]] <- site_p0 |>
      dplyr::mutate(
        eta = log(-log(pmax(.data$p0_site, 1e-12)))
      ) |>
      dplyr::select(dplyr::all_of(site_keys), eta)
    
    reduced_data <- reduced_data |>
      dplyr::select(-dplyr::any_of("p0_site")) |>
      dplyr::left_join(site_p0, by = site_keys) |>
      dplyr::mutate(
        p0_site = tidyr::replace_na(.data$p0_site, 1)
      )
    
    z_merge <- reduced_data |>
      dplyr::left_join(
        psi_long[, c(site_col, otu_col, "psi_prob")],
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
      posterior_z <- pmin(pmax(posterior_z, 1e-4), 1 - 1e-4)
      
      reduced_data$z_sim[zero_indices] <- stats::rbinom(
        n = length(zero_indices),
        size = 1,
        prob = posterior_z
      )
    }
    
    reduced_data$z_sim[reduced_data$z_obs == 1] <- 1
    
    diagnostic_AIC <- rbind(
      diagnostic_AIC,
      data.frame(
        iteration = i,
        occupancy_AIC = tryCatch(stats::AIC(model_occupancy), error = function(e) NA_real_),
        capture_AIC = tryCatch(stats::AIC(cap_fit), error = function(e) NA_real_),
        abundance_AIC = tryCatch(stats::AIC(model_abundance), error = function(e) NA_real_),
        abundance_logLik = tryCatch(as.numeric(stats::logLik(model_abundance)), error = function(e) NA_real_)
      )
    )
    
    if (verbose) {
      message(
        "Iteration ", i, " finished in ",
        round(difftime(Sys.time(), t_iter, units = "secs"), 2),
        " seconds"
      )
    }
  }
  
  keep <- seq.int(burn_in + 1, n_iter)
  
  summarise_link <- function(lst, link_type = c("logit", "log", "cloglog")) {
    
    link_type <- match.arg(link_type)
    
    df <- dplyr::bind_rows(lst)
    
    if (nrow(df) == 0) {
      return(data.frame())
    }
    
    inv_link <- switch(
      link_type,
      logit = stats::plogis,
      log = exp,
      cloglog = function(x) 1 - exp(-exp(x))
    )
    
    df |>
      dplyr::mutate(value = inv_link(.data$eta)) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(
        mean = mean(.data$value, na.rm = TRUE),
        median = stats::median(.data$value, na.rm = TRUE),
        lwr = stats::quantile(.data$value, 0.025, na.rm = TRUE),
        upr = stats::quantile(.data$value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  summarise_capture <- function(lst) {
    
    df <- dplyr::bind_rows(lst)
    
    if (nrow(df) == 0) {
      return(data.frame())
    }
    
    df |>
      dplyr::mutate(capture = stats::plogis(.data$eta)) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
      dplyr::summarise(
        capture_mean = mean(.data$capture, na.rm = TRUE),
        capture_median = stats::median(.data$capture, na.rm = TRUE),
        capture_lwr = stats::quantile(.data$capture, 0.025, na.rm = TRUE),
        capture_upr = stats::quantile(.data$capture, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  summarise_capture_site <- function(capture_summary) {
    
    if (nrow(capture_summary) == 0) {
      return(data.frame())
    }
    
    capture_summary |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(
        capture_mean = mean(.data$capture_mean, na.rm = TRUE),
        capture_median = stats::median(.data$capture_median, na.rm = TRUE),
        capture_lwr = stats::quantile(.data$capture_lwr, 0.025, na.rm = TRUE),
        capture_upr = stats::quantile(.data$capture_upr, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  psi_summary <- summarise_link(psi_list[keep], "logit") |>
    dplyr::rename_with(~ paste0("psi_", .x), -dplyr::all_of(site_keys))
  
  lambda_summary <- summarise_link(lambda_list[keep], "log") |>
    dplyr::rename_with(~ paste0("lambda_", .x), -dplyr::all_of(site_keys))
  
  p_detect_summary <- summarise_link(p_detect_list[keep], "cloglog") |>
    dplyr::rename_with(~ paste0("p_detect_", .x), -dplyr::all_of(site_keys))
  
  capture_summary <- summarise_capture(capture_list[keep])
  capture_site_summary <- summarise_capture_site(capture_summary)
  
  final_summary <- psi_summary |>
    dplyr::left_join(capture_site_summary, by = site_keys) |>
    dplyr::left_join(lambda_summary, by = site_keys) |>
    dplyr::left_join(p_detect_summary, by = site_keys)
  
  lv_sites_combined <- dplyr::bind_rows(lv_sites_list[keep])
  lv_species_combined <- dplyr::bind_rows(lv_species_list[keep])
  
  mean_lv_sites <- if (nrow(lv_sites_combined) > 0) {
    lv_sites_combined |>
      dplyr::group_by(.data[[site_col]]) |>
      dplyr::summarise(
        dplyr::across(dplyr::starts_with("LV"), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )
  } else {
    data.frame()
  }
  
  mean_lv_species <- if (nrow(lv_species_combined) > 0) {
    lv_species_combined |>
      dplyr::group_by(.data[[otu_col]]) |>
      dplyr::summarise(
        dplyr::across(dplyr::starts_with("LV"), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )
  } else {
    data.frame()
  }
  
  list(
    summary = final_summary,
    capture = capture_summary,
    capture_site = capture_site_summary,
    
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
    
    diagnostic_AIC = diagnostic_AIC,
    
    note = paste(
      "Three-level GLLVM-GLMM model.",
      "Hierarchy: site occupancy Z estimated by gllvm, biological-sample capture A estimated by GLMM,",
      "and replicate-level abundance Y estimated by GLMM."
    )
  )
}
