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
#' @param n_sim Number of simulations used for uncertainty propagation.
#'
#' @return A list containing:
#' \describe{
#'   \item{psi}{Occupancy probability summaries (\code{psi_mean}, \code{psi_lwr}, etc.).}
#'   \item{capture}{Capture probability summaries.}
#'   \item{lambda}{Abundance (expected counts) summaries.}
#'   \item{p_detect}{Detection probability derived from abundance.}
#'   \item{psi_list, capture_list, lambda_list, p_detect_list}{Per-iteration values (eta and SE).}
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
#' For a Poisson model, this simplifies to:
#'
#' \deqn{
#' p_{\text{detect}} = 1 - \exp(-\lambda)
#' }
#'
#' For NB / ZIP / ZINB models, the probability of zero is computed using the
#' corresponding distribution (including zero-inflation when applicable).
#'
#' \strong{Important:}
#' \itemize{
#'   \item \code{p_detect} is a deterministic function of the abundance model.
#'   \item It is not an independently estimated parameter.
#' }
#'
#' @section EM-like algorithm:
#' Each iteration performs:
#' \enumerate{
#'   \item Fit occupancy model using current \eqn{Z}
#'   \item Fit capture model conditional on \eqn{Z = 1}
#'   \item Fit abundance model conditional on \eqn{A = 1}
#'   \item Update \eqn{A} using capture + abundance probabilities
#'   \item Update \eqn{Z} using capture + abundance probabilities
#' }
#'
#' @section Uncertainty estimation:
#' Uncertainty is propagated by simulation:
#'
#' \itemize{
#'   \item Linear predictors (\eqn{\eta}) are treated as:
#'   \deqn{\eta \sim \mathcal{N}(\hat{\eta}, \text{SE}^2)}
#'
#'   \item Simulated values are transformed to the natural scale
#'   \item Empirical summaries (mean, median, 95\% intervals) are computed
#' }
#'
#' This approach avoids reliance on asymptotic (Wald) intervals and correctly
#' accounts for link-function nonlinearity.
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
#'   \item This is an approximate EM-like method, not a full joint likelihood model
#'   \item AIC values are component-wise diagnostics only
#'   \item Detection from abundance depends on distributional assumptions
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
    n_sim = 200
) {
  
  abundance_family <- match.arg(abundance_family)
  
  if (is.null(otu_col)) stop("Please specify otu_col.")
  if (is.null(count_col)) stop("Please specify count_col.")
  if (burn_in >= n_iter) stop("burn_in must be < n_iter.")
  
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
  # Type handling
  # ------------------------------------------------------------
  
  long_df[[site_col]]   <- as.factor(long_df[[site_col]])
  long_df[[otu_col]]    <- as.factor(long_df[[otu_col]])
  long_df[[sample_col]] <- as.factor(long_df[[sample_col]])
  long_df[[count_col]]  <- as.numeric(long_df[[count_col]])
  
  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    long_df[[replicate_col]] <- as.factor(long_df[[replicate_col]])
  }
  
  long_df <- long_df[!is.na(long_df[[count_col]]), ]
  
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
  
  long_df[[otu_col]] <- droplevels(long_df[[otu_col]])
  
  if (verbose) {
    message("OTUs before filtering: ", dplyr::n_distinct(otu_stats[[otu_col]]))
    message("OTUs after filtering: ", dplyr::n_distinct(long_df[[otu_col]]))
  }
  
  # ------------------------------------------------------------
  # Keys
  # ------------------------------------------------------------
  
  site_keys <- c(site_col, otu_col)
  sample_keys <- c(site_col, sample_col, otu_col)
  
  pcr_keys <- c(site_col, sample_col, otu_col)
  
  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    pcr_keys <- c(site_col, sample_col, replicate_col, otu_col)
  }
  
  # ------------------------------------------------------------
  # Build site-level Z and sample-level A
  # ------------------------------------------------------------
  
  site_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
    dplyr::summarise(
      z_obs = as.integer(any(.data[[count_col]] > abundance_threshold)),
      .groups = "drop"
    ) |>
    dplyr::mutate(z_sim = .data$z_obs)
  
  sample_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
    dplyr::summarise(
      a_obs = as.integer(any(.data[[count_col]] > abundance_threshold)),
      .groups = "drop"
    ) |>
    dplyr::mutate(a_sim = .data$a_obs)
  
  # ------------------------------------------------------------
  # Family setup
  # ------------------------------------------------------------
  
  fam <- switch(
    abundance_family,
    poisson = poisson(),
    nbinom  = nbinom2(),
    zip     = poisson(),
    zinb    = nbinom2()
  )
  
  
  zi_formula <- if (abundance_family %in% c("zip", "zinb")) ~1 else ~0
  
  # ------------------------------------------------------------
  # Helper: p_detect from simulated abundance eta
  # ------------------------------------------------------------
  
  compute_pdetect <- function(eta_sim, fit, family_name, zi_prob = 0) {
    
    mu <- exp(eta_sim)
    
    p0_cond <- switch(
      family_name,
      poisson = stats::dpois(0, lambda = mu),
      zip     = stats::dpois(0, lambda = mu),
      
      nbinom = {
        theta <- tryCatch(glmmTMB::sigma(fit), error = function(e) NA_real_)
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
    
    p0 <- zi_prob + (1 - zi_prob) * p0_cond
    1 - pmin(pmax(p0, 1e-12), 1)
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
    # 1. Occupancy model: Z
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
      se.fit = TRUE,
      allow.new.levels = TRUE
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
    # 2. Capture model: A | Z = 1
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
      se.fit = TRUE,
      allow.new.levels = TRUE
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
    # 3. Abundance model: Y | A = 1
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
      se.fit = TRUE,
      allow.new.levels = TRUE
    )
    
    eta_lambda <- as.numeric(pred_lambda$fit)
    se_lambda <- as.numeric(pred_lambda$se.fit)
    
    lambda_list[[i]] <- data.frame(
      long_df[pcr_keys],
      eta = eta_lambda,
      se = se_lambda
    )
    
    # ----------------------------------------------------------
    # 4. p_detect with uncertainty propagated from eta_lambda
    # ----------------------------------------------------------
    
    zi_prob <- rep(0, length(eta_lambda))
    
    if (abundance_family %in% c("zip", "zinb")) {
      zi_prob <- tryCatch(
        as.numeric(predict(
          abund_fit,
          type = "zprob",
          newdata = long_df,
          allow.new.levels = TRUE
        )),
        error = function(e) rep(0, length(eta_lambda))
      )
    }
    
    p_detect_sim <- lapply(seq_along(eta_lambda), function(j) {
      
      mu_j <- eta_lambda[j]
      se_j <- se_lambda[j]
      zi_j <- zi_prob[j]
      
      if (is.na(se_j) || se_j <= 0) {
        eta_sim <- mu_j
      } else {
        eta_sim <- stats::rnorm(n_sim, mean = mu_j, sd = se_j)
      }
      
      compute_pdetect(
        eta_sim = eta_sim,
        fit = abund_fit,
        family_name = abundance_family,
        zi_prob = zi_j
      )
    })
    
    p_detect_mean <- vapply(p_detect_sim, mean, numeric(1), na.rm = TRUE)
    p_detect_sd   <- vapply(p_detect_sim, stats::sd, numeric(1), na.rm = TRUE)
    
    p_detect_list[[i]] <- data.frame(
      long_df[pcr_keys],
      eta = p_detect_mean,
      se = p_detect_sd
    )
    
    abundance_models[[i]] <- abund_fit
    
    # ----------------------------------------------------------
    # 5. Collapse replicate probabilities to biological sample
    # ----------------------------------------------------------
    
    p0_pcr <- 1 - p_detect_mean
    
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
    # 6. Update A
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
    # 7. Collapse to site level and update Z
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
  # Summaries
  # ------------------------------------------------------------
  
  summarise_link <- function(lst, link_name, prefix) {
    
    inv_link <- switch(
      link_name,
      logit    = stats::plogis,
      log      = exp,
      identity = function(x) x,
      stop("Unknown link_name: ", link_name)
    )
    
    df <- dplyr::bind_rows(lst)
    
    if (nrow(df) == 0) return(data.frame())
    
    keys <- setdiff(names(df), c("eta", "se"))
    
    out <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::summarise(
        
        mean = {
          sim_vals <- unlist(
            lapply(seq_along(eta), function(j) {
              if (is.na(se[j]) || se[j] <= 0) {
                inv_link(eta[j])
              } else {
                inv_link(rnorm(n_sim, eta[j], se[j]))
              }
            })
          )
          mean(sim_vals, na.rm = TRUE)
        },
        
        median = {
          sim_vals <- unlist(
            lapply(seq_along(eta), function(j) {
              if (is.na(se[j]) || se[j] <= 0) {
                inv_link(eta[j])
              } else {
                inv_link(rnorm(n_sim, eta[j], se[j]))
              }
            })
          )
          stats::median(sim_vals, na.rm = TRUE)
        },
        
        lwr = {
          sim_vals <- unlist(
            lapply(seq_along(eta), function(j) {
              if (is.na(se[j]) || se[j] <= 0) {
                inv_link(eta[j])
              } else {
                inv_link(rnorm(n_sim, eta[j], se[j]))
              }
            })
          )
          as.numeric(stats::quantile(sim_vals, 0.025, na.rm = TRUE))
        },
        
        upr = {
          sim_vals <- unlist(
            lapply(seq_along(eta), function(j) {
              if (is.na(se[j]) || se[j] <= 0) {
                inv_link(eta[j])
              } else {
                inv_link(rnorm(n_sim, eta[j], se[j]))
              }
            })
          )
          as.numeric(stats::quantile(sim_vals, 0.975, na.rm = TRUE))
        },
        
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
    p_detect = summarise_link(p_detect_list[keep], "identity", "p_detect"),
    
    psi_list = psi_list[keep],
    capture_list = capture_list[keep],
    lambda_list = lambda_list[keep],
    p_detect_list = p_detect_list[keep],
    
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
      "Hierarchy: site occupancy Z -> biological-sample capture A -> replicate-level read counts Y.",
      "Intervals are based on simulation from eta uncertainty and transformed to the natural scale.",
      "p_detect uncertainty is propagated from the abundance model linear predictor.",
      "AIC values are exploratory component-model diagnostics, not a formal joint-likelihood criterion."
    )
  )
}
