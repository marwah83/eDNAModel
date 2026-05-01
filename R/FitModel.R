#' Fit a Hierarchical Occupancy–Detection–Abundance Model for Microbiome Data
#'
#' Implements an EM-like iterative framework using repeated GLMM fitting
#' to jointly estimate microbial occupancy (presence/absence), detection
#' probability, and abundance for OTU-level microbiome data stored in a
#' `phyloseq` object.
#'
#' @param phyloseq A `phyloseq` object containing OTU table, sample data, and taxonomy.
#' @param site_col Character. Name of the column representing sampling sites.
#' @param sample_col Character. Name of the column representing samples.
#' @param replicate_col Character (optional). Name of the column representing replicates.
#' @param occupancy_formula Formula for occupancy model (must use response `z_sim`).
#' @param capture_formula Formula for detection model (must use response `a_sim`).
#' @param abundance_formula Formula for abundance model (must use response matching `count_col`).
#' @param otu_col Character. Name of the OTU column.
#' @param count_col Character. Name of the count column.
#' @param min_species_sum Integer. Minimum total count for an OTU to be retained. Default: 10.
#' @param min_detection_replicates Integer. Minimum number of detections per OTU. Default: 1.
#' @param abundance_threshold Integer. Threshold to define presence. Default: 0.
#' @param n_iter Integer. Number of EM-like iterations. Default: 50.
#' @param burn_in Integer. Number of burn-in iterations to discard. Default: 10.
#' @param abundance_family Character. One of `"poisson"`, `"nbinom"`, `"zip"`, `"zinb"`. Default: `"poisson"`.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{psi}{Posterior summaries of occupancy probabilities}
#'   \item{capture}{Posterior summaries of detection probabilities}
#'   \item{lambda}{Posterior summaries of abundance}
#'   \item{p_detect}{Posterior summaries of detection probability implied by abundance}
#'   \item{psi_list}{List of occupancy linear predictors per iteration}
#'   \item{capture_list}{List of detection linear predictors per iteration}
#'   \item{lambda_list}{List of abundance linear predictors per iteration}
#'   \item{p_detect_list}{List of detection probabilities from abundance model}
#'   \item{occupancy_models}{Fitted occupancy GLMMs}
#'   \item{capture_models}{Fitted detection GLMMs}
#'   \item{abundance_models}{Fitted abundance GLMMs}
#'   \item{site_data}{Final site-level latent data}
#'   \item{sample_data}{Final sample-level latent data}
#'   \item{long_df}{Processed long-format data}
#'   \item{filter_summary}{Summary of OTU filtering}
#'   \item{diagnostic_AIC}{AIC values across iterations}
#' }
#'
#' @details
#' This function decomposes microbial occurrence into three hierarchical processes:
#'
#' 1. **Occupancy** (\eqn{\psi}): probability that an OTU is present at a site,
#'    modeled via a binomial GLMM.
#'
#' 2. **Detection (capture)** (\eqn{p}): probability of detecting an OTU given presence,
#'    modeled via a binomial GLMM conditional on occupancy.
#'
#' 3. **Abundance** (\eqn{\lambda}): counts conditional on detection,
#'    modeled using Poisson, Negative Binomial, or zero-inflated variants.
#'
#' The model is fitted using an EM-like iterative algorithm:
#' \enumerate{
#'   \item Fit occupancy model to latent variable \eqn{z_sim}.
#'   \item Fit detection model conditional on \eqn{z_sim = 1}.
#'   \item Fit abundance model conditional on \eqn{a_sim = 1}.
#'   \item Update latent detection \eqn{a_sim} using abundance predictions.
#'   \item Update latent occupancy \eqn{z_sim} using detection and abundance.
#' }
#'
#' Posterior summaries are computed from iterations after burn-in.
#'
#' @section Model Features:
#' \itemize{
#'   \item Supports complex random-effects structures (e.g., `(1 | OTU)`, `(1 | Sample / OTU)`).
#'   \item Automatically handles sequencing depth via `offset(log(total_reads))`.
#'   \item Robust filtering of low-abundance OTUs.
#'   \item EM-like latent state updates for occupancy and detection.
#'   \item Compatible with zero-inflated count models.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item AIC values are provided for diagnostic purposes only and do not correspond
#'         to a joint likelihood across all model components.
#'   \item Overly complex random-effects structures may lead to convergence issues.
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
#' head(fit$lambda)
#' }
#'
#' @importFrom glmmTMB glmmTMB nbinom2
#' @importFrom stats plogis predict rbinom quantile median
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
  
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    min_species_sum = 0,
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
  long_df[[otu_col]]    <- as.character(long_df[[otu_col]])
  long_df[[sample_col]] <- as.character(long_df[[sample_col]])
  long_df[[count_col]]  <- as.numeric(long_df[[count_col]])
  
  if (!is.null(replicate_col) && replicate_col %in% names(long_df)) {
    long_df[[replicate_col]] <- as.character(long_df[[replicate_col]])
  }
  
  long_df <- long_df[!is.na(long_df[[count_col]]), ]
  
  # ------------------------------------------------------------
  # Robust total_reads from cleaned long data
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
    if (verbose) message("Offset detected → cleaning zero-read samples...")
    
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
  
  if (nrow(long_df) == 0) stop("No OTUs remain after filtering.")
  
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
  
  site_keys <- c(site_col, otu_col)
  sample_keys <- c(site_col, sample_col, otu_col)
  
  # ------------------------------------------------------------
  # Build latent state data
  # ------------------------------------------------------------
  
  site_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
    dplyr::summarise(
      z_obs = as.integer(any(.data[[count_col]] > abundance_threshold)),
      .groups = "drop"
    ) |>
    dplyr::mutate(z_sim = z_obs)
  
  sample_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
    dplyr::summarise(
      a_obs = as.integer(any(.data[[count_col]] > abundance_threshold)),
      .groups = "drop"
    ) |>
    dplyr::mutate(a_sim = a_obs)
  
  # ------------------------------------------------------------
  # Validate formulas
  # ------------------------------------------------------------
  
  validate_formula <- function(formula, data, expected_response) {
    if (!inherits(formula, "formula")) stop("Formula must be a valid formula.")
    
    response <- all.vars(formula)[1]
    
    if (response != expected_response) {
      stop(
        "Formula response must be '", expected_response,
        "', but got '", response, "'."
      )
    }
    
    vars <- all.vars(formula)
    missing_vars <- setdiff(vars, names(data))
    
    if (length(missing_vars) > 0) {
      stop("Missing variables in formula: ", paste(missing_vars, collapse = ", "))
    }
    
    invisible(TRUE)
  }
  
  validate_formula(occupancy_formula, site_data, "z_sim")
  validate_formula(capture_formula, sample_data, "a_sim")
  validate_formula(abundance_formula, long_df, count_col)
  
  # ------------------------------------------------------------
  # Family setup
  # ------------------------------------------------------------
  
  fam <- switch(
    abundance_family,
    poisson = poisson(),
    nbinom  = glmmTMB::nbinom2(),
    zip     = poisson(),
    zinb    = glmmTMB::nbinom2()
  )
  
  zi_formula <- if (abundance_family %in% c("zip", "zinb")) ~1 else ~0
  
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
      family = binomial()
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
    # 2. Capture model conditional on z_sim = 1
    # ----------------------------------------------------------
    
    sample_data <- sample_data |>
      dplyr::select(-dplyr::any_of("z_sim")) |>
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
      family = binomial()
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
    # 3. Abundance model conditional on a_sim = 1
    # ----------------------------------------------------------
    
    abund_data <- long_df |>
      dplyr::left_join(
        sample_data |>
          dplyr::select(dplyr::all_of(c(sample_keys, "a_sim"))),
        by = sample_keys
      ) |>
      dplyr::filter(.data$a_sim == 1)
    
    if (nrow(abund_data) == 0) {
      stop("No rows available for abundance model at iteration ", i)
    }
    
    abund_fit <- glmmTMB::glmmTMB(
      formula = abundance_formula,
      data = abund_data,
      family = fam,
      ziformula = zi_formula
    )
    
    eta_lambda_abund <- as.numeric(predict(
      abund_fit,
      type = "link",
      newdata = abund_data,
      allow.new.levels = TRUE
    ))
    
    lambda_abund <- exp(eta_lambda_abund)
    
    lambda_list[[i]] <- data.frame(
      abund_data[sample_keys],
      eta = eta_lambda_abund
    )
    
    # implied P(Y > 0) = 1 - exp(-lambda)
    p_detect_eta <- log(-log(pmax(exp(-lambda_abund), 1e-12)))
    
    p_detect_list[[i]] <- data.frame(
      abund_data[sample_keys],
      eta = p_detect_eta
    )
    
    abundance_models[[i]] <- abund_fit
    
    # ----------------------------------------------------------
    # 4. Update a_sim using abundance feedback
    # ----------------------------------------------------------
    
    sample_prob <- abund_data |>
      dplyr::mutate(lambda_tmp = lambda_abund) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
      dplyr::summarise(
        p0_count = prod(exp(-.data$lambda_tmp)),
        .groups = "drop"
      )
    
    sample_data <- sample_data |>
      dplyr::select(-dplyr::any_of(c("p0_count", "capture_prob"))) |>
      dplyr::mutate(capture_prob = capture_prob) |>
      dplyr::left_join(sample_prob, by = sample_keys)
    
    sample_data$p0_count[is.na(sample_data$p0_count)] <- 1
    
    zero_sample <- which(sample_data$a_obs == 0 & sample_data$z_sim == 1)
    
    if (length(zero_sample) > 0) {
      
      numerator_a <- sample_data$capture_prob[zero_sample] *
        sample_data$p0_count[zero_sample]
      
      denominator_a <- (1 - sample_data$capture_prob[zero_sample]) +
        sample_data$capture_prob[zero_sample] *
        sample_data$p0_count[zero_sample]
      
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
    # 5. Update z_sim using capture + abundance feedback
    # ----------------------------------------------------------
    
    site_prob <- sample_data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(
        p0_site = prod(
          (1 - .data$capture_prob) + .data$capture_prob * .data$p0_count,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    site_data <- site_data |>
      dplyr::select(-dplyr::any_of("p0_site")) |>
      dplyr::left_join(site_prob, by = site_keys)
    
    site_data$p0_site[is.na(site_data$p0_site)] <- 1
    
    zero_site <- which(site_data$z_obs == 0)
    
    if (length(zero_site) > 0) {
      
      numerator_z <- psi[zero_site] * site_data$p0_site[zero_site]
      
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
        occupancy_AIC = AIC(occ_fit),
        capture_AIC = AIC(cap_fit),
        abundance_AIC = AIC(abund_fit)
      )
    )
  }
  
  # ------------------------------------------------------------
  # Summaries on probability / response scale
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
        mean   = mean(.data$value, na.rm = TRUE),
        median = stats::median(.data$value, na.rm = TRUE),
        lwr    = stats::quantile(.data$value, 0.025, na.rm = TRUE),
        upr    = stats::quantile(.data$value, 0.975, na.rm = TRUE),
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
      "Approximate iterative GLMM with EM-like latent-state updates.",
      "AIC values are exploratory for component models, not a formal joint-likelihood criterion."
    )
  )
}
