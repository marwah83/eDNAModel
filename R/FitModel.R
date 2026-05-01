#' Fit a Hierarchical Occupancy-Detection-Abundance Model for Microbiome Data
#'
#' Implements an EM-like iterative framework using repeated GLMM fitting
#' to jointly estimate microbial occupancy, detection probability, and abundance
#' for OTU-level microbiome data stored in a `phyloseq` object.
#'
#' Formula covariates are automatically preserved in the internally constructed
#' site-level and sample-level data sets. This allows covariates such as
#' `in_out`, `Samplingmonth`, treatment variables, interactions, and nested
#' random-effect terms to be used directly in the model formulas.
#'
#' @param phyloseq A `phyloseq` object containing an OTU table, sample data, and optionally taxonomy.
#' @param site_col Character. Name of the column representing sampling sites.
#' @param sample_col Character. Name of the column representing samples.
#' @param replicate_col Character or `NULL`. Name of the column representing replicates.
#' @param occupancy_formula Formula for the occupancy model. The response must be `z_sim`.
#' @param capture_formula Formula for the detection/capture model. The response must be `a_sim`.
#' @param abundance_formula Formula for the abundance model. The response must match `count_col`.
#' @param otu_col Character. Name of the OTU column.
#' @param count_col Character. Name of the count column.
#' @param min_species_sum Integer. Minimum total count for an OTU to be retained.
#' @param min_detection_replicates Integer. Minimum number of detections required for an OTU to be retained.
#' @param abundance_threshold Numeric. Count threshold above which an OTU is treated as detected.
#' @param n_iter Integer. Number of EM-like iterations.
#' @param burn_in Integer. Number of initial iterations discarded before summarising results.
#' @param abundance_family Character. One of `"poisson"`, `"nbinom"`, `"zip"`, or `"zinb"`.
#' @param verbose Logical. If `TRUE`, print progress messages.
#'
#' @return A list containing:
#' \describe{
#'   \item{psi}{Data frame of occupancy probability summaries with columns such as `psi_mean`, `psi_median`, `psi_lwr`, and `psi_upr`.}
#'   \item{capture}{Data frame of detection/capture probability summaries with columns such as `capture_mean`, `capture_median`, `capture_lwr`, and `capture_upr`.}
#'   \item{lambda}{Data frame of abundance summaries with columns such as `lambda_mean`, `lambda_median`, `lambda_lwr`, and `lambda_upr`.}
#'   \item{p_detect}{Data frame of abundance-implied detection probability summaries with columns such as `p_detect_mean`, `p_detect_median`, `p_detect_lwr`, and `p_detect_upr`.}
#'   \item{psi_list}{List of occupancy linear predictors retained after burn-in.}
#'   \item{capture_list}{List of detection linear predictors retained after burn-in.}
#'   \item{lambda_list}{List of abundance linear predictors retained after burn-in.}
#'   \item{p_detect_list}{List of abundance-implied detection linear predictors retained after burn-in.}
#'   \item{occupancy_models}{List of fitted occupancy `glmmTMB` models for all iterations.}
#'   \item{capture_models}{List of fitted detection/capture `glmmTMB` models for all iterations.}
#'   \item{abundance_models}{List of fitted abundance `glmmTMB` models for all iterations.}
#'   \item{site_data}{Final site-level latent-state data.}
#'   \item{sample_data}{Final sample-level latent-state data.}
#'   \item{long_df}{Processed long-format data used for model fitting.}
#'   \item{filter_summary}{List containing OTU filtering summaries and retained OTUs.}
#'   \item{diagnostic_AIC}{Data frame of component-model AIC values across iterations.}
#'   \item{note}{Character note describing the approximate EM-like fitting procedure.}
#' }
#'
#' @details
#' The model decomposes microbial occurrence into three linked processes:
#'
#' \enumerate{
#'   \item Occupancy: whether an OTU is present at a site, modeled using a binomial GLMM.
#'   \item Detection/capture: whether an OTU is detected in a sample conditional on occupancy, modeled using a binomial GLMM.
#'   \item Abundance: observed counts conditional on detection, modeled using Poisson, negative binomial, or zero-inflated count models.
#' }
#'
#' At each iteration, the function fits the occupancy, detection, and abundance
#' components, then updates the latent detection state `a_sim` and occupancy
#' state `z_sim` using the current fitted values. Summaries are computed from
#' iterations retained after burn-in.
#'
#' If `offset(log(total_reads))` is included in `abundance_formula`, total read
#' depth is computed automatically from the processed long-format data. Samples
#' with non-positive total read depth are removed before fitting the abundance
#' model because `log(0)` is undefined.
#'
#' @section Model features:
#' \itemize{
#'   \item Supports fixed effects, interactions, and random effects in all model components.
#'   \item Automatically preserves formula covariates in site-level and sample-level data.
#'   \item Supports offsets through `offset(log(total_reads))` in the abundance model.
#'   \item Supports Poisson, negative binomial, zero-inflated Poisson, and zero-inflated negative binomial abundance models.
#'   \item Stores fitted component models and AIC diagnostics for each iteration.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item The fitting procedure is an approximate EM-like algorithm, not a full joint maximum-likelihood or fully Bayesian sampler.
#'   \item AIC values are component-model diagnostics and should not be interpreted as formal joint-likelihood criteria for the full hierarchical model.
#'   \item Highly complex random-effect structures may be unstable for sparse microbiome or eDNA count data.
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
#'   capture_formula   = a_sim ~ 1 + (1 | OTU),
#'   abundance_formula = y ~ offset(log(total_reads)) + (1 | OTU),
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
#' @importFrom utils head
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
  
  # ------------------------------------------------------------
  # Helper: extract RHS variables from formula
  # ------------------------------------------------------------
  
  get_formula_vars <- function(formula, response) {
    vars <- all.vars(formula)
    setdiff(vars, response)
  }
  
  # ------------------------------------------------------------
  # Prepare long data
  # ------------------------------------------------------------
  
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
  # Compute total_reads from cleaned long data
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
  # Extract all variables needed by formulas
  # ------------------------------------------------------------
  
  occ_vars   <- get_formula_vars(occupancy_formula, "z_sim")
  cap_vars   <- get_formula_vars(capture_formula, "a_sim")
  abund_vars <- get_formula_vars(abundance_formula, count_col)
  
  all_formula_vars <- unique(c(occ_vars, cap_vars, abund_vars))
  all_formula_vars <- setdiff(all_formula_vars, c(otu_col, count_col, "total_reads"))
  
  missing_formula_vars <- setdiff(all_formula_vars, names(long_df))
  
  if (length(missing_formula_vars) > 0) {
    stop(
      "The following formula variables are missing from long_df: ",
      paste(missing_formula_vars, collapse = ", ")
    )
  }
  
  # Convert character formula variables to factors
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
  
  site_keys <- c(site_col, otu_col)
  sample_keys <- c(site_col, sample_col, otu_col)
  
  # ------------------------------------------------------------
  # Variables to keep in site_data and sample_data
  # ------------------------------------------------------------
  
  site_keep_vars <- unique(setdiff(
    get_formula_vars(occupancy_formula, "z_sim"),
    c(site_keys, "z_sim")
  ))
  site_keep_vars <- intersect(site_keep_vars, names(long_df))
  
  sample_keep_vars <- unique(setdiff(
    get_formula_vars(capture_formula, "a_sim"),
    c(sample_keys, "a_sim")
  ))
  sample_keep_vars <- intersect(sample_keep_vars, names(long_df))
  
  # ------------------------------------------------------------
  # Build site-level data
  # ------------------------------------------------------------
  
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
  
  # ------------------------------------------------------------
  # Build sample-level data
  # ------------------------------------------------------------
  
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
    
    # -----------------------------
    # Occupancy model
    # -----------------------------
    
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
    
    # -----------------------------
    # Capture model conditional on z_sim = 1
    # -----------------------------
    
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
    
    # -----------------------------
    # Abundance model conditional on a_sim = 1
    # -----------------------------
    
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
    
    p_detect_eta <- log(-log(pmax(exp(-lambda_abund), 1e-12)))
    
    p_detect_list[[i]] <- data.frame(
      abund_data[sample_keys],
      eta = p_detect_eta
    )
    
    abundance_models[[i]] <- abund_fit
    
    # -----------------------------
    # Update a_sim
    # -----------------------------
    
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
    
    # -----------------------------
    # Update z_sim
    # -----------------------------
    
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
  # Summaries
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
      "Formula covariates are automatically preserved in site/sample data.",
      "AIC values are exploratory for component models, not a formal joint-likelihood criterion."
    )
  )
}
