#' Fit a Multispecies 3‑Level Occupancy–Detection–Abundance Model using GLLVM and GLMM
#'
#' This function fits a hierarchical three‑level model for eDNA, microbiome, or other
#' community count data stored as a \code{phyloseq} object. The model explicitly
#' separates the ecological processes of occupancy, capture, and abundance.
#'
#' \itemize{
#'   \item \strong{Site level (Z)} — true presence/absence of each taxon
#'   \item \strong{Biological sample level (A)} — intermediate capture/detection process
#'   \item \strong{Replicate level (Y)} — observed read counts or replicate‑level measurements
#' }
#'
#' Occupancy is modelled using a binomial Generalized Linear Latent Variable Model (GLLVM),
#' while capture and abundance are modelled using GLMMs via \code{glmmTMB}.
#'
#' @param phyloseq A \code{phyloseq} object containing taxa/OTU count data and metadata.
#' @param site_col Character name of the site or location identifier at the occupancy level.
#' @param species_col Character name of the species/taxon column (formerly hardcoded as \code{"OTU"}).
#' @param response_col Character name of the abundance or response variable
#'   (formerly hardcoded as \code{"y"} in count data).
#' @param sample_col Character name of the biological sample identifier (capture level).
#' @param replicate_col Character name or \code{NULL}; identifies technical replicates
#'   (PCR or sequencing replicates) for the abundance model.
#' @param abundance_rhs Right‑hand‑side expression for the abundance model, e.g.
#'   \code{(1 | species)} or \code{offset(log(total_reads)) + (1 | species)}.
#' @param capture_formula Formula for the sample‑level capture model;
#'   its response must be \code{a_sim}.
#' @param occupancy_covars Optional character vector naming site‑level covariates
#'   to include in the GLLVM occupancy model.
#' @param min_species_sum Minimum total count required for a species (filtering threshold).
#' @param min_detection_replicates Minimum number of detections required to retain the species.
#' @param abundance_threshold Integer threshold defining a positive detection (default = 1).
#' @param n_iter Number of EM‑like iterations of the hierarchical fitting procedure.
#' @param burn_in Number of initial iterations to discard as burn‑in.
#' @param abundance_family Distribution family for the abundance GLMM;
#'   one of \code{"poisson"}, \code{"nbinom"}, \code{"zip"}, or \code{"zinb"}.
#' @param num_lv_c Number of latent variables (\code{num.lv}) in the GLLVM occupancy model.
#' @param verbose Logical; if \code{TRUE}, print progress at each iteration.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{Posterior summaries at the site–species level, including
#'     occupancy (\code{psi_*}), capture (\code{capture_*}),
#'     abundance (\code{lambda_*}), and detection (\code{p_detect_*}).}
#'   \item{capture}{Sample‑level capture summaries.}
#'   \item{capture_site}{Capture summaries aggregated by site.}
#'   \item{psi_list, capture_list, lambda_list, p_detect_list}{Stored linear predictors per iteration.}
#'   \item{occupancy_models}{List of fitted GLLVM models at the site‑level.}
#'   \item{capture_models}{List of fitted capture GLMMs.}
#'   \item{abundance_models}{List of fitted abundance GLMMs.}
#'   \item{reduced_data, sample_data, long_df}{Processed data used in each hierarchical stage.}
#'   \item{lv_sites, lv_species}{Latent variable coordinates from the GLLVM.}
#'   \item{mean_lv_sites, mean_lv_species}{Posterior mean latent variable scores.}
#'   \item{filter_summary}{Information on species retained after filtering.}
#' }
#'
#' @details
#' The model represents a three‑level hierarchy:
#'
#' \deqn{
#' Z_{i,m} \sim \mathrm{Bernoulli}(\psi_{i,m})
#' }
#' \deqn{
#' A_{i,j,m} \mid Z_{i,m} \sim \mathrm{Bernoulli}(p_{i,j,m} \cdot Z_{i,m})
#' }
#' \deqn{
#' Y_{i,j,k,m} \mid A_{i,j,m} \sim \mathrm{Count}(\lambda_{i,j,k,m} \cdot A_{i,j,m})
#' }
#'
#' where:
#' \itemize{
#'   \item \(i\): site,
#'   \item \(j\): biological sample,
#'   \item \(k\): technical replicate,
#'   \item \(m\): species (from \code{species_col}).
#' }
#'
#' The observed occupancy indicator is defined as
#' \deqn{z^{obs}_{i,m} = I\left(\max_{j,k} Y_{i,j,k,m} > c\right)}
#' where \(c =\) \code{abundance_threshold}.
#'
#' @section EM‑like algorithm:
#' Each iteration performs:
#' \enumerate{
#'   \item Fit GLLVM occupancy model for \(Z\)
#'   \item Fit GLMM capture model for \(A \mid Z = 1\)
#'   \item Fit GLMM abundance model for \(Y \mid A = 1\)
#'   \item Update \(A\) and \(Z\) using posterior probabilities
#' }
#'
#' @section Model features:
#' \itemize{
#'   \item Fully parameterized: \code{species_col} and \code{response_col} are user‑defined
#'   \item Three‑level hierarchical structure (Z → A → Y)
#'   \item Latent ecological gradients captured via GLLVM
#'   \item Flexible abundance and detection families (\code{Poisson}, \code{NB}, \code{ZIP}, \code{ZINB})
#'   \item Support for offsets such as \code{offset(log(total_reads))}
#' }
#'
#' @section Caveats:
#' \itemize{
#'   \item EM‑like approximation; not full joint likelihood
#'   \item Detection (\code{p_detect}) depends on abundance model assumptions
#'   \item Large \(\lambda\) may saturate detection near 1
#' }
#'
#' @examples
#' \dontrun{
#' fit <- FitModel_gllvm(
#'   phyloseq = ps,
#'   site_col = "site_month",
#'   species_col = "Taxon",
#'   response_col = "count",
#'   sample_col = "SampleID",
#'   replicate_col = "Replicate",
#'   abundance_rhs = (1 | Taxon),
#'   capture_formula = a_sim ~ 1 + (1 | Taxon),
#'   abundance_family = "nbinom",
#'   n_iter = 10, burn_in = 2
#' )
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
  species_col = "OTU",        # replaces hardcoded 'OTU'
  response_col = "y",         # replaces hardcoded 'y'
  sample_col = "Name",
  replicate_col = NULL,
  abundance_rhs,
  capture_formula = a_sim ~ 1 + (1 | species_col),
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

  if (burn_in >= n_iter) stop("burn_in must be smaller than n_iter.")

  # ------------------------------------------------------------
  # Handle abundance formula
  # ------------------------------------------------------------
  abundance_rhs_expr <- substitute(abundance_rhs)
  abundance_rhs_txt  <- paste(deparse(abundance_rhs_expr), collapse = " ")

  abundance_formula <- stats::as.formula(
    paste(response_col, "~", abundance_rhs_txt),
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
  # Prepare long-format data
  # ------------------------------------------------------------
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col,
    nested_cols = unique(c(sample_col, replicate_col))
  )
  long_df <- prep$long_df

  needed_cols <- c(site_col, sample_col, species_col, response_col)
  missing_cols <- setdiff(needed_cols, names(long_df))
  if (length(missing_cols) > 0)
    stop("Missing columns in long_df: ", paste(missing_cols, collapse = ", "))

  # ------------------------------------------------------------
  # Type handling
  # ------------------------------------------------------------
  long_df[[site_col]]    <- as.character(long_df[[site_col]])
  long_df[[sample_col]]  <- as.character(long_df[[sample_col]])
  long_df[[species_col]] <- as.character(long_df[[species_col]])
  long_df[[response_col]] <- as.numeric(long_df[[response_col]])

  if (!is.null(replicate_col) && replicate_col %in% names(long_df))
    long_df[[replicate_col]] <- as.character(long_df[[replicate_col]])

  long_df <- long_df[!is.na(long_df[[response_col]]), ]

  # ------------------------------------------------------------
  # Compute total_reads
  # ------------------------------------------------------------
  total_reads_df <- long_df |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      total_reads = sum(.data[[response_col]], na.rm = TRUE),
      .groups = "drop"
    )
  long_df <- dplyr::left_join(long_df, total_reads_df, by = sample_col)

  uses_offset <- any(grepl("offset", deparse(abundance_formula)))
  if (uses_offset) {
    zero_samples <- unique(long_df[[sample_col]][long_df$total_reads <= 0])
    if (length(zero_samples) > 0)
      long_df <- dplyr::filter(long_df, !(.data[[sample_col]] %in% zero_samples))
    if (any(long_df$total_reads <= 0, na.rm = TRUE))
      stop("total_reads must be positive when using offset(log(total_reads)).")
  }

  # ------------------------------------------------------------
  # Filter low-abundance species
  # ------------------------------------------------------------
  sp_stats <- long_df |>
    dplyr::group_by(.data[[species_col]]) |>
    dplyr::summarise(
      total_count = sum(.data[[response_col]], na.rm = TRUE),
      detected_replicates = sum(.data[[response_col]] > abundance_threshold, na.rm = TRUE),
      .groups = "drop"
    )

  keep_species <- sp_stats |>
    dplyr::filter(
      .data$total_count >= min_species_sum,
      .data$detected_replicates >= min_detection_replicates
    ) |>
    dplyr::pull(.data[[species_col]])

  long_df <- long_df |>
    dplyr::filter(.data[[species_col]] %in% keep_species)

  if (nrow(long_df) == 0)
    stop("No species remain after filtering.")

  if (verbose) {
    message("Species before filtering: ", dplyr::n_distinct(sp_stats[[species_col]]))
    message("Species after filtering: ", dplyr::n_distinct(long_df[[species_col]]))
  }

  # ------------------------------------------------------------
  # Factorize columns
  # ------------------------------------------------------------
  factor_cols <- unique(c(
    site_col, sample_col, replicate_col, species_col,
    occupancy_covars, all.vars(capture_formula), all.vars(abundance_formula)
  ))
  factor_cols <- setdiff(factor_cols, c(response_col, "total_reads", "z_sim", "a_sim"))
  long_df <- to_factor_cols(long_df, factor_cols)

  # ------------------------------------------------------------
  # Reference species for relevel
  # ------------------------------------------------------------
  sp_abundance <- long_df |>
    dplyr::group_by(.data[[species_col]]) |>
    dplyr::summarise(total_count = sum(.data[[response_col]], na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$total_count))
  top_species <- sp_abundance[[species_col]][1]
  long_df[[species_col]] <- stats::relevel(factor(long_df[[species_col]]), ref = as.character(top_species))

  # ------------------------------------------------------------
  # Define key identifiers
  # ------------------------------------------------------------
  site_keys <- c(site_col, species_col)
  sample_keys <- c(site_col, sample_col, species_col)
  pcr_keys <- if (!is.null(replicate_col) && replicate_col %in% names(long_df))
    c(site_col, sample_col, replicate_col, species_col)
  else
    c(site_col, sample_col, species_col)

  # ------------------------------------------------------------
  # Build site-level occupancy (Z)
  # ------------------------------------------------------------
  reduced_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
    dplyr::summarise(
      z_obs = as.integer(any(.data[[response_col]] > abundance_threshold)),
      dplyr::across(-dplyr::all_of(response_col), dplyr::first),
      .groups = "drop"
    ) |>
    dplyr::mutate(z_sim = .data$z_obs)

  reduced_data <- to_factor_cols(reduced_data, factor_cols)

  # ------------------------------------------------------------
  # Build biological sample-level capture data (A)
  # ------------------------------------------------------------
  cap_vars <- get_formula_vars(capture_formula, "a_sim")
  sample_keep_vars <- setdiff(cap_vars, c(sample_keys, "a_sim"))
  sample_keep_vars <- intersect(sample_keep_vars, names(long_df))

  sample_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
    dplyr::summarise(
      a_obs = as.integer(any(.data[[response_col]] > abundance_threshold)),
      dplyr::across(all_of(sample_keep_vars), dplyr::first),
      .groups = "drop"
    ) |>
    dplyr::mutate(a_sim = .data$a_obs)
  sample_data <- to_factor_cols(sample_data, factor_cols)

  # ------------------------------------------------------------
  # Validate formulas
  # ------------------------------------------------------------
  validate_formula <- function(formula, data, expected_response, model_name) {
    response <- all.vars(formula)[1]
    if (response != expected_response)
      stop(model_name, " formula response must be '", expected_response, "'. Got '", response, "'.")
    missing_vars <- setdiff(all.vars(formula), names(data))
    if (length(missing_vars) > 0)
      stop("Missing variables in ", model_name, " formula: ", paste(missing_vars, collapse = ", "))
  }
  validate_formula(capture_formula, sample_data, "a_sim", "capture")
  validate_formula(abundance_formula, long_df, response_col, "abundance")

  # ------------------------------------------------------------
  # Family definitions
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
  # Helper: zero probability
  # ------------------------------------------------------------
  get_p0_abundance <- function(fit, newdata, family_name) {
    mu <- as.numeric(stats::predict(fit, type = "response", newdata = newdata, allow.new.levels = TRUE))
    zi <- rep(0, length(mu))
    if (family_name %in% c("zip", "zinb")) {
      zi <- tryCatch(
        as.numeric(stats::predict(fit, type = "zprob", newdata = newdata, allow.new.levels = TRUE)),
        error = function(e) rep(0, length(mu))
      )
    }
    p0_cond <- switch(
      family_name,
      poisson = stats::dpois(0, lambda = mu),
      zip     = stats::dpois(0, lambda = mu),
      nbinom  = {
        theta <- tryCatch(stats::sigma(fit), error = function(e) NA_real_)
        if (!is.finite(theta) || theta <= 0) stats::dpois(0, lambda = mu)
        else stats::dnbinom(0, mu = mu, size = theta)
      },
      zinb = {
        theta <- tryCatch(stats::sigma(fit), error = function(e) NA_real_)
        if (!is.finite(theta) || theta <= 0) stats::dpois(0, lambda = mu)
        else stats::dnbinom(0, mu = mu, size = theta)
      }
    )
    p0 <- zi + (1 - zi) * p0_cond
    pmin(pmax(p0, 1e-12), 1)
  }

  # ------------------------------------------------------------
  # Storage
  # ------------------------------------------------------------
  psi_list <- capture_list <- lambda_list <- p_detect_list <- list()
  occupancy_models <- capture_models <- abundance_models <- list()
  lv_sites_list <- lv_species_list <- list()

  # ------------------------------------------------------------
  # EM-like iterations
  # ------------------------------------------------------------
  for (i in seq_len(n_iter)) {
    if (verbose) message("Iteration ", i)

    reduced_data[[site_col]] <- as.character(reduced_data[[site_col]])
    long_df[[site_col]]      <- as.character(long_df[[site_col]])

    # 1. Occupancy model ---------------------------------------
    z_matrix <- reshape2::acast(
      reduced_data,
      stats::as.formula(paste(site_col, "~", species_col)),
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
    X_cov <- if (!is.null(occupancy_covars) && length(occupancy_covars) > 0)
      stats::model.matrix(~ ., data = cov_df[, occupancy_covars, drop = FALSE])[, -1, drop = FALSE]
    else NULL

    model_occupancy <- gllvm::gllvm(y = z_matrix, X = X_cov, family = "binomial", num.lv = num_lv_c)
    occupancy_models[[i]] <- model_occupancy

    # latent vars
    lv_sites <- as.data.frame(model_occupancy$lvs)
    colnames(lv_sites) <- paste0("LV", seq_len(ncol(lv_sites)))
    lv_sites$Iteration <- i
    lv_sites$Site      <- rownames(model_occupancy$lvs)
    lv_species <- as.data.frame(model_occupancy$params$theta)
    colnames(lv_species) <- paste0("LV", seq_len(ncol(lv_species)))
    lv_species$Iteration <- i
    lv_species[[species_col]] <- rownames(model_occupancy$params$theta)
    lv_sites_list[[i]]   <- lv_sites
    lv_species_list[[i]] <- lv_species

    psi_prob <- stats::predict(model_occupancy, type = "response")
    rownames(psi_prob) <- rownames(z_matrix)
    psi_long <- reshape2::melt(psi_prob, varnames = c(site_col, species_col), value.name = "psi_prob")
    psi_long$psi_prob <- pmin(pmax(psi_long$psi_prob, 1e-6), 1 - 1e-6)
    psi_long$eta <- stats::qlogis(psi_long$psi_prob)
    psi_list[[i]] <- psi_long[, c(site_col, species_col, "eta")]

    # 2. Capture model -----------------------------------------
    sample_data <- sample_data |>
      dplyr::select(-dplyr::any_of(c("z_sim", "p0_sample", "capture_prob"))) |>
      dplyr::left_join(
        reduced_data |> dplyr::select(dplyr::all_of(c(site_col, species_col, "z_sim"))),
        by = c(site_col, species_col)
      )
    capture_fit_data <- dplyr::filter(sample_data, .data$z_sim == 1)
    if (nrow(capture_fit_data) == 0)
      stop("No rows available for capture model at iteration ", i)

    cap_fit <- glmmTMB::glmmTMB(formula = capture_formula,
                                data = capture_fit_data,
                                family = stats::binomial())
    capture_models[[i]] <- cap_fit

    eta_capture <- as.numeric(stats::predict(cap_fit, type = "link",
                                             newdata = sample_data, allow.new.levels = TRUE))
    capture_prob <- stats::plogis(eta_capture)

    capture_list[[i]] <- data.frame(sample_data[sample_keys], eta = eta_capture)

    # 3. Abundance model ---------------------------------------
    abundance_data <- long_df |>
      dplyr::left_join(sample_data |> dplyr::select(dplyr::all_of(c(sample_keys, "a_sim"))),
                       by = sample_keys) |>
      dplyr::filter(.data$a_sim == 1)
    if (nrow(abundance_data) == 0)
      stop("No rows available for abundance model at iteration ", i)

    abundance_data <- to_factor_cols(abundance_data, factor_cols)
    model_abundance <- glmmTMB::glmmTMB(formula = abundance_formula,
                                        data = abundance_data,
                                        family = abundance_glmm_family,
                                        ziformula = zi_formula)
    abundance_models[[i]] <- model_abundance

    lambda_eta <- as.numeric(stats::predict(model_abundance, type = "link",
                                            newdata = long_df, allow.new.levels = TRUE))
    lambda_list[[i]] <- data.frame(long_df[site_keys], eta = lambda_eta) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(eta = mean(.data$eta, na.rm = TRUE), .groups = "drop")

    # 4. p0 estimates ------------------------------------------
    pred_data_all <- long_df |>
      dplyr::left_join(sample_data |> dplyr::select(dplyr::all_of(c(sample_keys, "a_sim"))),
                       by = sample_keys)
    pred_data_all <- to_factor_cols(pred_data_all, factor_cols)
    p0_pcr <- get_p0_abundance(model_abundance, pred_data_all, abundance_family)
    sample_p0 <- pred_data_all |>
      dplyr::mutate(p0_pcr = p0_pcr) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(sample_keys))) |>
      dplyr::summarise(p0_sample = prod(.data$p0_pcr, na.rm = TRUE), .groups = "drop")

    sample_data <- sample_data |>
      dplyr::select(-dplyr::any_of(c("p0_sample", "capture_prob"))) |>
      dplyr::mutate(capture_prob = capture_prob) |>
      dplyr::left_join(sample_p0, by = sample_keys)
    sample_data$p0_sample[is.na(sample_data$p0_sample)] <- 1

    # latent A update
    zero_sample <- which(sample_data$a_obs == 0 & sample_data$z_sim == 1)
    if (length(zero_sample) > 0) {
      numerator_a <- sample_data$capture_prob[zero_sample] * sample_data$p0_sample[zero_sample]
      denominator_a <- (1 - sample_data$capture_prob[zero_sample]) +
        sample_data$capture_prob[zero_sample] * sample_data$p0_sample[zero_sample]
      posterior_a <- numerator_a / pmax(denominator_a, 1e-12)
      posterior_a <- pmin(pmax(posterior_a, 0.001), 0.999)
      sample_data$a_sim[zero_sample] <- stats::rbinom(length(zero_sample), 1, posterior_a)
    }
    sample_data$a_sim[sample_data$a_obs == 1] <- 1
    sample_data$a_sim[sample_data$z_sim == 0] <- 0

    # 5. Collapse to site level --------------------------------
    site_p0 <- sample_data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(site_keys))) |>
      dplyr::summarise(p0_site = prod((1 - .data$capture_prob) +
                                         .data$capture_prob * .data$p0_sample,
                                      na.rm = TRUE),
                       .groups = "drop")
    p_detect_list[[i]] <- site_p0 |>
      dplyr::mutate(eta = log(-log(pmax(.data$p0_site, 1e-12)))) |>
      dplyr::select(dplyr::all_of(site_keys), eta)

    # 6. Update latent Z ---------------------------------------
    reduced_data <- reduced_data |>
      dplyr::select(-dplyr::any_of("p0_site")) |>
      dplyr::left_join(site_p0, by = site_keys) |>
      dplyr::mutate(p0_site = tidyr::replace_na(.data$p0_site, 1))
    z_merge <- reduced_data |>
      dplyr::left_join(psi_long[, c(site_col, species_col, "psi_prob")], by = site_keys)
    zero_idx <- which(z_merge$z_obs == 0)
    if (length(zero_idx) > 0) {
      num_z <- z_merge$psi_prob[zero_idx] * z_merge$p0_site[zero_idx]
      den_z <- (1 - z_merge$psi_prob[zero_idx]) +
        z_merge$psi_prob[zero_idx] * z_merge$p0_site[zero_idx]
      post_z <- num_z / pmax(den_z, 1e-12)
      post_z <- pmin(pmax(post_z, 0.001), 0.999)
      reduced_data$z_sim[zero_idx] <- stats::rbinom(length(zero_idx), 1, post_z)
    }
    reduced_data$z_sim[reduced_data$z_obs == 1] <- 1
  }

  # ------------------------------------------------------------
  # Summaries
  # ------------------------------------------------------------
  summarise_link <- function(lst, link_type = c("logit", "log", "cloglog")) {
    link_type <- match.arg(link_type)
    df <- dplyr::bind_rows(lst)
    if (nrow(df) == 0) return(data.frame())
    inv_link <- switch(link_type,
                       logit = stats::plogis,
                       log = exp,
                       cloglog = function(x) 1 - exp(-exp(x)))
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
    if (nrow(df) == 0) return(data.frame())
    df |> dplyr::mutate(capture = stats::plogis(.data$eta)) |>
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
    if (nrow(capture_summary) == 0) return(data.frame())
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

  keep <- seq.int(burn_in + 1, n_iter)
  psi_summary <- summarise_link(psi_list[keep], "logit") |>
    dplyr::rename_with(~paste0("psi_", .x), -dplyr::all_of(site_keys))
  lambda_summary <- summarise_link(lambda_list[keep], "log") |>
    dplyr::rename_with(~paste0("lambda_", .x), -dplyr::all_of(site_keys))
  p_detect_summary <- summarise_link(p_detect_list[keep], "cloglog") |>
    dplyr::rename_with(~paste0("p_detect_", .x), -dplyr::all_of(site_keys))
  capture_summary <- summarise_capture(capture_list[keep])
  capture_site_summary <- summarise_capture_site(capture_summary)
  final_summary <- psi_summary |>
    dplyr::left_join(capture_site_summary, by = site_keys) |>
    dplyr::left_join(lambda_summary, by = site_keys) |>
    dplyr::left_join(p_detect_summary, by = site_keys)

  # ------------------------------------------------------------
  # Latent variable summaries
  # ------------------------------------------------------------
  lv_sites_combined <- dplyr::bind_rows(lv_sites_list[keep])
  lv_species_combined <- dplyr::bind_rows(lv_species_list[keep])
  mean_lv_sites <- lv_sites_combined |>
    dplyr::group_by(.data$Site) |>
    dplyr::summarise(dplyr::across(dplyr::starts_with("LV"), mean, na.rm = TRUE), .groups = "drop")
  mean_lv_species <- lv_species_combined |>
    dplyr::group_by(.data[[species_col]]) |>
    dplyr::summarise(dplyr::across(dplyr::starts_with("LV"), mean, na.rm = TRUE), .groups = "drop")

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
      species_stats = sp_stats,
      kept_species = keep_species,
      min_species_sum = min_species_sum,
      min_detection_replicates = min_detection_replicates,
      abundance_threshold = abundance_threshold
    ),
    note = paste(
      "Three-level flexible GLLVM–GLMM model.",
      "Hierarchy: site occupancy (Z) via gllvm, sample capture (A) via GLMM,",
      "replicate-level abundance (Y) via GLMM; species_col and response_col are user-defined."
    )
  )
}
