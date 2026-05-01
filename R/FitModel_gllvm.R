#' Fit a Multispecies Occupancy–Detection Model using GLLVM and GLMM
#'
#' This function implements a hierarchical Bayesian occupancy–detection model
#' using Generalized Linear Latent Variable Models (GLLVM) for occupancy and 
#' GLMMs for abundance, iterating between them in a Monte Carlo scheme.
#' It is specifically designed for analyzing environmental DNA (eDNA) data
#' structured as a `phyloseq` object.
#'
#' @param phyloseq A `phyloseq` object containing OTU count data, sample metadata, and taxonomy.
#' @param site_col Character string specifying the name of the column in sample metadata to use as the site identifier.
#' @param abundance_rhs A one-sided formula specifying the right-hand side of the abundance model (e.g., `(1 | Samplingmonth/OTU)`).
#' @param occupancy_covars Optional character vector specifying covariates to use in the occupancy model.
#' @param min_species_sum Integer specifying the minimum total abundance required for an OTU to be retained (default = 50).
#' @param abundance_threshold Integer; a threshold of detections above which a taxon is considered present (default = 1).
#' @param n_iter Integer; number of Monte Carlo iterations to run (default = 50).
#' @param burn_in Integer; number of burn-in iterations to discard (default = 10).
#' @param abundance_family Distribution to use for abundance model. One of `"poisson"`, `"nbinom"`, `"zip"`, or `"zinb"`.
#' @param num_lv_c Number of constrained latent variables for the occupancy model. Defaults to 2. Must be <= number of species.
#' 
#' @return A list with the following components:
#' \describe{
#'   \item{summary}{A data frame with posterior summaries of occupancy (`psi_`), abundance (`lambda_`), and detection probability (`p_detect_`).}
#'   \item{psi_list}{List of occupancy probability estimates across iterations.}
#'   \item{lambda_list}{List of abundance estimates across iterations (log-scale).}
#'   \item{p_detect_list}{List of detection probabilities (1 - exp(-lambda)) across iterations.}
#'   \item{occupancy_models}{List of fitted GLLVM models for occupancy.}
#'   \item{abundance_models}{List of fitted GLMMs for abundance.}
#'   \item{reduced_data}{Data frame containing presence/absence data and site/taxon information.}
#'   \item{lv_sites}{Latent variable coordinates for sites from GLLVM across iterations.}
#'   \item{lv_species}{Latent variable coordinates for species (OTUs) from GLLVM across iterations.}
#'   \item{mean_lv_sites}{Averaged site positions in latent variable space (used in biplots).}
#'   \item{mean_lv_species}{Averaged species positions in latent variable space (used in biplots).}
#' }
#' @description
#' This function implements a hierarchical two-part multispecies occupancy model using **gllvm** for modeling species presence (occupancy).
#'
#' The function decomposes species occurrence into two ecological processes:
#'
#' 1. **Occupancy** (presence/absence):
#'    - Modeled using a **binomial Generalized Linear Latent Variable Model (GLLVM)** with optional environmental covariates (`occupancy_covars`) and two common latent variables.
#'    - Latent variables capture residual structure (e.g., species co-occurrence or unmeasured environmental gradients).
#'
#' 2. **Abundance** (counts conditional on presence):
#'    - Modeled using a **Generalized Linear Mixed Model (GLMM)** via `glmmTMB`, allowing for Poisson, Negative Binomial, Zero-Inflated Poisson (ZIP), or Zero-Inflated Negative Binomial (ZINB) families.
#'
#' At each iteration:
#' - Presence data (`z_sim`) is simulated based on detection and occupancy probabilities.
#' - A GLLVM is fitted to estimate occupancy-related parameters.
#' - Conditional abundance is modeled with the specified GLMM formula (`abundance_rhs`).
#' - Detection probabilities and occupancy are updated using posterior predictions.
#'
#' After the burn-in period, posterior summaries (occupancy, abundance, detection probabilities) are computed using **inverse-variance weighted averaging**.
#'
#' In addition, latent variable scores are averaged across iterations to allow **ordination-style visualizations** (biplots) of species and site structure.
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#' # Replace with real site_col and covariates relevant to your dataset
#' physeq <- subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Mock"))
#' out <- FitModel_gllvm(
#' phyloseq = physeq,
#' site_col = "SampleType",
#' abundance_rhs = (1 | OTU),
#' occupancy_covars = c("SampleType"),
#' abundance_family = "poisson",
#' n_iter = 10,
#' burn_in = 2, 
#' num_lv_c = 2
#')
#' head(out$summary)
#' }
#'
#' @seealso \code{\link[gllvm]{gllvm}}, \code{\link[glmmTMB]{glmmTMB}}
#' @import dplyr
#' @import tidyr
#' @importFrom stats reformulate predict rbinom
#' @importFrom gllvm gllvm
#' @importFrom glmmTMB glmmTMB
#' @importFrom phyloseq taxa_sums
#' @importFrom reshape2 melt acast
#' @importFrom utils head
#' @importFrom phyloseq sample_data sample_data<-
#' @export
FitModel_gllvm <- function(
  phyloseq,
  site_col,
  abundance_rhs,
  occupancy_covars = NULL,
  min_species_sum = 50,
  abundance_threshold = 1,
  n_iter = 50,
  burn_in = 10,
  abundance_family = "poisson",
  num_lv_c = 2
) {

  if (burn_in >= n_iter) {
    stop("burn_in must be smaller than n_iter.")
  }

  abundance_rhs_expr <- substitute(abundance_rhs)
  abundance_rhs_txt <- paste(deparse(abundance_rhs_expr), collapse = " ")
  abundance_formula <- stats::as.formula(
    paste("y ~", abundance_rhs_txt),
    env = parent.frame()
  )

  valid_families <- c("poisson", "nbinom", "zip", "zinb")
  if (!(abundance_family %in% valid_families)) {
    stop(
      "'abundance_family' must be one of: ",
      paste(valid_families, collapse = ", ")
    )
  }

  to_factor_cols <- function(df, cols) {
    for (col in cols) {
      if (!is.null(col) && col %in% names(df)) {
        df[[col]] <- as.factor(as.character(df[[col]]))
      }
    }
    df
  }

  # ------------------------------------------------------------
  # Prepare long data: NO filtering here
  # ------------------------------------------------------------
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col,
    nested_cols = NULL
  )

  long_df <- prep$long_df

  needed_cols <- c(site_col, "OTU", "y")
  missing_cols <- setdiff(needed_cols, names(long_df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in long_df: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # ------------------------------------------------------------
  # Explicit OTU filtering inside FitModel_gllvm
  # ------------------------------------------------------------
  otu_stats <- long_df |>
    dplyr::group_by(.data$OTU) |>
    dplyr::summarise(
      total_count = sum(.data$y, na.rm = TRUE),
      detected_replicates = sum(.data$y > abundance_threshold, na.rm = TRUE),
      .groups = "drop"
    )

  keep_otus <- otu_stats |>
    dplyr::filter(.data$total_count >= min_species_sum) |>
    dplyr::pull(.data$OTU)

  long_df <- long_df |>
    dplyr::filter(.data$OTU %in% keep_otus)

  if (nrow(long_df) == 0) {
    stop("No OTUs remain after filtering.")
  }

  message("OTUs before filtering: ", dplyr::n_distinct(otu_stats$OTU))
  message("OTUs after filtering: ", dplyr::n_distinct(long_df$OTU))

  # ------------------------------------------------------------
  # Type handling
  # ------------------------------------------------------------
  factor_cols <- unique(c(site_col, "OTU", "Name", "Samplingmonth"))
  long_df <- to_factor_cols(long_df, factor_cols)

  otu_abundances <- long_df |>
    dplyr::group_by(.data$OTU) |>
    dplyr::summarise(total_count = sum(.data$y, na.rm = TRUE), .groups = "drop")

  top_otu <- otu_abundances |>
    dplyr::arrange(dplyr::desc(.data$total_count)) |>
    dplyr::slice(1) |>
    dplyr::pull(.data$OTU)

  long_df$OTU <- stats::relevel(
    factor(long_df$OTU),
    ref = as.character(top_otu)
  )

  # ------------------------------------------------------------
  # Site × OTU occupancy data
  # ------------------------------------------------------------
  reduced_data <- long_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(site_col, "OTU")))) |>
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
  # Storage
  # ------------------------------------------------------------
  psi_list <- list()
  lambda_list <- list()
  p_detect_list <- list()

  occupancy_models <- list()
  abundance_models <- list()

  lv_sites_list <- list()
  lv_species_list <- list()

  psi_cols <- c(site_col, "OTU")

  # ------------------------------------------------------------
  # Helper: zero probability for abundance model
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
  # Iterations
  # ------------------------------------------------------------
  for (i in seq_len(n_iter)) {

    message("Iteration ", i)

    reduced_data[[site_col]] <- as.character(reduced_data[[site_col]])
    long_df[[site_col]] <- as.character(long_df[[site_col]])

    # ----------------------------------------------------------
    # Occupancy matrix for GLLVM
    # ----------------------------------------------------------
    z_matrix <- reshape2::acast(
      reduced_data,
      stats::as.formula(paste(site_col, "~ OTU")),
      value.var = "z_sim",
      fill = 0
    )

    z_sites <- rownames(z_matrix)

    # ----------------------------------------------------------
    # Occupancy covariates
    # ----------------------------------------------------------
    cov_cols <- unique(c(site_col, occupancy_covars))

    cov_df <- reduced_data |>
      dplyr::select(dplyr::all_of(cov_cols)) |>
      dplyr::distinct()

    cov_df <- to_factor_cols(cov_df, occupancy_covars)

    cov_df <- cov_df[match(z_sites, cov_df[[site_col]]), , drop = FALSE]

    X_cov <- if (!is.null(occupancy_covars) && length(occupancy_covars) > 0) {
      stats::model.matrix(~ ., data = cov_df[, occupancy_covars, drop = FALSE])[, -1, drop = FALSE]
    } else {
      NULL
    }

    # ----------------------------------------------------------
    # GLLVM occupancy model
    # ----------------------------------------------------------
    model_occupancy <- gllvm::gllvm(
      y = z_matrix,
      X = X_cov,
      family = "binomial",
      num.lv.c = num_lv_c
    )

    # ----------------------------------------------------------
    # Store latent variables
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
    # Predict occupancy probabilities
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
    occupancy_models[[i]] <- model_occupancy

    # ----------------------------------------------------------
    # Abundance data conditional on Z = 1
    # ----------------------------------------------------------
    abundance_data <- long_df |>
      dplyr::left_join(
        reduced_data |>
          dplyr::select(dplyr::all_of(c(site_col, "OTU", "z_sim"))),
        by = c(site_col, "OTU")
      ) |>
      dplyr::filter(.data$z_sim == 1)

    if (nrow(abundance_data) == 0) {
      stop("No rows available for abundance model at iteration ", i)
    }

    abundance_data <- to_factor_cols(abundance_data, factor_cols)

    # ----------------------------------------------------------
    # GLMM family
    # ----------------------------------------------------------
    if (abundance_family == "poisson") {
      abundance_glmm_family <- stats::poisson()
      zi_formula <- ~0
    } else if (abundance_family == "nbinom") {
      abundance_glmm_family <- glmmTMB::nbinom2()
      zi_formula <- ~0
    } else if (abundance_family == "zip") {
      abundance_glmm_family <- stats::poisson()
      zi_formula <- ~1
    } else {
      abundance_glmm_family <- glmmTMB::nbinom2()
      zi_formula <- ~1
    }

    # ----------------------------------------------------------
    # Abundance model
    # This keeps nesting, e.g.:
    # y ~ (1 | OTU) + (1 | Name / OTU)
    # ----------------------------------------------------------
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
      se.fit = TRUE,
      allow.new.levels = TRUE
    )

    lambda_eta <- as.numeric(lambda_pred$fit)
    lambda_mu <- exp(lambda_eta)

    p0_abund <- get_p0_abundance(
      fit = model_abundance,
      newdata = abundance_data,
      family_name = abundance_family
    )

    p_detect_eta <- log(-log(pmax(p0_abund, 1e-12)))

    lambda_list[[i]] <- data.frame(
      abundance_data[psi_cols],
      eta = lambda_eta
    )

    p_detect_list[[i]] <- data.frame(
      abundance_data[psi_cols],
      eta = p_detect_eta
    )

    # ----------------------------------------------------------
    # Collapse detection probability to site × OTU level
    # ----------------------------------------------------------
    site_p0 <- abundance_data |>
      dplyr::mutate(p0_abund = p0_abund) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(psi_cols))) |>
      dplyr::summarise(
        p0_site = prod(.data$p0_abund, na.rm = TRUE),
        .groups = "drop"
      )

    reduced_data <- reduced_data |>
      dplyr::select(-dplyr::any_of("p0_site")) |>
      dplyr::left_join(site_p0, by = psi_cols) |>
      dplyr::mutate(
        p0_site = tidyr::replace_na(.data$p0_site, 1)
      )

    # ----------------------------------------------------------
    # Update latent occupancy Z
    # ----------------------------------------------------------
    z_merge <- reduced_data |>
      dplyr::left_join(
        psi_long[, c(site_col, "OTU", "psi_prob")],
        by = psi_cols
      )

    zero_indices <- which(z_merge$z_obs == 0)

    if (length(zero_indices) > 0) {

      numerator <- z_merge$psi_prob[zero_indices] *
        z_merge$p0_site[zero_indices]

      denominator <- (1 - z_merge$psi_prob[zero_indices]) +
        z_merge$psi_prob[zero_indices] *
        z_merge$p0_site[zero_indices]

      posterior_z <- numerator / pmax(denominator, 1e-12)
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
  # Posterior summaries
  # ------------------------------------------------------------
  bind_summary_link <- function(lst, link_type = c("logit", "log", "cloglog")) {

    link_type <- match.arg(link_type)

    df <- dplyr::bind_rows(lst)

    if (nrow(df) == 0) {
      return(data.frame())
    }

    df_summary <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(psi_cols))) |>
      dplyr::summarise(
        eta_mean = mean(.data$eta, na.rm = TRUE),
        eta_var = stats::var(.data$eta, na.rm = TRUE),
        .groups = "drop"
      )

    df <- dplyr::left_join(df, df_summary, by = psi_cols) |>
      dplyr::mutate(
        weight = 1 / ifelse(is.na(.data$eta_var) | .data$eta_var == 0, 1e-6, .data$eta_var)
      )

    out <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(psi_cols))) |>
      dplyr::summarise(
        eta_mean = sum(.data$eta * .data$weight, na.rm = TRUE) /
          sum(.data$weight, na.rm = TRUE),
        eta_var = stats::var(.data$eta, na.rm = TRUE),
        se = sqrt(.data$eta_var),
        lwr_eta = .data$eta_mean - 1.96 * .data$se,
        upr_eta = .data$eta_mean + 1.96 * .data$se,
        .groups = "drop"
      )

    inv_link <- switch(
      link_type,
      logit = stats::plogis,
      log = exp,
      cloglog = function(x) 1 - exp(-exp(x))
    )

    out |>
      dplyr::mutate(
        mean = inv_link(.data$eta_mean),
        lwr = inv_link(.data$lwr_eta),
        upr = inv_link(.data$upr_eta)
      ) |>
      dplyr::select(dplyr::all_of(psi_cols), mean, se, lwr, upr)
  }

  keep <- seq.int(burn_in + 1, n_iter)

  psi_summary <- bind_summary_link(psi_list[keep], "logit")
  lambda_summary <- bind_summary_link(lambda_list[keep], "log")
  p_detect_summary <- bind_summary_link(p_detect_list[keep], "cloglog")

  final_summary <- psi_summary |>
    dplyr::rename_with(~ paste0("psi_", .x), -dplyr::all_of(psi_cols)) |>
    dplyr::left_join(
      lambda_summary |>
        dplyr::rename_with(~ paste0("lambda_", .x), -dplyr::all_of(psi_cols)),
      by = psi_cols
    ) |>
    dplyr::left_join(
      p_detect_summary |>
        dplyr::rename_with(~ paste0("p_detect_", .x), -dplyr::all_of(psi_cols)),
      by = psi_cols
    )

  # ------------------------------------------------------------
  # Latent variables
  # ------------------------------------------------------------
  lv_sites_combined <- dplyr::bind_rows(lv_sites_list[keep])
  lv_species_combined <- dplyr::bind_rows(lv_species_list[keep])

  mean_lv_sites <- lv_sites_combined |>
    dplyr::group_by(.data$Site) |>
    dplyr::summarise(
      dplyr::across(dplyr::starts_with("LV"), mean, na.rm = TRUE),
      .groups = "drop"
    )

  mean_lv_species <- lv_species_combined |>
    dplyr::group_by(.data$OTU) |>
    dplyr::summarise(
      dplyr::across(dplyr::starts_with("LV"), mean, na.rm = TRUE),
      .groups = "drop"
    )

  return(list(
    summary = final_summary,
    psi_list = psi_list[keep],
    lambda_list = lambda_list[keep],
    p_detect_list = p_detect_list[keep],
    occupancy_models = occupancy_models,
    abundance_models = abundance_models,
    reduced_data = reduced_data,
    lv_sites = lv_sites_combined,
    lv_species = lv_species_combined,
    mean_lv_sites = mean_lv_sites,
    mean_lv_species = mean_lv_species,
    filter_summary = list(
      otu_stats = otu_stats,
      kept_otus = keep_otus,
      min_species_sum = min_species_sum,
      abundance_threshold = abundance_threshold
    )
  ))
}
