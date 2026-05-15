#' Fit a Multispecies Occupancy--Detection Model using GLLVM and GLMM
#'
#' This function fits a hierarchical multispecies occupancy--abundance model
#' for eDNA data stored as a `phyloseq` object. Occupancy is modeled at the
#' site-by-OTU level using a binomial Generalized Linear Latent Variable Model
#' (GLLVM), while read abundance conditional on simulated occupancy is modeled
#' using a GLMM via `glmmTMB`.
#'
#' @param phyloseq A `phyloseq` object containing OTU count data and sample metadata.
#' @param site_col Character string giving the sample metadata column used as the site identifier.
#' @param abundance_rhs Right-hand side of the abundance model, for example `(1 | OTU)` or
#'   `offset(log(total_reads)) + (1 | OTU) + (1 | Name / OTU)`.
#' @param occupancy_covars Optional character vector of site-level covariates used in the GLLVM occupancy model.
#' @param min_species_sum Integer; minimum total read count required for an OTU to be retained. Default is 50.
#' @param min_detection_replicates Integer; minimum number of replicate rows with reads greater than
#'   `abundance_threshold` required for an OTU to be retained. Default is 1.
#' @param abundance_threshold Numeric; read-count threshold above which an OTU is treated as detected.
#'   Default is 1.
#' @param n_iter Integer; number of Monte Carlo EM-style iterations. Default is 50.
#' @param burn_in Integer; number of initial iterations discarded before summarising. Default is 10.
#' @param abundance_family Character string giving the abundance family. One of `"poisson"`,
#'   `"nbinom"`, `"zip"`, or `"zinb"`.
#' @param num_lv_c Integer; number of latent variables in the GLLVM occupancy model. Default is 2.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{summary}{Data frame with posterior summaries for occupancy probability (`psi_`),
#'   abundance intensity (`lambda_`), and site-level detection probability (`p_detect_`).}
#'   \item{psi_list}{List of occupancy linear predictors across post-burn-in iterations.}
#'   \item{lambda_list}{List of abundance linear predictors across post-burn-in iterations,
#'   aggregated to the site-by-OTU level.}
#'   \item{p_detect_list}{List of site-level detection linear predictors across post-burn-in iterations.}
#'   \item{occupancy_models}{List of fitted GLLVM occupancy models from all iterations.}
#'   \item{abundance_models}{List of fitted `glmmTMB` abundance models from all iterations.}
#'   \item{reduced_data}{Final site-by-OTU data frame containing observed and simulated occupancy states.}
#'   \item{lv_sites}{Site latent variable coordinates from post-burn-in GLLVM fits.}
#'   \item{lv_species}{Species/OTU latent variable coordinates from post-burn-in GLLVM fits.}
#'   \item{mean_lv_sites}{Post-burn-in mean site positions in latent variable space.}
#'   \item{mean_lv_species}{Post-burn-in mean OTU positions in latent variable space.}
#'   \item{filter_summary}{List containing OTU filtering information.}
#' }
#'
#' @details
#' The model separates the eDNA observation process into two components.
#' First, latent occupancy is defined at the site-by-OTU level. The observed
#' detection indicator is
#' \deqn{
#'   z_{s,m}^{obs} = I\left(\max_j y_{s,j,m} > c\right),
#' }
#' where \(y_{s,j,m}\) is the read count for OTU \(m\) in replicate/sample row
#' \(j\) at site \(s\), and \(c\) is `abundance_threshold`.
#'
#' Occupancy is modeled with a binomial GLLVM. The latent variables represent
#' residual site structure and OTU associations not explained by the supplied
#' occupancy covariates.
#'
#' Second, read abundance is modeled conditional on the current simulated
#' occupancy state, \(Z_{s,m}=1\), using `glmmTMB`. The abundance model may use
#' Poisson, negative binomial, zero-inflated Poisson, or zero-inflated negative
#' binomial families.
#'
#' At each iteration, the function:
#' \enumerate{
#'   \item fits a binomial GLLVM to the current simulated occupancy matrix;
#'   \item fits a GLMM abundance model to rows whose current simulated occupancy is one;
#'   \item predicts the probability of observing all zero reads for every site-by-OTU
#'   combination conditional on occupancy;
#'   \item updates unobserved occupancy states using
#'   \deqn{
#'   P(Z=1 \mid \text{all observed reads are zero}) =
#'   \frac{\psi p_0}{(1-\psi) + \psi p_0}.
#'   }
#' }
#'
#' Post-burn-in summaries are computed from the retained iterations. Intervals
#' are empirical 2.5% and 97.5% quantiles on the link scale and transformed to
#' the natural scale.
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#'
#' physeq <- subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Mock"))
#'
#' out <- FitModel_gllvm(
#'   phyloseq = physeq,
#'   site_col = "SampleType",
#'   abundance_rhs = (1 | OTU),
#'   occupancy_covars = c("SampleType"),
#'   min_species_sum = 50,
#'   min_detection_replicates = 1,
#'   abundance_threshold = 1,
#'   abundance_family = "poisson",
#'   n_iter = 10,
#'   burn_in = 2,
#'   num_lv_c = 2
#' )
#'
#' head(out$summary)
#' }
#'
#' @seealso \code{\link[gllvm]{gllvm}}, \code{\link[glmmTMB]{glmmTMB}}
#'
#' @import dplyr
#' @import tidyr
#' @importFrom gllvm gllvm
#' @importFrom glmmTMB glmmTMB nbinom2
#' @importFrom reshape2 melt acast
#' @importFrom utils head
#' @export
FitModel_gllvm <- function(
  phyloseq,
  site_col,
  abundance_rhs,
  occupancy_covars = NULL,
  min_species_sum = 50,
  min_detection_replicates = 1,
  abundance_threshold = 1,
  n_iter = 50,
  burn_in = 10,
  abundance_family = "poisson",
  num_lv_c = 2
) {

  if (burn_in >= n_iter) {
    stop("burn_in must be smaller than n_iter.")
  }

  valid_families <- c("poisson", "nbinom", "zip", "zinb")
  if (!(abundance_family %in% valid_families)) {
    stop(
      "'abundance_family' must be one of: ",
      paste(valid_families, collapse = ", ")
    )
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

  # ------------------------------------------------------------
  # Prepare long data
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
  # Explicit OTU filtering
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

  message("OTUs before filtering: ", dplyr::n_distinct(otu_stats$OTU))
  message("OTUs after filtering: ", dplyr::n_distinct(long_df$OTU))

  # ------------------------------------------------------------
  # Type handling
  # ------------------------------------------------------------
  factor_cols <- unique(c(site_col, "OTU", "Name", "Samplingmonth"))
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
  # Main EM-like iterations
  # ------------------------------------------------------------
  for (i in seq_len(n_iter)) {

    message("Iteration ", i)

    reduced_data[[site_col]] <- as.character(reduced_data[[site_col]])
    long_df[[site_col]] <- as.character(long_df[[site_col]])

    # ----------------------------------------------------------
    # Occupancy matrix: site × OTU
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
      stats::model.matrix(
        ~ .,
        data = cov_df[, occupancy_covars, drop = FALSE]
      )[, -1, drop = FALSE]
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
      num.lv = num_lv_c
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
    # Abundance data conditional on current Z = 1
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
    # ----------------------------------------------------------
    model_abundance <- glmmTMB::glmmTMB(
      formula = abundance_formula,
      data = abundance_data,
      family = abundance_glmm_family,
      ziformula = zi_formula
    )

    abundance_models[[i]] <- model_abundance

    # ----------------------------------------------------------
    # Lambda predictions on abundance_data
    # ----------------------------------------------------------
    lambda_pred <- stats::predict(
      model_abundance,
      type = "link",
      newdata = abundance_data,
      se.fit = TRUE,
      allow.new.levels = TRUE
    )

    lambda_eta <- as.numeric(lambda_pred$fit)

    lambda_list[[i]] <- data.frame(
      abundance_data[psi_cols],
      eta = lambda_eta
    ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(psi_cols))) |>
      dplyr::summarise(
        eta = mean(.data$eta, na.rm = TRUE),
        .groups = "drop"
      )

    # ----------------------------------------------------------
    # IMPORTANT:
    # Predict p0 for ALL long_df rows, not only z_sim == 1 rows
    # ----------------------------------------------------------
    pred_data_all <- long_df |>
      dplyr::left_join(
        reduced_data |>
          dplyr::select(dplyr::all_of(c(site_col, "OTU", "z_sim"))),
        by = c(site_col, "OTU")
      )

    pred_data_all <- to_factor_cols(pred_data_all, factor_cols)

    p0_all <- get_p0_abundance(
      fit = model_abundance,
      newdata = pred_data_all,
      family_name = abundance_family
    )

    site_p0 <- pred_data_all |>
      dplyr::mutate(p0_abund = p0_all) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(psi_cols))) |>
      dplyr::summarise(
        p0_site = prod(.data$p0_abund, na.rm = TRUE),
        .groups = "drop"
      )

    # Site-level detection probability:
    # p_detect = 1 - P(all PCR/sample reads are zero | Z = 1)
    p_detect_list[[i]] <- site_p0 |>
      dplyr::mutate(
        eta = log(-log(pmax(.data$p0_site, 1e-12)))
      ) |>
      dplyr::select(dplyr::all_of(psi_cols), eta)

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

    out <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(psi_cols))) |>
      dplyr::summarise(
        eta_mean = mean(.data$eta, na.rm = TRUE),
        eta_sd = stats::sd(.data$eta, na.rm = TRUE),
        lwr_eta = stats::quantile(.data$eta, 0.025, na.rm = TRUE),
        upr_eta = stats::quantile(.data$eta, 0.975, na.rm = TRUE),
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
        sd = .data$eta_sd,
        lwr = inv_link(.data$lwr_eta),
        upr = inv_link(.data$upr_eta)
      ) |>
      dplyr::select(dplyr::all_of(psi_cols), mean, sd, lwr, upr)
  }

  keep <- seq.int(burn_in + 1, n_iter)

  psi_summary <- bind_summary_link(psi_list[keep], "logit")
  lambda_summary <- bind_summary_link(lambda_list[keep], "log")
  p_detect_summary <- bind_summary_link(p_detect_list[keep], "cloglog")

  final_summary <- psi_summary |>
    dplyr::rename_with(
      ~ paste0("psi_", .x),
      -dplyr::all_of(psi_cols)
    ) |>
    dplyr::left_join(
      lambda_summary |>
        dplyr::rename_with(
          ~ paste0("lambda_", .x),
          -dplyr::all_of(psi_cols)
        ),
      by = psi_cols
    ) |>
    dplyr::left_join(
      p_detect_summary |>
        dplyr::rename_with(
          ~ paste0("p_detect_", .x),
          -dplyr::all_of(psi_cols)
        ),
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
      min_detection_replicates = min_detection_replicates,
      abundance_threshold = abundance_threshold
    )
  ))
}
                          
