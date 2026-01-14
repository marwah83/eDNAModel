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
#'
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#' result <- FitModel_gllvm(GlobalPatterns,
#'                          site_col = "SampleID",
#'                          abundance_rhs = ~ SampleType,
#'                          occupancy_covars = c("SampleType", "X.SampleID"))
#' }
#'
#' @seealso \code{\link[gllvm]{gllvm}}, \code{\link[glmmTMB]{glmmTMB}}
#' @import dplyr
#' @import tidyr
#' @importFrom stats reformulate predict rbinom
#' @importFrom gllvm gllvm
#' @importFrom glmmTMB glmmTMB
#' @importFrom reshape2 melt acast
#' @importFrom utils head
#' @export
FitModel_gllvm <- function(phyloseq,
                     site_col,
                     abundance_rhs,
                     occupancy_covars = NULL,
                     min_species_sum = 50,
                     abundance_threshold = 1,
                     n_iter = 50,
                     burn_in = 10,
                     abundance_family = "poisson") {
  
  
  abundance_rhs <- substitute(abundance_rhs)
  
  # Helper: Convert columns to factor
  to_factor_cols <- function(df, cols) {
    for (col in cols) {
      if (col %in% colnames(df)) {
        df[[col]] <- as.factor(as.character(df[[col]]))
      }
    }
    return(df)
  }
  
  # Check abundance family
  valid_families <- c("poisson", "nbinom", "zip", "zinb")
  if (!abundance_family %in% valid_families) {
    stop(" 'abundance_family' must be one of: ", paste(valid_families, collapse = ", "))
  }
  
  # Prepare data
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    min_species_sum = min_species_sum,
    site_col = site_col
  )
  long_df <- prep$long_df
  long_df <- to_factor_cols(long_df, c(site_col, "OTU", "Name", "Samplingmonth"))
  
  # Set reference OTU
  otu_abundances <- taxa_sums(prep$physeq_filtered)
  top_otu <- names(sort(otu_abundances, decreasing = TRUE))[1]
  long_df$OTU <- relevel(factor(long_df$OTU, levels = unique(long_df$OTU)), ref = top_otu)
  
  # z_obs
  reduced_data <- long_df %>%
    group_by(across(all_of(c(site_col, "OTU")))) %>%
    summarise(z_obs = as.integer(sum(y > 0) > abundance_threshold),
              across(-y, first), .groups = "drop") %>%
    mutate(z_sim = z_obs)
  reduced_data <- to_factor_cols(reduced_data, c(site_col, "OTU", "Name", "Samplingmonth"))
  
  # Initialize
  psi_list <- list()
  lambda_list <- list()
  p_detect_list <- list()
  occupancy_models <- list()
  abundance_models <- list()
  lv_sites_list <- list()
  lv_species_list <- list()
  
  for (i in 1:n_iter) {
    message(" Iteration ", i)
    
    reduced_data[[site_col]] <- as.character(reduced_data[[site_col]])
    long_df[[site_col]] <- as.character(long_df[[site_col]])
    
    # z matrix
    z_matrix <- acast(reduced_data, formula = paste(site_col, "~ OTU"), value.var = "z_sim", fill = 0)
    z_sites <- rownames(z_matrix)
    
    # Covariates
    cov_df <- reduced_data %>%
      select(all_of(c(site_col, occupancy_covars))) %>%
      distinct()
    cov_df <- to_factor_cols(cov_df, occupancy_covars)
    cov_df <- cov_df[match(z_sites, cov_df[[site_col]]), ]
    
    X_cov <- if (!is.null(occupancy_covars) && length(occupancy_covars) > 0) {
      model.matrix(~ ., data = cov_df)[, -1, drop = FALSE]
    } else NULL
    
    # === GLLVM ===
    model_occupancy <- gllvm::gllvm(
      y = z_matrix,
      X = X_cov,
      family = "binomial",
      num.lv.c = 2
    )
    
    # Latent variables
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
    
    # Predict ψ
    psi_pred <- predict(model_occupancy, type = "response")
    rownames(psi_pred) <- rownames(z_matrix)
    psi_long <- melt(psi_pred, varnames = c(site_col, "OTU"), value.name = "eta")
    psi_list[[i]] <- psi_long
    occupancy_models[[i]] <- model_occupancy
    
    # Abundance model
    abundance_data <- long_df %>%
      left_join(reduced_data %>% select(all_of(c(site_col, "OTU", "z_sim"))),
                by = c(site_col, "OTU")) %>%
      filter(z_sim == 1)
    abundance_data <- to_factor_cols(abundance_data, c(site_col, "OTU", "Name", "Samplingmonth"))
    
    # GLMM settings
    if (abundance_family == "poisson") {
      abundance_glmm_family <- poisson
      zi_formula <- ~0
    } else if (abundance_family == "nbinom") {
      abundance_glmm_family <- nbinom2
      zi_formula <- ~0
    } else if (abundance_family == "zip") {
      abundance_glmm_family <- poisson
      zi_formula <- ~1
    } else {
      abundance_glmm_family <- nbinom2
      zi_formula <- ~1
    }
    
    model_abundance <- glmmTMB::glmmTMB(
      formula = reformulate(deparse(abundance_rhs), response = "y"),
      data = abundance_data,
      family = abundance_glmm_family,
      ziformula = zi_formula
    )
    
    lambda_pred <- predict(model_abundance, type = "link", se.fit = TRUE)
    lambda <- exp(lambda_pred$fit)
    
    psi_cols <- c(site_col, "OTU")
    lambda_list[[i]] <- data.frame(abundance_data[psi_cols], eta = lambda_pred$fit)
    p_detect_list[[i]] <- data.frame(abundance_data[psi_cols], eta = 1 - exp(-lambda))
    abundance_models[[i]] <- model_abundance
    
    lambda_total <- abundance_data %>%
      mutate(lambda_pred = lambda) %>%
      group_by(across(all_of(c(site_col, "OTU")))) %>%
      summarise(lambda_prod = 1 - prod(1 - exp(-lambda_pred)), .groups = "drop")
    
    reduced_data <- reduced_data %>%
      select(-any_of("lambda_prod")) %>%
      left_join(lambda_total, by = c(site_col, "OTU")) %>%
      mutate(lambda_prod = replace_na(lambda_prod, 0))
    
    z_merge <- merge(reduced_data, psi_long, by = c(site_col, "OTU"))
    zero_indices <- which(z_merge$z_obs == 0)
    adjusted_prob <- z_merge$eta[zero_indices] * z_merge$lambda_prod[zero_indices]
    reduced_data$z_sim[zero_indices] <- rbinom(length(zero_indices), 1,
                                               pmin(pmax(adjusted_prob, 0.001), 0.999))
  }
  
  # Posterior summary
  bind_summary_link <- function(lst, link_type = c("logit", "log", "cloglog")) {
    link_type <- match.arg(link_type)
    df <- bind_rows(lst)
    df_summary <- df %>%
      group_by(across(all_of(psi_cols))) %>%
      summarise(
        eta_mean = mean(eta),
        eta_var = var(eta),
        .groups = "drop"
      )
    df <- left_join(df, df_summary, by = psi_cols) %>%
      mutate(weight = 1 / ifelse(is.na(eta_var) | eta_var == 0, 1e-6, eta_var))
    df %>%
      group_by(across(all_of(psi_cols))) %>%
      summarise(
        eta_mean = sum(eta * weight) / sum(weight),
        eta_var = var(eta),
        se = sqrt(eta_var),
        lwr = eta_mean - 1.96 * se,
        upr = eta_mean + 1.96 * se,
        mean = switch(link_type,
                      logit = plogis(eta_mean),
                      log = exp(eta_mean),
                      cloglog = 1 - exp(-exp(eta_mean))),
        lwr = switch(link_type,
                     logit = plogis(lwr),
                     log = exp(lwr),
                     cloglog = 1 - exp(-exp(upr))),
        upr = switch(link_type,
                     logit = plogis(upr),
                     log = exp(upr),
                     cloglog = 1 - exp(-exp(lwr))),
        .groups = "drop"
      ) %>%
      select(all_of(psi_cols), mean, se, lwr, upr)
  }
  
  psi_summary <- bind_summary_link(psi_list[-seq_len(burn_in)], "logit")
  lambda_summary <- bind_summary_link(lambda_list[-seq_len(burn_in)], "log")
  p_detect_summary <- bind_summary_link(p_detect_list[-seq_len(burn_in)], "cloglog")
  
  final_summary <- psi_summary %>%
    rename_with(~paste0("psi_", .), -all_of(psi_cols)) %>%
    left_join(rename_with(lambda_summary, ~paste0("lambda_", .), -all_of(psi_cols)), by = psi_cols) %>%
    left_join(rename_with(p_detect_summary, ~paste0("p_detect_", .), -all_of(psi_cols)), by = psi_cols)
  
  # Biplot latent variables
  lv_sites_combined <- bind_rows(lv_sites_list[-seq_len(burn_in)])
  lv_species_combined <- bind_rows(lv_species_list[-seq_len(burn_in)])
  
  mean_lv_sites <- lv_sites_combined %>%
    group_by(Site) %>%
    summarise(across(starts_with("LV"), mean), .groups = "drop")
  
  mean_lv_species <- lv_species_combined %>%
    group_by(OTU) %>%
    summarise(across(starts_with("LV"), mean), .groups = "drop")
  
  return(list(
    summary = final_summary,
    psi_list = psi_list[-seq_len(burn_in)],
    lambda_list = lambda_list[-seq_len(burn_in)],
    p_detect_list = p_detect_list[-seq_len(burn_in)],
    occupancy_models = occupancy_models,
    abundance_models = abundance_models,
    reduced_data = reduced_data,
    lv_sites = lv_sites_combined,
    lv_species = lv_species_combined,
    mean_lv_sites = mean_lv_sites,
    mean_lv_species = mean_lv_species
  ))
}

