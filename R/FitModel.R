#' Fit Site-Level Occupancy Model with Detection and Abundance Components
#'
#' @param phyloseq A `phyloseq` object with OTU table and sample data.
#' @param poisson_rhs Right-hand side formula (quoted) for the Poisson (abundance) model.
#' @param binomial_rhs Right-hand side formula (quoted) for the binomial (occupancy) model.
#' @param min_species_sum Minimum OTU abundance threshold to retain OTUs.
#' @param sampletype_keep Sample type to filter for (e.g., "biologicalsample").
#' @param abundance_threshold Integer threshold for detecting presence/absence.
#' @param treatment_exclude Optional character to exclude a treatment level (e.g., "control").
#' @param n_iter Number of MCMC iterations.
#' @param burn_in Number of burn-in iterations to discard.
#'
#' @return A list with summaries and iteration-level posterior samples.
#' @export
FitModel <- function(phyloseq,
                     poisson_rhs,
                     binomial_rhs,
                     min_species_sum = 50,
                     sampletype_keep = NULL,
                     abundance_threshold = 1,
                     treatment_exclude = NULL,
                     n_iter = 50,
                     burn_in = 10) {

  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    min_species_sum = min_species_sum,
    sampletype_keep = sampletype_keep
  )
  long_df <- prep$long_df

  message(" Calculating most abundant OTU...")
  otu_abundances <- taxa_sums(prep$physeq_filtered)
  top_otu <- names(sort(otu_abundances, decreasing = TRUE))[1]
  long_df$OTU <- factor(long_df$OTU, levels = unique(long_df$OTU))
  long_df$OTU <- relevel(long_df$OTU, ref = top_otu)

  psi_list <- list()
  lambda_list <- list()
  p_detect_list <- list()
  binomial_models <- list()
  poisson_models  <- list()

  reduced_data <- long_df %>%
    group_by(Site, OTU) %>%
    summarise(z_obs = as.integer(sum(y > 0) > abundance_threshold),
              across(-y, dplyr::first), .groups = "drop")

  if (!is.null(treatment_exclude)) {
    reduced_data <- reduced_data %>% filter(treatment != treatment_exclude)
  }

  reduced_data <- reduced_data %>%
    mutate(treatment = droplevels(factor(treatment)),
           z_sim = z_obs)

  if (nlevels(reduced_data$treatment) < 2) {
    stop("Error: 'treatment' must have at least 2 levels.")
  }

  for (i in 1:n_iter) {
    message(" Iteration ", i)

    model_binomial <- glmmTMB::glmmTMB(
      formula = reformulate(deparse(binomial_rhs), response = "z_sim"),
      data = reduced_data,
      family = binomial
    )

    pred_link <- predict(model_binomial, type = "link", se.fit = TRUE, newdata = reduced_data)
    logit_draw <- rnorm(nrow(reduced_data), mean = pred_link$fit, sd = pred_link$se.fit)
    psi_pred <- plogis(logit_draw)

    psi_list[[i]] <- data.frame(Site = reduced_data$Site, OTU = reduced_data$OTU,
                                treatment = reduced_data$treatment, eta = logit_draw)
    binomial_models[[i]] <- model_binomial

    poisson_data <- long_df %>%
      filter(is.null(treatment_exclude) | treatment != treatment_exclude) %>%
      mutate(treatment = droplevels(factor(treatment))) %>%
      left_join(reduced_data[, c("Site", "OTU", "z_sim")], by = c("Site", "OTU")) %>%
      filter(z_sim == 1)

    model_poisson <- glmmTMB::glmmTMB(
      formula = reformulate(deparse(poisson_rhs), response = "y"),
      data = poisson_data,
      family = poisson
    )

    lambda_pred <- predict(model_poisson, type = "link", se.fit = TRUE, newdata = poisson_data)
    lambda <- exp(lambda_pred$fit)

    lambda_list[[i]] <- data.frame(Site = poisson_data$Site, OTU = poisson_data$OTU,
                                   treatment = poisson_data$treatment, eta = lambda_pred$fit)

    p_detect_list[[i]] <- data.frame(
      Site = poisson_data$Site,
      OTU = poisson_data$OTU,
      treatment = poisson_data$treatment,
      eta = log(-log(1 - (1 - exp(-lambda))))  # cloglog link
    )

    poisson_models[[i]] <- model_poisson

    lambda_total <- poisson_data %>%
      mutate(lambda_pred = lambda) %>%
      group_by(Site, OTU) %>%
      summarise(lambda_prod = 1 - prod(1 - exp(-lambda_pred)), .groups = "drop")

    reduced_data <- reduced_data %>%
      select(-any_of("lambda_prod")) %>%
      left_join(lambda_total, by = c("Site", "OTU")) %>%
      mutate(lambda_prod = tidyr::replace_na(lambda_prod, 0))

    zero_indices <- which(reduced_data$z_obs == 0)
    adjusted_prob <- psi_pred[zero_indices] * reduced_data$lambda_prod[zero_indices]
    reduced_data$z_sim[zero_indices] <- rbinom(length(zero_indices), 1,
                                               pmin(pmax(adjusted_prob, 0.001), 0.999))
  }

  bind_summary_link <- function(lst, link_type = c("logit", "log", "cloglog")) {
    link_type <- match.arg(link_type)
    df <- bind_rows(lst)

    df_summary <- df %>%
      group_by(Site, OTU, treatment) %>%
      summarise(eta_mean = mean(eta), eta_var = var(eta), .groups = "drop")

    df <- df %>%
      left_join(df_summary, by = c("Site", "OTU", "treatment")) %>%
      mutate(weight = 1 / ifelse(is.na(eta_var) | eta_var == 0, 1e-6, eta_var))

    df %>%
      group_by(Site, OTU, treatment) %>%
      summarise(
        eta_mean = sum(eta * weight) / sum(weight),
        eta_var = sum(weight * (eta - eta_mean)^2) /
          ((sum(weight)^2 - sum(weight^2)) / sum(weight)),
        se_eta = sqrt(eta_var),
        lwr_eta = eta_mean - 1.96 * se_eta,
        upr_eta = eta_mean + 1.96 * se_eta,
        mean = switch(link_type,
                      logit = plogis(eta_mean),
                      log = exp(eta_mean),
                      cloglog = 1 - exp(-exp(eta_mean))),
        lwr = switch(link_type,
                     logit = plogis(lwr_eta),
                     log = exp(lwr_eta),
                     cloglog = 1 - exp(-exp(upr_eta))),
        upr = switch(link_type,
                     logit = plogis(upr_eta),
                     log = exp(upr_eta),
                     cloglog = 1 - exp(-exp(lwr_eta))),
        se = se_eta,
        .groups = "drop"
      ) %>%
      select(Site, OTU, treatment, mean, se, lwr, upr)
  }

  psi_summary      <- bind_summary_link(psi_list[-seq_len(burn_in)], "logit")
  lambda_summary   <- bind_summary_link(lambda_list[-seq_len(burn_in)], "log")
  p_detect_summary <- bind_summary_link(p_detect_list[-seq_len(burn_in)], "cloglog")

  final_summary <- psi_summary %>%
    rename_with(~paste0("psi_", .), -c(Site, OTU, treatment)) %>%
    left_join(rename_with(lambda_summary, ~paste0("lambda_", .), -c(Site, OTU, treatment)),
              by = c("Site", "OTU", "treatment")) %>%
    left_join(rename_with(p_detect_summary, ~paste0("p_detect_", .), -c(Site, OTU, treatment)),
              by = c("Site", "OTU", "treatment"))

  return(list(
    summary = final_summary,
    psi_list = psi_list,
    lambda_list = lambda_list,
    p_detect_list = p_detect_list,
    binomial_models = binomial_models,
    poisson_models = poisson_models,
    reduced_data = reduced_data
  ))
}
