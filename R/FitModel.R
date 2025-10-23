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
  # ──────────────────────────────────────
  # 1. Prepare data
  # ──────────────────────────────────────
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    min_species_sum = min_species_sum,
    sampletype_keep = sampletype_keep
  )
  long_df <- prep$long_df

  message("🔍 Calculating most abundant OTU from filtered phyloseq object...")
  otu_abundances <- taxa_sums(prep$physeq_filtered)
  top_otu <- names(sort(otu_abundances, decreasing = TRUE))[1]
  message("📌 Using most abundant OTU as reference: ", top_otu)

  long_df$OTU <- factor(long_df$OTU, levels = unique(long_df$OTU))
  long_df$OTU <- relevel(long_df$OTU, ref = top_otu)

  # Initialize storage lists
  psi_list <- list()
  lambda_list <- list()
  p_detect_list <- list()
  binomial_models <- list()
  poisson_models  <- list()

  # ──────────────────────────────────────
  # 2. Create occupancy (z_obs)
  # ──────────────────────────────────────
  reduced_data <- long_df %>%
    group_by(Site, OTU) %>%
    summarise(
      z_obs = as.integer(sum(y > 0) > abundance_threshold),
      across(-y, dplyr::first),
      .groups = "drop"
    )

  if (!is.null(treatment_exclude)) {
    reduced_data <- reduced_data %>% filter(treatment != treatment_exclude)
  }

  reduced_data <- reduced_data %>%
    mutate(
      treatment = droplevels(factor(treatment)),
      z_sim = z_obs
    )

  if (nlevels(reduced_data$treatment) < 2) {
    stop("❌ Error: 'treatment' must have at least 2 levels.")
  }

  # ──────────────────────────────────────
  # 3. Iterative Model Fitting
  # ──────────────────────────────────────
  for (i in 1:n_iter) {
    message("🔁 Iteration ", i)

    # --- Binomial model ---
    binomial_formula <- reformulate(termlabels = deparse(binomial_rhs), response = "z_sim")
    message("📌 Binomial model formula: ", deparse(binomial_formula))

    model_binomial <- glmmTMB::glmmTMB(
      formula = binomial_formula,
      data = reduced_data,
      family = binomial
    )

    pred_link <- predict(model_binomial, type = "link", se.fit = TRUE, newdata = reduced_data)
    logit_draw <- rnorm(n = nrow(reduced_data), mean = pred_link$fit, sd = pred_link$se.fit)
    psi_pred <- plogis(logit_draw)

    psi_list[[i]] <- data.frame(Site = reduced_data$Site, OTU = reduced_data$OTU,
                                treatment = reduced_data$treatment, eta = logit_draw)
    binomial_models[[i]] <- model_binomial

    # --- Poisson model ---
    full_data <- long_df
    if (!is.null(treatment_exclude)) {
      full_data <- full_data %>% filter(treatment != treatment_exclude)
    }

    full_data <- full_data %>%
      mutate(treatment = droplevels(factor(treatment))) %>%
      select(-any_of("z_sim")) %>%
      left_join(reduced_data[, c("Site", "OTU", "z_sim")], by = c("Site", "OTU"))

    poisson_data <- full_data %>% filter(z_sim == 1)

    poisson_formula <- reformulate(termlabels = deparse(poisson_rhs), response = "y")
    message("📌 Poisson model formula: ", deparse(poisson_formula))

    model_poisson <- glmmTMB::glmmTMB(
      formula = poisson_formula,
      data = poisson_data,
      family = poisson
    )

    lambda_pred <- predict(model_poisson, type = "link", se.fit = TRUE, newdata = poisson_data)
    lambda <- exp(lambda_pred$fit)

    lambda_list[[i]] <- data.frame(Site = poisson_data$Site, OTU = poisson_data$OTU,
                                   treatment = poisson_data$treatment, eta = lambda_pred$fit)

    p_detect_list[[i]] <- data.frame(Site = poisson_data$Site, OTU = poisson_data$OTU,
                                     treatment = poisson_data$treatment, eta = lambda_pred$fit)

    poisson_models[[i]] <- model_poisson

    # Update z_sim
    lambda_total <- data.frame(Site = poisson_data$Site, OTU = poisson_data$OTU, lambda_pred = lambda) %>%
      group_by(Site, OTU) %>%
      summarise(lambda_prod = 1 - prod(1 - exp(-lambda_pred)), .groups = "drop")

    reduced_data <- reduced_data %>%
      select(-any_of("lambda_prod")) %>%
      left_join(lambda_total, by = c("Site", "OTU")) %>%
      mutate(lambda_prod = tidyr::replace_na(lambda_prod, 0))

    zero_indices <- which(reduced_data$z_obs == 0)
    adjusted_prob <- psi_pred[zero_indices] * reduced_data$lambda_prod[zero_indices]
    adjusted_prob <- pmin(pmax(adjusted_prob, 0.001), 0.999)
    reduced_data$z_sim[zero_indices] <- rbinom(length(zero_indices), 1, adjusted_prob)
  }

  # ──────────────────────────────────────
  # 4. Posterior Summaries
  # ──────────────────────────────────────
  bind_summary_link <- function(lst, link_type = c("logit", "log", "cloglog")) {
    link_type <- match.arg(link_type)
    df <- bind_rows(lst)

    df_summary <- df %>%
      group_by(Site, OTU, treatment) %>%
      summarise(eta_mean = mean(eta, na.rm = TRUE), eta_var = var(eta, na.rm = TRUE), .groups = "drop")

    df <- df %>%
      left_join(df_summary, by = c("Site", "OTU", "treatment")) %>%
      mutate(weight = 1 / ifelse(is.na(eta_var) | eta_var == 0, 1e-6, eta_var))

    weighted_summary <- df %>%
      group_by(Site, OTU, treatment) %>%
      summarise(
        eta_mean = sum(eta * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE),
        eta_var = sum(weight * (eta - eta_mean)^2, na.rm = TRUE) /
          ((sum(weight)^2 - sum(weight^2)) / sum(weight)),
        se_eta = sqrt(eta_var),
        .groups = "drop"
      ) %>%
      mutate(
        lwr_eta = eta_mean - 1.96 * se_eta,
        upr_eta = eta_mean + 1.96 * se_eta,
        mean = case_when(
          link_type == "logit" ~ plogis(eta_mean),
          link_type == "log" ~ exp(eta_mean),
          link_type == "cloglog" ~ 1 - exp(-exp(eta_mean))
        ),
        lwr = case_when(
          link_type == "logit" ~ plogis(lwr_eta),
          link_type == "log" ~ exp(lwr_eta),
          link_type == "cloglog" ~ 1 - exp(-exp(upr_eta))
        ),
        upr = case_when(
          link_type == "logit" ~ plogis(upr_eta),
          link_type == "log" ~ exp(upr_eta),
          link_type == "cloglog" ~ 1 - exp(-exp(lwr_eta))
        ),
        se = se_eta
      ) %>%
      select(Site, OTU, treatment, mean, se, lwr, upr)

    return(weighted_summary)
  }

  # Apply burn-in
  psi_list_burned      <- psi_list[-seq_len(burn_in)]
  lambda_list_burned   <- lambda_list[-seq_len(burn_in)]
  p_detect_list_burned <- p_detect_list[-seq_len(burn_in)]

  # Final summaries
  psi_summary      <- bind_summary_link(psi_list_burned,      link_type = "logit")
  lambda_summary   <- bind_summary_link(lambda_list_burned,   link_type = "log")
  p_detect_summary <- bind_summary_link(p_detect_list_burned, link_type = "cloglog")

  final_summary <- psi_summary %>%
    rename_with(~paste0("psi_", .), -c(Site, OTU, treatment)) %>%
    left_join(rename_with(lambda_summary, ~paste0("lambda_", .), -c(Site, OTU, treatment)),
              by = c("Site", "OTU", "treatment")) %>%
    left_join(rename_with(p_detect_summary, ~paste0("p_detect_", .), -c(Site, OTU, treatment)),
              by = c("Site", "OTU", "treatment"))

  return(list(
    summary = final_summary,
    psi_list = psi_list_burned,
    lambda_list = lambda_list_burned,
    p_detect_list = p_detect_list_burned,
    binomial_models = binomial_models,
    poisson_models = poisson_models,
    reduced_data = reduced_data
  ))
}
