#' Fit Hierarchical Occupancy-Detection Model using GLMMs
#'
#' Automatically detects treatment, sample, and replicate columns from model formulae.
#' Uses a Gibbs-style loop for estimating detection (lambda), occupancy (psi), and 
#' detection probability (p_detect).
#'
#' @param phyloseq A phyloseq object with OTU table and sample metadata.
#' @param site_col Name of the column representing site-level grouping.
#' @param poisson_rhs Right-hand-side of Poisson model (as a formula inside quote()).
#' @param binomial_rhs Right-hand-side of Binomial model (as a formula inside quote()).
#' @param min_species_sum Minimum total abundance required to retain a taxon.
#' @param abundance_threshold Minimum detection threshold to consider presence at site level.
#' @param n_iter Number of Gibbs iterations.
#' @param burn_in Number of iterations to discard as burn-in.
#'
#' @return A list with `summary`, `psi_list`, `lambda_list`, `p_detect_list`, models, and `reduced_data`.
#'
#' @importFrom dplyr group_by summarise mutate select filter left_join arrange ungroup across
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats predict rnorm reformulate plogis var
#' @export
FitModel <- function(phyloseq,
                     site_col,
                     poisson_rhs,
                     binomial_rhs,
                     min_species_sum = 50,
                     abundance_threshold = 1,
                     n_iter = 50,
                     burn_in = 10) {

  # --- Detect column names from model terms ---
  rhs_text <- paste(deparse(poisson_rhs), collapse = " ")
  rhs_text <- gsub("`", "", rhs_text)  # Remove backticks
  rhs_terms <- unique(all.vars(poisson_rhs))

  treatment_candidates <- rhs_terms[sapply(rhs_terms, function(x)
    grepl(paste0("\\b", x, "\\b\\s*\\*\\s*OTU"), rhs_text) |
    grepl(paste0("OTU\\s*\\*\\s*\\b", x, "\\b"), rhs_text))]

  nested_candidates <- rhs_terms[sapply(rhs_terms, function(x)
    grepl(paste0("\\b", x, "\\b\\s*/\\s*OTU"), rhs_text) |
    grepl(paste0("OTU\\s*/\\s*\\b", x, "\\b"), rhs_text))]

  #if (length(treatment_candidates) == 0)
   # stop(" No treatment variable (e.g., Samplingmonth * OTU) found.")
  #if (length(nested_candidates) < 2)
    #stop(" Need at least two nested variables (e.g., Sample / OTU, Replicate / OTU).")

  treatment_col <- treatment_candidates[1]
  sample_col    <- nested_candidates[1]
  replicate_col <- nested_candidates[2]

  message(" Detected columns:",
          "\n   Treatment: ", treatment_col,
          "\n   Sample: ", sample_col,
          "\n   Replicate: ", replicate_col)

  # --- Prepare data ---
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    min_species_sum = min_species_sum,
    site_col = site_col
  )

  long_df <- prep$long_df

  # --- Setup OTU factors ---
  message(" Calculating most abundant OTU...")
  otu_abundances <- taxa_sums(prep$physeq_filtered)
  top_otu <- names(sort(otu_abundances, decreasing = TRUE))[1]
  long_df$OTU <- factor(long_df$OTU, levels = unique(long_df$OTU))
  long_df$OTU <- relevel(long_df$OTU, ref = top_otu)

  # --- Prepare storage ---
  psi_list <- list()
  lambda_list <- list()
  p_detect_list <- list()
  binomial_models <- list()
  poisson_models <- list()

  treatment_levels <- levels(factor(long_df[[treatment_col]]))
  message(" Treatment levels: ", paste(treatment_levels, collapse = ", "))

  # --- Site-level reduction ---
  reduced_data <- long_df %>%
    group_by(.data[[site_col]], OTU) %>%
    summarise(z_obs = as.integer(sum(y > 0) > abundance_threshold),
              across(-y, dplyr::first), .groups = "drop") %>%
    mutate(
      !!treatment_col := factor(.data[[treatment_col]], levels = treatment_levels),
      z_sim = z_obs
    )

  if (nlevels(reduced_data[[treatment_col]]) < 2) {
    stop(" treatment_col must have ≥ 2 levels.")
   }

  # --- Iterative estimation ---
  for (i in 1:n_iter) {
    message(" Iteration ", i)

    model_binomial <- glmmTMB::glmmTMB(
      formula = reformulate(deparse(binomial_rhs), response = "z_sim"),
      data = reduced_data,
      family = binomial
    )

    pred_link <- predict(model_binomial, type = "link", se.fit = TRUE)
    logit_draw <- rnorm(nrow(reduced_data), mean = pred_link$fit, sd = pred_link$se.fit)
    psi_pred <- plogis(logit_draw)

    psi_list[[i]] <- setNames(
      data.frame(reduced_data[[site_col]], reduced_data$OTU, reduced_data[[treatment_col]], eta = logit_draw),
      c(site_col, "OTU", treatment_col, "eta")
    )
    binomial_models[[i]] <- model_binomial

    poisson_data <- long_df %>%
      mutate(!!treatment_col := droplevels(factor(.data[[treatment_col]]))) %>%
      left_join(reduced_data %>% select(all_of(c(site_col, "OTU", "z_sim"))), by = c(site_col, "OTU")) %>%
      filter(z_sim == 1)

    model_poisson <- glmmTMB::glmmTMB(
      formula = reformulate(deparse(poisson_rhs), response = "y"),
      data = poisson_data,
      family = poisson
    )

    lambda_pred <- predict(model_poisson, type = "link", se.fit = TRUE)
    lambda <- exp(lambda_pred$fit)

    lambda_list[[i]] <- setNames(
      data.frame(poisson_data[[site_col]], poisson_data$OTU, poisson_data[[treatment_col]], eta = lambda_pred$fit),
      c(site_col, "OTU", treatment_col, "eta")
    )

    p_detect_list[[i]] <- setNames(
      data.frame(poisson_data[[site_col]], poisson_data$OTU, poisson_data[[treatment_col]], eta = 1 - exp(-lambda)),
      c(site_col, "OTU", treatment_col, "eta")
    )

    poisson_models[[i]] <- model_poisson

    lambda_total <- poisson_data %>%
      mutate(lambda_pred = lambda) %>%
      group_by(.data[[site_col]], OTU) %>%
      summarise(lambda_prod = 1 - prod(1 - exp(-lambda_pred)), .groups = "drop")

    reduced_data <- reduced_data %>%
      select(-any_of("lambda_prod")) %>%
      left_join(lambda_total, by = c(site_col, "OTU")) %>%
      mutate(lambda_prod = tidyr::replace_na(lambda_prod, 0))

    zero_indices <- which(reduced_data$z_obs == 0)
    adjusted_prob <- psi_pred[zero_indices] * reduced_data$lambda_prod[zero_indices]
    reduced_data$z_sim[zero_indices] <- rbinom(length(zero_indices), 1,
                                               pmin(pmax(adjusted_prob, 0.001), 0.999))
  }

  # --- Summary wrapper ---
  bind_summary_link <- function(lst, link_type = c("logit", "log", "cloglog")) {
    link_type <- match.arg(link_type)
    df <- bind_rows(lst)

    df_summary <- df %>%
      group_by(.data[[site_col]], OTU, .data[[treatment_col]]) %>%
      summarise(eta_mean = mean(eta), eta_var = var(eta), .groups = "drop")

    df <- df %>%
      left_join(df_summary, by = c(site_col, "OTU", treatment_col)) %>%
      mutate(weight = 1 / ifelse(is.na(eta_var) | eta_var == 0, 1e-6, eta_var))

    df %>%
      group_by(.data[[site_col]], OTU, .data[[treatment_col]]) %>%
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
      select(all_of(c(site_col, "OTU", treatment_col, "mean", "se", "lwr", "upr")))
  }

  # --- Final summary ---
  psi_summary      <- bind_summary_link(psi_list[-seq_len(burn_in)], "logit")
  lambda_summary   <- bind_summary_link(lambda_list[-seq_len(burn_in)], "log")
  p_detect_summary <- bind_summary_link(p_detect_list[-seq_len(burn_in)], "cloglog")

  final_summary <- psi_summary %>%
    rename_with(~paste0("psi_", .), -all_of(c(site_col, "OTU", treatment_col))) %>%
    left_join(rename_with(lambda_summary, ~paste0("lambda_", .), -all_of(c(site_col, "OTU", treatment_col))),
              by = c(site_col, "OTU", treatment_col)) %>%
    left_join(rename_with(p_detect_summary, ~paste0("p_detect_", .), -all_of(c(site_col, "OTU", treatment_col))),
              by = c(site_col, "OTU", treatment_col))

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
