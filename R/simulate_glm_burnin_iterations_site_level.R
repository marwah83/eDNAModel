simulate_glm_burnin_iterations_site_level <- function(data_glm,
                                                      poisson_formula,
                                                      binomial_formula,
                                                      num_iterations = 100,
                                                      burn_in = 50,
                                                      site_var = "Site") {
  # STEP 0: Initialize z_sim in full data
  data_glm$z_sim <- ifelse(data_glm$y > 0, 1, rbinom(nrow(data_glm), 1, 0.5))
  
  # Storage for models
  poisson_models <- list()
  binomial_models <- list()
  
  for (iter in 1:num_iterations) {
    
    # STEP 1: Aggregate to site-level
    data_binom <- data_glm %>%
      dplyr::group_by(.data[[site_var]]) %>%
      dplyr::summarise(
        z_sim = as.integer(any(y > 0)),
        across(-y, dplyr::first),
        .groups = "drop"
      )
    
    # STEP 2: Fit binomial model on site-level data
    model_binomial <- glmmTMB::glmmTMB(binomial_formula, data = data_binom, family = binomial)
    
    # STEP 3: Predict site-level occupancy probabilities
    psi_pred <- predict(model_binomial, type = "response", newdata = data_binom)
    
    # STEP 4: Simulate new z_sim for site-level data
    data_binom$z_sim <- rbinom(nrow(data_binom), 1, psi_pred)
    
    # STEP 5: Write z_sim back to full data
    data_glm <- data_glm %>%
      dplyr::left_join(data_binom[, c(site_var, "z_sim")], by = site_var, suffix = c("", ".new")) %>%
      dplyr::mutate(z_sim = .data$z_sim.new) %>%
      dplyr::select(-z_sim.new)
    
    # STEP 6: Subset full data where site was occupied
    Q <- data_glm[data_glm$z_sim == 1, ]
    
    # STEP 7: Fit Poisson model on full (filtered) data
    model_poisson <- glmmTMB::glmmTMB(poisson_formula, family = poisson, data = Q)
    
    # STEP 8: Store models after burn-in
    if (iter > burn_in) {
      binomial_models[[iter - burn_in]] <- model_binomial
      poisson_models[[iter - burn_in]] <- model_poisson
    }
  }
  
  return(list(
    poisson_models = poisson_models,
    binomial_models = binomial_models
  ))
}
