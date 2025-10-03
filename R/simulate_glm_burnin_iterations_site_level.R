simulate_glm_burnin_iterations_site_level <- function(data_glm,
                                                      poisson_formula,
                                                      binomial_formula,
                                                      num_iterations = 100,
                                                      burn_in = 50,
                                                      site_var = "Site") {
  # Step 1: Aggregate data to one row per Site
  data_binom <- data_glm %>%
    group_by(.data[[site_var]]) %>%
    summarise(
      z_sim = as.integer(any(y > 0)),  # at least one detection per site
      across(-y, first),               # grab one value per covariate (assuming same across site)
      .groups = "drop"
    )

  poisson_models <- list()
  binomial_models <- list()

  for (iter in 1:num_iterations) {
    # Step 2: Fit binomial model on reduced (site-level) data
    model_binom <- glmmTMB::glmmTMB(binomial_formula,
                                    data = data_binom,
                                    family = binomial)

    # Step 3: Predict occupancy probabilities (P) for each site
    P_i <- predict(model_binom, type = "response", newdata = data_binom)

    # Step 4: Update z_sim (latent occupancy for each site)
    lambda_dummy <- rep(0, nrow(data_binom))  # λ not used yet here
    prob_z1_given_y0 <- P_i * exp(-lambda_dummy) / (P_i * exp(-lambda_dummy) + (1 - P_i))
    prob_z1_given_y0 <- pmin(pmax(prob_z1_given_y0, 1e-6), 1 - 1e-6)

    data_binom$z_sim[data_binom$z_sim == 0] <-
      rbinom(sum(data_binom$z_sim == 0), 1, prob_z1_given_y0[data_binom$z_sim == 0])

    # Step 5: Write site-level z_sim back to full data
    data_glm <- data_glm %>%
      left_join(data_binom %>% select(.data[[site_var]], z_sim), by = site_var)

    # Step 6: Subset full data where z_sim == 1
    Q <- data_glm[data_glm$z_sim == 1, ]

    # Step 7: Fit Poisson model on full data
    model_pois <- glmmTMB::glmmTMB(poisson_formula,
                                   family = poisson,
                                   data = Q)

    # Step 8: Save post-burn-in models
    if (iter > burn_in) {
      poisson_models[[iter - burn_in]] <- model_pois
      binomial_models[[iter - burn_in]] <- model_binom
    }
  }

  return(list(
    poisson_models = poisson_models,
    binomial_models = binomial_models,
    data_glm = data_glm,
    data_binom = data_binom
  ))
}
