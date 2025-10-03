simulate_glm_burnin_iterations <- function(data_glm,
                                           poisson_formula,
                                           binomial_formula,
                                           num_iterations = 100,
                                           burn_in = 50,
                                           species_var = "OTU") {
  # Step 1: Initialize latent occupancy variable
  data_glm$z_sim <- ifelse(data_glm$y > 0, 1, rbinom(nrow(data_glm), 1, 0.5))
  
  poisson_models <- list()
  binomial_models <- list()
  
  for (iter in 1:num_iterations) {
    # Subset data where occupancy = 1
    Q <- data_glm[data_glm$z_sim == 1, ]
    
    # Step 2: Fit Poisson and Binomial models
    model_poisson <- glmmTMB::glmmTMB(poisson_formula, family = poisson, data = Q)
    model_binomial <- glmmTMB::glmmTMB(binomial_formula, family = binomial, data = data_glm)
    
    # Step 3: Predict lambda (abundance) and P (occupancy)
    lambda_i <- predict(model_poisson, type = "response", newdata = data_glm)
    P_i <- predict(model_binomial, type = "response", newdata = data_glm)
    
    # Step 4: Update z_sim based on posterior probability
    prob_Z1_given_y0 <- P_i * exp(-lambda_i) / (P_i * exp(-lambda_i) + (1 - P_i))
    prob_Z1_given_y0 <- pmin(pmax(prob_Z1_given_y0, 1e-6), 1 - 1e-6)
    data_glm$z_sim[data_glm$y == 0] <-
      rbinom(sum(data_glm$y == 0), 1, prob_Z1_given_y0[data_glm$y == 0])
    
    # Step 5: Store models post burn-in
    if (iter > burn_in) {
      poisson_models[[iter - burn_in]] <- model_poisson
      binomial_models[[iter - burn_in]] <- model_binomial
    }
  }
  
  return(list(
    poisson_models = poisson_models,
    binomial_models = binomial_models
  ))
}
