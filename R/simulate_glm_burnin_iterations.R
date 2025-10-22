#' Simulate Latent Occupancy-Abundance Models via GLMMs
#'
#' Performs iterative fitting of Poisson (abundance) and Binomial (occupancy) 
#' generalized linear mixed models (GLMMs) using the `glmmTMB` package. 
#' The method introduces a latent binary variable (`z_sim`) to separate 
#' detection from true presence and iteratively refines model estimates 
#' across multiple iterations, with burn-in.
#'
#' @param data_glm A long-format data frame prepared from a `phyloseq` object. Must include species, count (`y`), and covariate columns.
#' @param poisson_formula A `formula` specifying the abundance (Poisson) model.
#' @param binomial_formula A `formula` specifying the occupancy (Binomial) model.
#' @param num_iterations Total number of iterations to run (default: 100).
#' @param burn_in Number of burn-in iterations to discard from the results (default: 50).
#' @param species_var Name of the column representing the species (default: "OTU").
#'
#' @return A list containing:
#' \describe{
#'   \item{poisson_models}{A list of fitted Poisson GLMMs post burn-in.}
#'   \item{binomial_models}{A list of fitted Binomial GLMMs post burn-in.}
#' }
#'
#' @details This function implements a form of multi-species occupancy-abundance
#' modeling using latent variables. At each iteration, species with `y > 0` are
#' assumed present, while those with `y == 0` are assigned occupancy probabilistically.
#'
#' @importFrom glmmTMB glmmTMB
#' @export
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
