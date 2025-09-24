#' Simulate ZIP (Zero-Inflated Poisson) Models with Burn-in
#'
#' This function performs iterative fitting of a zero-inflated Poisson model using
#' separate Poisson and Binomial components. It automatically constructs model formulas
#' by multiplying all fixed and random terms by the species variable (default: `"i"`),
#' while keeping offsets unmodified.
#'
#' @param data_glm A data frame in long format with count column `y`, OTU ID column (e.g., `i`),
#'        treatment, site, sample, replicate, and optionally an offset column.
#' @param poisson_rhs A quoted formula object representing the right-hand side of the Poisson model,
#'        e.g., `quote((1 | Site) + treatment + offset(log_total_reads))`.
#' @param binomial_rhs A quoted formula object representing the right-hand side of the Binomial model.
#' @param num_iterations Total number of simulation iterations (default: 100).
#' @param burn_in Number of burn-in iterations to discard before saving models (default: 50).
#' @param species_var A string name of the species (OTU) column to interact with (default: "i").
#'
#' @return A list with:
#' \describe{
#'   \item{poisson_models}{List of fitted Poisson `glmmTMB` model objects after burn-in.}
#'   \item{binomial_models}{List of fitted Binomial `glmmTMB` model objects after burn-in.}
#' }
#'
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats predict rbinom
#' @export
simulate_glm_burnin_iterations <- function(data_glm,
                                           poisson_rhs,
                                           binomial_rhs,
                                           num_iterations = 100,
                                           burn_in = 50,
                                           species_var = "i") {
  # Build full formulas
  poisson_formula <- build_formula_with_species_interaction(poisson_rhs, response = "y", species_var = species_var)
  binomial_formula <- build_formula_with_species_interaction(binomial_rhs, response = "z_sim", species_var = species_var)

  message("📌 Poisson model formula: ", deparse(poisson_formula))
  message("📌 Binomial model formula: ", deparse(binomial_formula))

  # Initialize latent occupancy variable
  data_glm$z_sim <- ifelse(data_glm$y > 0, 1, rbinom(nrow(data_glm), 1, 0.5))

  poisson_models <- list()
  binomial_models <- list()

  for (iter in 1:num_iterations) {
    # Subset data where occupancy = 1
    Q <- data_glm[data_glm$z_sim == 1, ]

    # Fit models
    model_poisson <- glmmTMB::glmmTMB(poisson_formula, family = poisson, data = Q)
    model_binomial <- glmmTMB::glmmTMB(binomial_formula, family = binomial, data = data_glm)

    # Predict lambda and occupancy
    lambda_i <- predict(model_poisson, type = "response", newdata = data_glm)
    P_i <- predict(model_binomial, type = "response", newdata = data_glm)

    # Update latent occupancy
    prob_Z1_given_y0 <- P_i * exp(-lambda_i) / (P_i * exp(-lambda_i) + (1 - P_i))
    prob_Z1_given_y0 <- pmin(pmax(prob_Z1_given_y0, 1e-6), 1 - 1e-6)
    data_glm$z_sim[data_glm$y == 0] <-
      rbinom(sum(data_glm$y == 0), 1, prob_Z1_given_y0[data_glm$y == 0])

    # Store models after burn-in
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
