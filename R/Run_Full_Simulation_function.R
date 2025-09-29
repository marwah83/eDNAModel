# ✅ Function to Run Full Simulation
run_full_simulation <- function(species, sites, replicates, lambda, ZIP) {
  set.seed(123)  # Ensure reproducibility

  # ✅ Step 1: Create Zero-Inflation Matrix
  zero_inflation_prob_matrix <- create_zero_inflation_matrix(species, sites, ZIP)

  # ✅ Step 2: Simulate Data
  simulated_data <- simulate_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)

  # ✅ Step 3: Check Zero Proportion
  zero <- mean(simulated_data == 0)

  # ✅ Step 4: Process Data and Fit Model
  fit_results <- process_and_fit_data(simulated_data)
  if (is.null(fit_results)) {
    return(data.frame(
      Species = species, Sites = sites, Replicates = replicates,
      Lambda = lambda, ZIP = ZIP,
      Mean_Occ_Prob = NA, Mean_Prob_Detect = NA,
      Bias_Occupancy = NA, RMSE_Occupancy = NA,
      Bias_Detection = NA, RMSE_Detection = NA,
      zero = zero
    ))
  }

  # ✅ Step 5: Compute Probabilities
  probabilities <- compute_probabilities_wrapper(fit_results)

  # ✅ Step 6: Compute Bias and RMSE
  true_occ_prob <- 1 - zero_inflation_prob_matrix
  est_occ_prob <- probabilities$occ.prob

  true_prob_detect <- 1 - true_occ_prob * exp(-lambda)
  est_prob_detect <- probabilities$prob.detect

  bias_rmse_results <- compute_bias_rmse(true_occ_prob, est_occ_prob, true_prob_detect, est_prob_detect)

  # ✅ Step 7: Return Simulation Results
  return(data.frame(
    Species = species, Sites = sites, Replicates = replicates,
    Lambda = lambda, ZIP = ZIP,
    Mean_Occ_Prob = mean(est_occ_prob, na.rm = TRUE),
    Mean_Prob_Detect = mean(est_prob_detect, na.rm = TRUE),
    Bias_Occupancy = bias_rmse_results$bias_occ, RMSE_Occupancy = bias_rmse_results$rmse_occ,
    Bias_Detection = bias_rmse_results$bias_detect, RMSE_Detection = bias_rmse_results$rmse_detect,
    zero = zero
  ))
}
