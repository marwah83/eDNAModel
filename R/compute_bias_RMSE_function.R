# ✅ Function to Compute Bias and RMSE
compute_bias_rmse <- function(true_occ_prob, est_occ_prob, true_prob_detect, est_prob_detect) {

  # ✅ Fix Mismatched Dimensions
  if (!all(dim(true_occ_prob) == dim(est_occ_prob))) {
    est_occ_prob <- t(est_occ_prob)
  }
  if (!all(dim(true_prob_detect) == dim(est_prob_detect))) {
    est_prob_detect <- t(est_prob_detect)
  }

  # ✅ Compute Bias and RMSE
  bias_occ <- tryCatch(compute_bias(as.vector(true_occ_prob), as.vector(est_occ_prob)), error = function(e) NA)
  rmse_occ <- tryCatch(compute_rmse(as.vector(true_occ_prob), as.vector(est_occ_prob)), error = function(e) NA)

  bias_detect <- tryCatch(compute_bias(as.vector(true_prob_detect), as.vector(est_prob_detect)), error = function(e) NA)
  rmse_detect <- tryCatch(compute_rmse(as.vector(true_prob_detect), as.vector(est_prob_detect)), error = function(e) NA)

  return(list(
    bias_occ = bias_occ, rmse_occ = rmse_occ,
    bias_detect = bias_detect, rmse_detect = rmse_detect
  ))
}
