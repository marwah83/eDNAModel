# ✅ Function to Compute RMSE
compute_rmse <- function(actual, estimated) {
  sqrt(mean((actual - estimated) ^ 2, na.rm = TRUE))
}
