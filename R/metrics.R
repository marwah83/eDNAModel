compute_metrics <- function(obs_df, pred_df, metric = "log_score") {

  if (!all(c("y") %in% names(obs_df))) {
    stop("Observed data must contain column 'y'")
  }

  if (!("lambda_mean" %in% names(pred_df))) {
    stop("Prediction must contain 'lambda_mean'")
  }

  y <- obs_df$y
  mu <- pred_df$lambda_mean

  if (metric == "log_score") {
    # Poisson log predictive density
    return(mean(stats::dpois(y, lambda = mu, log = TRUE), na.rm = TRUE))
  }

  if (metric == "rmse") {
    return(sqrt(mean((y - mu)^2, na.rm = TRUE)))
  }

  stop("Unknown metric")
}
