compute_metrics <- function(test_df, pred_df, metric = "log_score") {

  y <- test_df$y
  lambda <- pred_df$lambda_mean

  if (metric == "log_score") {
    return(mean(dpois(y, lambda = lambda, log = TRUE), na.rm = TRUE))
  }

  if (metric == "rmse") {
    return(sqrt(mean((y - lambda)^2, na.rm = TRUE)))
  }

  if (metric == "mae") {
    return(mean(abs(y - lambda), na.rm = TRUE))
  }

  if (metric == "brier") {
    y_bin <- as.integer(y > 0)
    p <- 1 - exp(-lambda)
    return(mean((y_bin - p)^2, na.rm = TRUE))
  }

  stop("Unknown metric.")
}
