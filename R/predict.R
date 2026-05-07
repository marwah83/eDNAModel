predict_FitModel <- function(fit, newdata) {

  # Merge predictions with newdata keys
  pred <- fit$lambda

  keys <- intersect(names(pred), names(newdata))

  out <- dplyr::left_join(newdata, pred, by = keys)

  # fallback if missing
  if (!("lambda_mean" %in% names(out))) {
    stop("lambda_mean missing in predictions")
  }

  out
}
