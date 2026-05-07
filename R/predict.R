predict_FitModel <- function(object, newdata) {

  if (is.null(object$lambda)) {
    stop("FitModel object must contain lambda predictions.")
  }

  lambda <- object$lambda

  # Join keys (adjust if needed)
  keys <- intersect(names(newdata), names(lambda))

  pred <- dplyr::left_join(newdata, lambda, by = keys)

  # fallback
  if ("lambda_mean" %in% names(pred)) {
    pred$lambda_mean[is.na(pred$lambda_mean)] <-
      mean(lambda$lambda_mean, na.rm = TRUE)
  }

  pred
}
