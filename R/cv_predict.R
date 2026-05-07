#' Predict from FitModel object
#'
#' @param fit Output from FitModel
#' @param newdata Data frame for prediction
#'
#' @return Data frame with predictions
#' @export
predict_FitModel <- function(fit, newdata) {

  pred <- fit$lambda

  keys <- intersect(names(pred), names(newdata))

  out <- dplyr::left_join(newdata, pred, by = keys)

  if (!("lambda_mean" %in% names(out))) {
    stop("lambda_mean missing in predictions.")
  }

  out
}
