#' @title Predict from TMB Model with Confidence Intervals (Safe Version)
#'
#' @param model Output from run_full_TMB()
#' @param newX New data.frame of covariates
#' @param formula The same formula used for abundance or occupancy
#' @param type "abundance" or "occupancy"
#' @param level Confidence level for CI
#'
#' @return Data frame with predictions and CIs
#' @export
predict_TMB <- function(model, newX, formula, type = c("abundance", "occupancy"), level = 0.95) {
  type <- match.arg(type)

  # Determine which coefficient matrix to use
  coef_matrix <- switch(type,
                        "abundance" = model$TMBobj$env$parList()$Ba,
                        "occupancy" = model$TMBobj$env$parList()$Bo
  )
  coef_name <- ifelse(type == "abundance", "Ba", "Bo")

  # Compute species-averaged coefficients
  beta_mean <- rowMeans(coef_matrix)

  # Ensure all factor levels in newX match those in training
  # This avoids mismatch when factor levels differ
  training_data <- model$TMBobj$env$data
  for (var in names(newX)) {
    if (is.factor(newX[[var]]) && var %in% names(training_data)) {
      ref_levels <- levels(as.factor(training_data[[var]]))
      newX[[var]] <- factor(newX[[var]], levels = ref_levels)
    }
  }

  # Construct design matrix
  Xmat <- tryCatch(
    model.matrix(formula, data = newX),
    error = function(e) stop("Design matrix creation failed: ", e$message)
  )

  # Ensure dimensions match
  if (ncol(Xmat) != length(beta_mean)) {
    stop(sprintf("Design matrix (%d columns) does not match coefficient length (%d).",
                 ncol(Xmat), length(beta_mean)))
  }

  # Prediction point estimates
  prediction <- as.vector(Xmat %*% beta_mean)

  # Confidence intervals from standard errors
  sdr <- TMB::sdreport(model$TMBobj)
  coef_rows <- which(rownames(summary(sdr)) == coef_name)
  coef_vars <- summary(sdr)[coef_rows, 2]^2

  if (length(coef_vars) < length(beta_mean)) {
    stop("Standard errors not available for all coefficients.")
  }

  # Compute standard error of predictions
  pred_se <- sqrt(rowSums((Xmat^2) %*% diag(coef_vars[1:length(beta_mean)])))

  # Confidence intervals
  z <- qnorm(1 - (1 - level) / 2)
  lower <- prediction - z * pred_se
  upper <- prediction + z * pred_se

  return(data.frame(
    prediction = prediction,
    se = pred_se,
    lower = lower,
    upper = upper
  ))
}



