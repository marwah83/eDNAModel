#' @title Predict from a fitted eDNAModel
#' @description Generate predicted values for abundance or occupancy using a fitted model.
#'
#' @param model A fitted model object from \code{\link{run_full_TMB}}.
#' @param newX A new covariate data frame for prediction.
#' @param formula The right-hand side of the formula used for the model (e.g., \code{~ Site}).
#' @param which Either \code{"abundance"} or \code{"occupancy"} indicating which model to predict from.
#' @param type Prediction scale: one of \code{"response"}, \code{"link"}, or \code{"terms"}.
#' @param level Confidence level for prediction intervals.
#' @param se Logical; should standard errors and intervals be returned?
#'
#' @return A data frame with predicted values and, if \code{se = TRUE}, standard errors and confidence intervals.
#' @export
predict_TMB <- function(model, newX, formula, which = c("abundance", "occupancy"),
                        type = c("response", "link", "terms"), level = 0.95, se = TRUE) {
  which <- match.arg(which)
  type <- match.arg(type)

  # Select coefficient matrix and parameter name
  coefs <- if (which == "abundance") {
    model$TMBobj$env$parList()$Ba
  } else {
    model$TMBobj$env$parList()$Bo
  }
  coef_name <- if (which == "abundance") "Ba" else "Bo"

  beta_mean <- rowMeans(coefs)

  # Build design matrix
  Xmat <- tryCatch(
    model.matrix(formula, newX),
    error = function(e) stop("Design matrix creation failed: ", e$message)
  )

  if (ncol(Xmat) != length(beta_mean)) {
    stop(sprintf("Design matrix (%d columns) does not match coefficient length (%d).",
                 ncol(Xmat), length(beta_mean)))
  }

  # Compute linear predictor
  linpred <- as.vector(Xmat %*% beta_mean)

  # Option: Return terms
  if (type == "terms") {
    term_matrix <- sweep(Xmat, 2, beta_mean, `*`)
    colnames(term_matrix) <- colnames(Xmat)
    return(as.data.frame(term_matrix))
  }

  # Inverse link function
  inv_link <- switch(which,
                     "abundance" = exp,
                     "occupancy" = function(x) 1 - plogis(x))

  # Response scale
  prediction <- if (type == "response") inv_link(linpred) else linpred

  # If no SE requested
  if (!se) return(data.frame(prediction = prediction))

  # Get standard errors via sdreport
  sdr <- TMB::sdreport(model$TMBobj)
  summary_sdr <- summary(sdr)

  idx <- which(rownames(summary_sdr) == coef_name)
  V <- summary_sdr[idx, 2]^2

  if (length(V) < length(beta_mean)) stop("SEs missing for some coefficients.")

  pred_se_link <- sqrt(rowSums((Xmat^2) %*% diag(V[1:length(beta_mean)])))

  # CI on link or response scale
  z <- qnorm(1 - (1 - level)/2)
  lower_link <- linpred - z * pred_se_link
  upper_link <- linpred + z * pred_se_link

  if (type == "response") {
    lower <- inv_link(lower_link)
    upper <- inv_link(upper_link)
  } else {
    lower <- lower_link
    upper <- upper_link
  }

  return(data.frame(
    prediction = prediction,
    se = pred_se_link,
    lower = lower,
    upper = upper
  ))
}


 

