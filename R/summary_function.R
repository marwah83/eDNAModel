#' @title Summary of TMB Model Coefficients
#' @description Summarizes fixed-effect and random-effect parameters with standard deviations.
#'
#' @param model Output from \code{\link{run_full_TMB}}
#'
#' @return Data frame with parameter names, estimates, and standard deviations
#' @aliases summary.eDNAModel summary
#' @method summary eDNAModel
#' @export
summary.eDNAModel <- function(model) {
  # Get standard deviation report
  sdr <- TMB::sdreport(model$TMBobj)

  # Extract full summary matrix (estimates + sd + CI)
  full_summary <- summary(sdr)

  # Get parameter names from the rownames
  param_names <- rownames(full_summary)

  if (length(param_names) != nrow(full_summary)) {
    stop("Parameter names missing or not matching number of estimates.")
  }

  # Build summary table
  summary_df <- data.frame(
    parameter = param_names,
    estimate = full_summary[, 1],
    sd = full_summary[, 2],
    stringsAsFactors = FALSE
  )

  # Filter for model parameters of interest
  summary_df <- summary_df[grepl("Ba|Bo|Ua|Uo|logphi|logsda|logsdo", summary_df$parameter), ]

  class(summary_df) <- "summary.eDNAModel"

  # Return clean result
  return(summary_df)
}
