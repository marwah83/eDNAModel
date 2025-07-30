#' @title Summary of TMB Model Coefficients
#' @description Summarizes fixed-effect and random-effect parameters with standard deviations.
#'
#' @param object Output from \code{\link{run_full_TMB}} (an object of class \code{eDNAModel})
#' @param getJointPrecision \code{logical}, defaults to FALSE.
#' @param ... Additional arguments passed to TMB's summary function.
#'
#' @return A data frame with parameter names, estimates, and standard deviations.
#' @aliases summary.eDNAModel summary
#' @method summary eDNAModel
#' @export
summary.eDNAModel <- function(object, ...) {
  model <- object

  # Get standard deviation report
  sdr <- TMB::sdreport(model$TMBobj, ...)

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

  # Add explanation table as attribute
  explanation <- data.frame(
    parameter_type = c("Ba", "Bo", "Ua", "Uo", "logphi", "logsda", "logsdo"),
    meaning = c(
      "Abundance model fixed effects",
      "Occupancy model fixed effects",
      "Random effects for abundance (e.g. Replicates)",
      "Random effects for occupancy (if any)",
      "Species-specific dispersion (ZINB models)",
      "Log standard deviation of abundance random effects",
      "Log standard deviation of occupancy random effects"
    ),
    stringsAsFactors = FALSE
  )

  attr(summary_df, "explanation") <- explanation

  # Set dual class
  class(summary_df) <- c("summary.eDNAModel", "data.frame")

  return(summary_df)
}
