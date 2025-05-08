#' @title Summary of TMB Model Coefficients
#' @description Summarizes fixed-effect and random-effect parameters with standard deviations.
#'
#' @param model Output from \code{\link{run_full_TMB}}
#' @param explain Logical; if TRUE (default), prints descriptions of each parameter type.
#'
#' @return A data frame with parameter names, estimates, and standard deviations
#' @aliases summary.eDNAModel summary
#' @method summary eDNAModel
#' @export
summary.eDNAModel <- function(model,explain = TRUE) {
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

  if (explain) {
    cat("\n🧾 Parameter Description:\n")
    cat("---------------------------------------------------\n")
    cat("Ba      - Fixed effects for abundance (log scale)\n")
    cat("Bo      - Fixed effects for occupancy (logit/probit scale)\n")
    cat("Ua      - Random effects for abundance (e.g., replicate-level)\n")
    cat("Uo      - Random effects for occupancy\n")
    cat("logphi  - Dispersion parameter (ZINB only)\n")
    cat("logsda  - Log standard deviation of abundance random effects\n")
    cat("logsdo  - Log standard deviation of occupancy random effects\n")
    cat("---------------------------------------------------\n\n")
  }


  # Return clean result
  return(summary_df)
}
