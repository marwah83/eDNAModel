# ✅ Function to Compute Probabilities
compute_probabilities_wrapper <- function(fit_results) {
  if (is.null(fit_results)) return(NULL)

  tryCatch(
    compute_probabilities(fit_results),
    error = function(e) return(NULL)
  )
}
