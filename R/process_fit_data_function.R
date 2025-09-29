# ✅ Function to Process and Fit Data
process_and_fit_data <- function(simulated_data) {
  if (is.null(simulated_data)) return(NULL)

  Y <- to2D(simulated_data)  # Convert to 2D format

  fit_results <- tryCatch(
    prepare_tmb_data_fit(Y),
    error = function(e) return(NULL)
  )

  return(fit_results)
}
