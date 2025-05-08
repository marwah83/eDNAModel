#' @export
print.summary.eDNAModel <- function(x, ...) {
  # Print the main summary table
  cat("Model Coefficient Summary:\n")
  print(as.data.frame(x), row.names = FALSE)
  
  # Extract and print the explanation
  explanation <- attr(x, "explanation")
  if (!is.null(explanation)) {
    cat("\nLegend:\n")
    print(explanation, row.names = FALSE)
  } else {
    cat("\n(No explanation metadata available.)\n")
  }
  invisible(x)
}
