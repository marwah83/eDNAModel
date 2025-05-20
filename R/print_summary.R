#' @title Print Method for Summary of eDNAModel
#' @description Print method for objects of class \code{summary.eDNAModel}.
#' @param x An object of class \code{summary.eDNAModel}.
#' @param ... Additional arguments (currently ignored).
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
