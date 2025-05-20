#' @title Print Method for summary.eDNAModel Objects
#' @description Custom print method for objects of class \code{summary.eDNAModel}.
#' Prints the parameter estimates with standard deviations and explains parameter types.
#'
#' @param x A \code{summary.eDNAModel} object returned by \code{summary.eDNAModel}.
#' @param ... Additional arguments (currently ignored).
#' @export
print.summary.eDNAModel <- function(x, ...) {
  cat("\nModel Parameter Summary:\n")
  print.data.frame(x, row.names = FALSE)

  if (!is.null(attr(x, "explanation"))) {
    cat("\nLegend for parameter types:\n")
    print(attr(x, "explanation"), row.names = FALSE)
  }

  invisible(x)
}
