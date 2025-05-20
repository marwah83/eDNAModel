#' @title Print Method for summary.eDNAModel Objects
#' @description Custom print method for summary results from \code{eDNAModel}.
#'
#' @param x A \code{summary.eDNAModel} object returned by \code{\link{summary.eDNAModel}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @seealso \code{\link{summary.eDNAModel}}
#' @method print summary.eDNAModel
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
