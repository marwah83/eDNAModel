#' Evaluate cross-validation results
#'
#' Summarises cross-validation results from
#' \code{compare_models_cv_gllvm()} or
#' \code{compare_models_cv()}.
#'
#' Computes:
#' \itemize{
#'   \item Mean occupancy AUC
#'   \item Mean capture AUC
#'   \item Mean predictive abundance log-likelihood
#' }
#'
#' Models are ranked by predictive abundance
#' log-likelihood (\code{lambda_log_score}),
#' where larger values are better.
#'
#' @param cv_res A list returned by
#' \code{compare_models_cv_gllvm()} or
#' \code{compare_models_cv()}.
#'
#' @return A data.frame summarising model performance.
#' @examples
#' \dontrun{
#'
#' results_table <- evaluate_cv(cv_res)
#'
#' print(results_table)
#'
#' }
#' @export
evaluate_cv <- function(cv_res) {

  if (!is.list(cv_res)) {
    stop("cv_res must be a list.")
  }

  results_table <- dplyr::bind_rows(

    lapply(
      names(cv_res),
      function(mod) {

        x <- cv_res[[mod]]

        data.frame(

          Model = mod,

          psi_AUC = mean(
            x$psi_AUC,
            na.rm = TRUE
          ),

          capture_AUC = mean(
            x$capture_AUC,
            na.rm = TRUE
          ),

          lambda_log_score = mean(
            x$lambda_log_score,
            na.rm = TRUE
          ),

          n_folds = nrow(x),

          n_valid_psi_AUC = sum(
            !is.na(x$psi_AUC)
          ),

          n_valid_capture_AUC = sum(
            !is.na(x$capture_AUC)
          )
        )
      }
    )
  )

  # ==========================================================
  # Ranking
  # ==========================================================

  results_table <- results_table |>

    dplyr::arrange(
      dplyr::desc(lambda_log_score)
    )

  rownames(results_table) <- NULL

  return(results_table)
}
