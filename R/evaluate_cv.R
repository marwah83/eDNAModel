#' Summarise cross-validation results
#'
#' Aggregates fold-level CV results (AUC and log score) into model-level summaries.
#'
#' @param cv_res A named list returned by \code{compare_models_cv_multilevel()}.
#'
#' @return A data.frame with one row per model and averaged metrics:
#' \itemize{
#'   \item psi_AUC
#'   \item capture_AUC
#'   \item lambda_log_score
#' }
#'
#' @export
evaluate_cv <- function(cv_res) {
  
  # -------------------------
  # Input check
  # -------------------------
  if (!is.list(cv_res) || length(cv_res) == 0) {
    stop("cv_res must be a non-empty list.")
  }
  
  # -------------------------
  # Loop over models
  # -------------------------
  out_list <- lapply(names(cv_res), function(model_name) {
    
    df <- cv_res[[model_name]]
    
    # -------------------------
    # Safety checks
    # -------------------------
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame(
        Model = model_name,
        psi_AUC = NA_real_,
        capture_AUC = NA_real_,
        lambda_log_score = NA_real_
      ))
    }
    
    # Ensure columns exist
    if (!all(c("psi_AUC", "capture_AUC", "lambda_log_score") %in% names(df))) {
      warning("Missing metrics for model: ", model_name)
      return(data.frame(
        Model = model_name,
        psi_AUC = NA_real_,
        capture_AUC = NA_real_,
        lambda_log_score = NA_real_
      ))
    }
    
    # -------------------------
    # Compute means (robust)
    # -------------------------
    data.frame(
      Model = model_name,
      
      psi_AUC = if (all(is.na(df$psi_AUC))) {
        NA_real_
      } else {
        mean(df$psi_AUC, na.rm = TRUE)
      },
      
      capture_AUC = if (all(is.na(df$capture_AUC))) {
        NA_real_
      } else {
        mean(df$capture_AUC, na.rm = TRUE)
      },
      
      lambda_log_score = if (all(is.na(df$lambda_log_score))) {
        NA_real_
      } else {
        mean(df$lambda_log_score, na.rm = TRUE)
      }
    )
  })
  
  # -------------------------
  # Combine
  # -------------------------
  out <- do.call(rbind, out_list)
  
  rownames(out) <- NULL
  
  return(out)
}
