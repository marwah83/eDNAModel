#' Extract Fixed Effects from a List of glmmTMB Models
#'
#' This function extracts the fixed effect coefficients from a list of `glmmTMB` model objects,
#' adding metadata such as the iteration index and model name.
#'
#' @param model_list A list of fitted `glmmTMB` models (e.g., from `simulate_glm_burnin_iterations`).
#' @param model_name A character string to label the model type (e.g., `"poisson"` or `"binomial"`).
#'
#' @return A data frame with fixed effect estimates across all iterations. Columns include:
#' \describe{
#'   \item{Estimate, Std. Error, z value, Pr(>|z|)}{Standard coefficient outputs from `summary.glmmTMB`.}
#'   \item{term}{Name of the fixed effect term.}
#'   \item{iter}{The iteration index from the model list.}
#'   \item{model}{The model name provided.}
#' }
#'
#' @importFrom glmmTMB summary.glmmTMB
#' @export
extract_fixed <- function(model_list, model_name) {
  do.call(rbind, lapply(seq_along(model_list), function(i) {
    sm <- as.data.frame(summary(model_list[[i]])$coefficients$cond)
    sm$term <- rownames(sm)
    sm$iter <- i
    sm$model <- model_name
    sm
  }))
}
