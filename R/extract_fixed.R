extract_fixed <- function(model_list, model_name) {
  do.call(rbind, lapply(seq_along(model_list), function(i) {
    sm <- as.data.frame(summary(model_list[[i]])$coefficients$cond)
    sm$term <- rownames(sm)
    sm$iter <- i
    sm$model <- model_name
    sm
  }))
}
