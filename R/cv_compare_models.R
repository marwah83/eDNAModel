#' Compare models using grouped K-fold cross-validation
#'
#' @param phyloseq A phyloseq object
#' @param models Named list of model functions
#' @param site_col Column defining grouping (e.g. site)
#' @param K Number of folds
#' @param seed Random seed
#' @param metric Evaluation metric
#' @param verbose Logical
#'
#' @return List with fold scores and summary
#' @export
compare_models_cv <- function(
  phyloseq,
  models,
  site_col,
  K = 5,
  seed = 123,
  metric = "log_score",
  verbose = TRUE
) {

  set.seed(seed)

  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col
  )

  long_df <- prep$long_df

  sites <- unique(long_df[[site_col]])
  folds <- sample(rep(1:K, length.out = length(sites)))

  fold_map <- data.frame(
    site = sites,
    fold = folds
  )
  names(fold_map)[1] <- site_col

  results <- list()

  for (m in seq_along(models)) {

    model_name <- names(models)[m]
    model_fun  <- models[[m]]

    if (verbose) message("\nModel: ", model_name)

    scores <- numeric(K)

    for (k in seq_len(K)) {

      test_sites <- fold_map[[site_col]][fold_map$fold == k]

      train_idx <- !(long_df[[site_col]] %in% test_sites)
      test_idx  <-  (long_df[[site_col]] %in% test_sites)

      test_df <- long_df[test_idx, ]

      fit <- model_fun(phyloseq)

      pred <- predict_FitModel(fit, test_df)

      scores[k] <- compute_metrics(test_df, pred, metric)

      if (verbose) {
        message("  Fold ", k, ": ", round(scores[k], 3))
      }
    }

    results[[model_name]] <- data.frame(
      Model = model_name,
      mean = mean(scores, na.rm = TRUE),
      sd   = sd(scores, na.rm = TRUE)
    )
  }

  dplyr::bind_rows(results) |>
    dplyr::arrange(dplyr::desc(mean))
}
