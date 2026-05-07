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

  if (!(site_col %in% names(long_df))) {
    stop("site_col not found in long_df.")
  }

  sites <- unique(long_df[[site_col]])

  if (length(sites) < K) {
    stop("Number of sites < K.")
  }

  folds <- sample(rep(seq_len(K), length.out = length(sites)))

  fold_map <- data.frame(site = sites, fold = folds)
  names(fold_map)[1] <- site_col

  results <- list()

  for (m in seq_along(models)) {

    model_name <- names(models)[m]
    model_fun  <- models[[m]]

    if (verbose) message("\nModel: ", model_name)

    fold_results <- list()

    for (k in seq_len(K)) {

      test_sites  <- fold_map[[site_col]][fold_map$fold == k]
      train_sites <- fold_map[[site_col]][fold_map$fold != k]

      train_phyloseq <- subset_phyloseq_by_sites(
        physeq = phyloseq,
        site_col = site_col,
        keep_sites = train_sites
      )

      test_df <- long_df[long_df[[site_col]] %in% test_sites, , drop = FALSE]

      fit <- model_fun(train_phyloseq)

      pred <- predict_FitModel(fit, test_df)

      score <- compute_metrics(test_df, pred, metric)

      fold_results[[k]] <- data.frame(
        Model = model_name,
        Fold  = k,
        Score = score
      )

      if (verbose) {
        message("  Fold ", k, ": ", round(score, 4))
      }
    }

    results[[model_name]] <- dplyr::bind_rows(fold_results)
  }

  all_results <- dplyr::bind_rows(results)

  summary <- all_results |>
    dplyr::group_by(Model) |>
    dplyr::summarise(
      mean_score = mean(Score, na.rm = TRUE),
      sd_score   = sd(Score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(mean_score))

  list(
    fold_scores = all_results,
    summary = summary,
    fold_map = fold_map,
    metric = metric,
    cv_type = "Grouped K-fold CV (site-level)"
  )
}
