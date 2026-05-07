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

  # ---------------------------------------
  # Prepare data
  # ---------------------------------------
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col
  )

  long_df <- prep$long_df

  if (!(site_col %in% names(long_df))) {
    stop("site_col not found in long_df")
  }

  sites <- unique(long_df[[site_col]])

  if (length(sites) < K) {
    stop("Number of sites < K folds")
  }

  # ---------------------------------------
  # Create grouped folds
  # ---------------------------------------
  folds <- sample(rep(seq_len(K), length.out = length(sites)))

  fold_map <- data.frame(
    site = sites,
    fold = folds
  )
  names(fold_map)[1] <- site_col

  results <- list()

  # ---------------------------------------
  # Loop over models
  # ---------------------------------------
  for (m in seq_along(models)) {

    model_name <- names(models)[m]
    model_fun  <- models[[m]]

    if (verbose) {
      message("\n==============================")
      message("Model: ", model_name)
      message("==============================")
    }

    fold_results <- list()

    # ---------------------------------------
    # Loop over folds
    # ---------------------------------------
    for (k in seq_len(K)) {

      test_sites  <- fold_map[[site_col]][fold_map$fold == k]
      train_sites <- fold_map[[site_col]][fold_map$fold != k]

      # -------------------------------
      # TRAIN DATA (phyloseq subset)
      # -------------------------------
      train_phyloseq <- subset_phyloseq_by_sites(
        physeq = phyloseq,
        site_col = site_col,
        keep_sites = train_sites
      )

      # -------------------------------
      # TEST DATA (long_df subset)
      # -------------------------------
      test_df <- long_df[long_df[[site_col]] %in% test_sites, , drop = FALSE]

      # -------------------------------
      # Fit model on TRAIN only
      # -------------------------------
      fit <- model_fun(train_phyloseq)

      # -------------------------------
      # Predict on TEST
      # -------------------------------
      pred <- predict_FitModel(fit, test_df)

      # -------------------------------
      # Compute metric
      # -------------------------------
      score <- compute_metrics(test_df, pred, metric)

      fold_results[[k]] <- data.frame(
        Model = model_name,
        Fold  = k,
        Score = score
      )

      if (verbose) {
        message("  Fold ", k, " → ", round(score, 4))
      }
    }

    results[[model_name]] <- dplyr::bind_rows(fold_results)
  }

  # ---------------------------------------
  # Combine results
  # ---------------------------------------
  all_results <- dplyr::bind_rows(results)

  summary <- all_results |>
    dplyr::group_by(Model) |>
    dplyr::summarise(
      mean_score = mean(Score, na.rm = TRUE),
      sd_score   = sd(Score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(mean_score))

  return(list(
    fold_scores = all_results,
    summary = summary,
    fold_map = fold_map,
    metric = metric,
    cv_type = "Grouped K-fold CV (site-level)"
  ))
}
