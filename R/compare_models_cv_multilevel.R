#' Cross-validation for multilevel eDNA models
#'
#' Performs K-fold cross-validation for hierarchical occupancy–detection–abundance models
#' fitted using `FitModel()`. Supports multiple CV levels: site, sample, replicate, and OTU.
#'
#' @param phyloseq A phyloseq object.
#' @param models Named list of model functions. Each must take a phyloseq object and return a FitModel object.
#' @param site_col Column name for site.
#' @param sample_col Column name for biological sample.
#' @param replicate_col Column name for replicate (optional).
#' @param otu_col Column name for OTU.
#' @param count_col Column name for counts.
#' @param cv_level Cross-validation level: "site", "sample", "replicate", or "OTU".
#' @param K Number of folds.
#' @param seed Random seed.
#'
#' @return A named list of data frames with CV metrics per model.
#'
#' @export
compare_models_cv_multilevel <- function(
  phyloseq,
  models,
  site_col,
  sample_col,
  replicate_col = NULL,
  otu_col,
  count_col,
  cv_level = c("site", "sample", "OTU", "replicate"),
  K = 5,
  seed = 123
) {

  cv_level <- match.arg(cv_level)
  set.seed(seed)

  # -------------------------
  # Metadata
  # -------------------------
  meta <- as.data.frame(phyloseq::sample_data(phyloseq))
  meta$.__sample_id__ <- rownames(meta)

  if (!(sample_col %in% names(meta))) {
    meta[[sample_col]] <- meta$.__sample_id__
  }

  if (!(site_col %in% names(meta))) {
    stop(site_col, " not found in sample_data.")
  }

  if (cv_level == "replicate") {
    if (is.null(replicate_col)) stop("replicate_col required for replicate CV")
    if (!(replicate_col %in% names(meta))) {
      stop(replicate_col, " not found in sample_data.")
    }
  }

  # -------------------------
  # Define groups
  # -------------------------
  if (cv_level == "site") {
    groups <- unique(as.character(meta[[site_col]]))
  } else if (cv_level == "sample") {
    groups <- unique(as.character(meta[[sample_col]]))
  } else if (cv_level == "replicate") {
    groups <- unique(as.character(meta[[replicate_col]]))
  } else if (cv_level == "OTU") {
    groups <- phyloseq::taxa_names(phyloseq)
  }

  groups <- groups[!is.na(groups)]

  if (length(groups) < 2) {
    stop("Not enough groups for CV.")
  }

  K_eff <- min(K, length(groups))

  fold_ids <- sample(rep(seq_len(K_eff), length.out = length(groups)))
  folds <- split(groups, fold_ids)

  results <- list()

  # ============================================================
  # LOOP MODELS
  # ============================================================
  for (model_name in names(models)) {

    model_fun <- models[[model_name]]
    fold_metrics <- vector("list", length(folds))

    for (k in seq_along(folds)) {

      test_groups <- folds[[k]]

      # -------------------------
      # Split data
      # -------------------------
      if (cv_level == "OTU") {

        all_otus <- phyloseq::taxa_names(phyloseq)

        test_otus  <- intersect(test_groups, all_otus)
        train_otus <- setdiff(all_otus, test_otus)

        if (length(test_otus) == 0 || length(train_otus) == 0) {
          fold_metrics[[k]] <- data.frame(
            psi_AUC = NA_real_,
            capture_AUC = NA_real_,
            lambda_log_score = NA_real_
          )
          next
        }

        ps_train <- phyloseq::prune_taxa(train_otus, phyloseq)
        ps_test  <- phyloseq::prune_taxa(test_otus, phyloseq)

      } else {

        if (cv_level == "site") {
          test_samples <- meta$.__sample_id__[meta[[site_col]] %in% test_groups]
        } else if (cv_level == "sample") {
          test_samples <- meta$.__sample_id__[meta[[sample_col]] %in% test_groups]
        } else if (cv_level == "replicate") {
          test_samples <- meta$.__sample_id__[meta[[replicate_col]] %in% test_groups]
        }

        train_samples <- setdiff(meta$.__sample_id__, test_samples)

        if (length(test_samples) == 0 || length(train_samples) == 0) {
          fold_metrics[[k]] <- data.frame(
            psi_AUC = NA_real_,
            capture_AUC = NA_real_,
            lambda_log_score = NA_real_
          )
          next
        }

        ps_train <- phyloseq::prune_samples(train_samples, phyloseq)
        ps_test  <- phyloseq::prune_samples(test_samples, phyloseq)
      }

      # -------------------------
      # Fit model
      # -------------------------
      fit <- tryCatch(model_fun(ps_train), error = function(e) NULL)

      if (is.null(fit)) {
        fold_metrics[[k]] <- data.frame(
          psi_AUC = NA_real_,
          capture_AUC = NA_real_,
          lambda_log_score = NA_real_
        )
        next
      }

      # -------------------------
      # Prepare test data
      # -------------------------
      test_long <- tryCatch(
        prepare_long_data(
          physeq_obj = ps_test,
          site_col = site_col,
          nested_cols = unique(c(sample_col, replicate_col))
        )$long_df,
        error = function(e) NULL
      )

      if (is.null(test_long) || nrow(test_long) == 0) {
        fold_metrics[[k]] <- data.frame(
          psi_AUC = NA_real_,
          capture_AUC = NA_real_,
          lambda_log_score = NA_real_
        )
        next
      }

      # -------------------------
      # Offset fix
      # -------------------------
      if (!("total_reads" %in% names(test_long))) {
        total_reads_df <- test_long |>
          dplyr::group_by(.data[[sample_col]]) |>
          dplyr::summarise(
            total_reads = sum(.data[[count_col]], na.rm = TRUE),
            .groups = "drop"
          )

        test_long <- dplyr::left_join(test_long, total_reads_df, by = sample_col)
      }

      test_long$total_reads <- pmax(test_long$total_reads, 1)

      # -------------------------
      # Truth
      # -------------------------
      site_eval <- test_long |>
        dplyr::group_by(.data[[site_col]], .data[[otu_col]]) |>
        dplyr::summarise(z = as.integer(any(.data[[count_col]] > 0)), .groups = "drop")

      sample_eval <- test_long |>
        dplyr::group_by(.data[[site_col]], .data[[sample_col]], .data[[otu_col]]) |>
        dplyr::summarise(a = as.integer(any(.data[[count_col]] > 0)), .groups = "drop")

      # -------------------------
      # Extract models
      # -------------------------
      occ_fit   <- tail(fit$occupancy_models, 1)[[1]]
      cap_fit   <- tail(fit$capture_models, 1)[[1]]
      abund_fit <- tail(fit$abundance_models, 1)[[1]]

      # -------------------------
      # Predictions
      # -------------------------
      site_eval$psi <- suppressWarnings(
        as.numeric(predict(occ_fit, site_eval, type = "response", allow.new.levels = TRUE))
      )

      sample_eval$capture <- suppressWarnings(
        as.numeric(predict(cap_fit, sample_eval, type = "response", allow.new.levels = TRUE))
      )

      test_long$lambda <- suppressWarnings(
        as.numeric(predict(abund_fit, test_long, type = "response", allow.new.levels = TRUE))
      )

      # -------------------------
      # Metrics
      # -------------------------
      psi_auc <- if (length(unique(site_eval$z)) > 1) {
        as.numeric(pROC::roc(site_eval$z, site_eval$psi, quiet = TRUE)$auc)
      } else NA_real_

      capture_auc <- if (length(unique(sample_eval$a)) > 1) {
        as.numeric(pROC::roc(sample_eval$a, sample_eval$capture, quiet = TRUE)$auc)
      } else NA_real_

      lambda_log_score <- mean(
        stats::dpois(
          test_long[[count_col]],
          lambda = pmax(test_long$lambda, 1e-8),
          log = TRUE
        ),
        na.rm = TRUE
      )

      fold_metrics[[k]] <- data.frame(
        psi_AUC = psi_auc,
        capture_AUC = capture_auc,
        lambda_log_score = lambda_log_score
      )
    }

    results[[model_name]] <- dplyr::bind_rows(fold_metrics)
  }

  results
}
