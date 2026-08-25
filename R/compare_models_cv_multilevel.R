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

  # ============================================================
  # Helper functions
  # ============================================================

  safe_auc <- function(truth, score) {

    ok <- is.finite(score) & !is.na(truth)

    truth <- truth[ok]
    score <- score[ok]

    if (length(truth) == 0 || length(unique(truth)) < 2) {
      return(NA_real_)
    }

    tryCatch(
      as.numeric(
        pROC::roc(
          response = truth,
          predictor = score,
          quiet = TRUE
        )$auc
      ),
      error = function(e) NA_real_
    )
  }

  # ============================================================
  # Metadata
  # ============================================================

  meta <- as.data.frame(
    phyloseq::sample_data(phyloseq)
  )

  meta$.__sample_id__ <- rownames(meta)

  if (!(sample_col %in% names(meta))) {
    meta[[sample_col]] <- meta$.__sample_id__
  }

  if (!(site_col %in% names(meta))) {
    stop(site_col, " not found in sample_data.")
  }

  if (cv_level == "replicate") {

    if (is.null(replicate_col)) {
      stop("replicate_col required for replicate CV.")
    }

    if (!(replicate_col %in% names(meta))) {
      stop(replicate_col, " not found in sample_data.")
    }
  }

  # ============================================================
  # Define CV groups
  # ============================================================

  if (cv_level == "site") {

    groups <- unique(
      as.character(meta[[site_col]])
    )

  } else if (cv_level == "sample") {

    groups <- unique(
      as.character(meta[[sample_col]])
    )

  } else if (cv_level == "replicate") {

    groups <- unique(
      as.character(meta[[replicate_col]])
    )

  } else if (cv_level == "OTU") {

    groups <- phyloseq::taxa_names(phyloseq)
  }

  groups <- groups[!is.na(groups)]

  if (length(groups) < 2) {
    stop("Not enough groups for CV.")
  }

  K_eff <- min(K, length(groups))

  fold_ids <- sample(
    rep(seq_len(K_eff), length.out = length(groups))
  )

  folds <- split(groups, fold_ids)

  results <- list()

  # ============================================================
  # LOOP OVER MODELS
  # ============================================================

  for (model_name in names(models)) {

    cat("\n==============================\n")
    cat("Model:", model_name, "\n")
    cat("==============================\n")

    model_fun <- models[[model_name]]

    fold_metrics <- vector(
      "list",
      length(folds)
    )

    # ==========================================================
    # LOOP OVER FOLDS
    # ==========================================================

    for (k in seq_along(folds)) {

      cat(
        "\n--- Fold",
        k,
        "/",
        length(folds),
        "---\n"
      )

      test_groups <- folds[[k]]

      # ========================================================
      # Split data
      # ========================================================

      if (cv_level == "OTU") {

        all_otus <- as.character(
          phyloseq::taxa_names(phyloseq)
        )

        test_otus <- intersect(
          as.character(test_groups),
          all_otus
        )

        train_otus <- setdiff(
          all_otus,
          test_otus
        )

        if (length(test_otus) == 0 ||
            length(train_otus) == 0) {

          fold_metrics[[k]] <- data.frame(
            fold = k,
            psi_AUC = NA_real_,
            capture_AUC = NA_real_,
            lambda_log_score = NA_real_
          )

          next
        }

        ps_train <- phyloseq::prune_taxa(
          as.character(train_otus),
          phyloseq
        )

        ps_test <- phyloseq::prune_taxa(
          as.character(test_otus),
          phyloseq
        )

      } else {

        if (cv_level == "site") {

          test_samples <- meta$.__sample_id__[
            as.character(meta[[site_col]]) %in%
              as.character(test_groups)
          ]

        } else if (cv_level == "sample") {

          test_samples <- meta$.__sample_id__[
            as.character(meta[[sample_col]]) %in%
              as.character(test_groups)
          ]

        } else if (cv_level == "replicate") {

          test_samples <- meta$.__sample_id__[
            as.character(meta[[replicate_col]]) %in%
              as.character(test_groups)
          ]
        }

        train_samples <- setdiff(
          meta$.__sample_id__,
          test_samples
        )

        if (length(test_samples) == 0 ||
            length(train_samples) == 0) {

          fold_metrics[[k]] <- data.frame(
            fold = k,
            psi_AUC = NA_real_,
            capture_AUC = NA_real_,
            lambda_log_score = NA_real_
          )

          next
        }

        ps_train <- phyloseq::prune_samples(
          train_samples,
          phyloseq
        )

        ps_test <- phyloseq::prune_samples(
          test_samples,
          phyloseq
        )
      }

      # ========================================================
      # Fit model
      # ========================================================

      fit <- tryCatch(
        model_fun(ps_train),
        error = function(e) {

          message(
            "Model failed in fold ",
            k,
            ": ",
            e$message
          )

          NULL
        }
      )

      if (is.null(fit)) {

        fold_metrics[[k]] <- data.frame(
          fold = k,
          psi_AUC = NA_real_,
          capture_AUC = NA_real_,
          lambda_log_score = NA_real_
        )

        next
      }

      # ========================================================
      # Prepare test data
      # ========================================================

      test_long <- tryCatch(
        prepare_long_data(
          physeq_obj = ps_test,
          site_col = site_col,
          nested_cols = unique(
            c(sample_col, replicate_col)
          )
        )$long_df,
        error = function(e) {

          message(
            "prepare_long_data failed in fold ",
            k,
            ": ",
            e$message
          )

          NULL
        }
      )

      if (is.null(test_long) ||
          nrow(test_long) == 0) {

        fold_metrics[[k]] <- data.frame(
          fold = k,
          psi_AUC = NA_real_,
          capture_AUC = NA_real_,
          lambda_log_score = NA_real_
        )

        next
      }

      # ========================================================
      # Type handling
      # ========================================================

      test_long[[site_col]] <- as.character(
        test_long[[site_col]]
      )

      test_long[[sample_col]] <- as.character(
        test_long[[sample_col]]
      )

      test_long[[otu_col]] <- as.character(
        test_long[[otu_col]]
      )

      test_long[[count_col]] <- as.numeric(
        test_long[[count_col]]
      )

      # ========================================================
      # Offset
      # ========================================================

      if (!("total_reads" %in% names(test_long))) {

        total_reads_df <- test_long |>
          dplyr::group_by(
            .data[[sample_col]]
          ) |>
          dplyr::summarise(
            total_reads = sum(
              .data[[count_col]],
              na.rm = TRUE
            ),
            .groups = "drop"
          )

        test_long <- dplyr::left_join(
          test_long,
          total_reads_df,
          by = sample_col
        )
      }

      test_long$total_reads <- pmax(
        test_long$total_reads,
        1
      )

      # ========================================================
      # Truth
      # ========================================================

      site_eval <- test_long |>
        dplyr::group_by(
          .data[[site_col]],
          .data[[otu_col]]
        ) |>
        dplyr::summarise(
          z = as.integer(
            any(.data[[count_col]] > 0)
          ),
          .groups = "drop"
        )

      sample_eval <- test_long |>
        dplyr::group_by(
          .data[[site_col]],
          .data[[sample_col]],
          .data[[otu_col]]
        ) |>
        dplyr::summarise(
          a = as.integer(
            any(.data[[count_col]] > 0)
          ),
          .groups = "drop"
        )

      # ========================================================
      # Extract final retained models
      # ========================================================

      occ_fit <- tryCatch(
        tail(
          Filter(Negate(is.null), fit$occupancy_models),
          1
        )[[1]],
        error = function(e) NULL
      )

      cap_fit <- tryCatch(
        tail(
          Filter(Negate(is.null), fit$capture_models),
          1
        )[[1]],
        error = function(e) NULL
      )

      abund_fit <- tryCatch(
        tail(
          Filter(Negate(is.null), fit$abundance_models),
          1
        )[[1]],
        error = function(e) NULL
      )

      if (is.null(occ_fit) ||
          is.null(cap_fit) ||
          is.null(abund_fit)) {

        message(
          "Missing fitted component in fold ",
          k
        )

        fold_metrics[[k]] <- data.frame(
          fold = k,
          psi_AUC = NA_real_,
          capture_AUC = NA_real_,
          lambda_log_score = NA_real_
        )

        next
      }

      # ========================================================
      # Predictions
      # ========================================================

      site_eval$psi <- tryCatch(
        as.numeric(
          predict(
            occ_fit,
            newdata = site_eval,
            type = "response",
            allow.new.levels = TRUE
          )
        ),
        error = function(e) {
          rep(NA_real_, nrow(site_eval))
        }
      )

      sample_eval$capture <- tryCatch(
        as.numeric(
          predict(
            cap_fit,
            newdata = sample_eval,
            type = "response",
            allow.new.levels = TRUE
          )
        ),
        error = function(e) {
          rep(NA_real_, nrow(sample_eval))
        }
      )

      test_long$lambda <- tryCatch(
        as.numeric(
          predict(
            abund_fit,
            newdata = test_long,
            type = "response",
            allow.new.levels = TRUE
          )
        ),
        error = function(e) {
          rep(NA_real_, nrow(test_long))
        }
      )

      # ========================================================
      # AUC metrics
      # ========================================================

      psi_auc <- safe_auc(
        site_eval$z,
        site_eval$psi
      )

      capture_auc <- safe_auc(
        sample_eval$a,
        sample_eval$capture
      )

      # ========================================================
      # Predictive abundance log-likelihood
      # ========================================================

      y_obs <- as.numeric(
        test_long[[count_col]]
      )

      mu <- pmax(
        as.numeric(test_long$lambda),
        1e-8
      )

      family_name <- tryCatch(
        abund_fit$modelInfo$family$family,
        error = function(e) NA_character_
      )

      if (is.na(family_name)) {
        stop(
          "Could not determine abundance family in fold ",
          k
        )
      }

      # --------------------------------------------------------
      # Zero inflation probability
      # --------------------------------------------------------

      zi_prob <- tryCatch(
        as.numeric(
          predict(
            abund_fit,
            newdata = test_long,
            type = "zprob",
            allow.new.levels = TRUE
          )
        ),
        error = function(e) {
          rep(0, length(mu))
        }
      )

      zi_prob <- pmin(
        pmax(
          zi_prob,
          0
        ),
        1 - 1e-12
      )

      # --------------------------------------------------------
      # Poisson / ZIP
      # --------------------------------------------------------

      if (identical(
        family_name,
        "poisson"
      )) {

        if (all(
          zi_prob < 1e-12,
          na.rm = TRUE
        )) {

          loglik_vec <- stats::dpois(
            y_obs,
            lambda = mu,
            log = TRUE
          )

        } else {

          p0_pois <- stats::dpois(
            0,
            lambda = mu
          )

          loglik_vec <- ifelse(
            y_obs == 0,

            log(
              pmax(
                zi_prob +
                  (1 - zi_prob) * p0_pois,
                1e-300
              )
            ),

            log(
              pmax(
                1 - zi_prob,
                1e-300
              )
            ) +
              stats::dpois(
                y_obs,
                lambda = mu,
                log = TRUE
              )
          )
        }

      # --------------------------------------------------------
      # NB1 / NB2 / ZINB
      # --------------------------------------------------------

      } else if (
        family_name %in%
        c("nbinom1", "nbinom2")
      ) {

        theta <- tryCatch(
          as.numeric(
            stats::sigma(abund_fit)
          ),
          error = function(e) NA_real_
        )

        if (!is.finite(theta) ||
            theta <= 0) {

          stop(
            "Invalid negative-binomial dispersion in fold ",
            k,
            ". theta = ",
            theta
          )
        }

        if (all(
          zi_prob < 1e-12,
          na.rm = TRUE
        )) {

          loglik_vec <- stats::dnbinom(
            y_obs,
            mu = mu,
            size = theta,
            log = TRUE
          )

        } else {

          p0_nb <- stats::dnbinom(
            0,
            mu = mu,
            size = theta
          )

          loglik_vec <- ifelse(
            y_obs == 0,

            log(
              pmax(
                zi_prob +
                  (1 - zi_prob) * p0_nb,
                1e-300
              )
            ),

            log(
              pmax(
                1 - zi_prob,
                1e-300
              )
            ) +
              stats::dnbinom(
                y_obs,
                mu = mu,
                size = theta,
                log = TRUE
              )
          )
        }

      } else {

        stop(
          "Unsupported abundance family in CV scoring: ",
          family_name
        )
      }

      # ========================================================
      # Final log score
      # ========================================================

      loglik_vec[
        !is.finite(loglik_vec)
      ] <- NA_real_

      lambda_log_score <- if (
        all(is.na(loglik_vec))
      ) {

        NA_real_

      } else {

        mean(
          loglik_vec,
          na.rm = TRUE
        )
      }

      # ========================================================
      # Diagnostics
      # ========================================================

      cat(
        "Fold diagnostics:",
        "\n  family =", family_name,
        "\n  test rows =", nrow(test_long),
        "\n  observed y min/median/max =",
        paste(
          round(
            c(
              min(y_obs, na.rm = TRUE),
              stats::median(y_obs, na.rm = TRUE),
              max(y_obs, na.rm = TRUE)
            ),
            3
          ),
          collapse = " / "
        ),
        "\n  predicted mu min/median/max =",
        paste(
          round(
            c(
              min(mu, na.rm = TRUE),
              stats::median(mu, na.rm = TRUE),
              max(mu, na.rm = TRUE)
            ),
            6
          ),
          collapse = " / "
        ),
        "\n  psi_AUC =",
        round(psi_auc, 4),
        "\n  capture_AUC =",
        round(capture_auc, 4),
        "\n  lambda_log_score =",
        round(lambda_log_score, 4),
        "\n"
      )

      # ========================================================
      # Store fold result
      # ========================================================

      fold_metrics[[k]] <- data.frame(
        fold = k,
        psi_AUC = psi_auc,
        capture_AUC = capture_auc,
        lambda_log_score = lambda_log_score
      )
    }

    results[[model_name]] <- dplyr::bind_rows(
      fold_metrics
    )
  }

  results
}
