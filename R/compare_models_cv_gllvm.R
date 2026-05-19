#' Cross-validation model comparison for FitModel_gllvm
#'
#' Performs K-fold cross-validation for hierarchical occupancy–
#' capture–abundance models fitted using \code{FitModel_gllvm()}.
#'
#' The function evaluates predictive performance using:
#' \itemize{
#'   \item Occupancy AUC (\code{psi_AUC})
#'   \item Capture AUC (\code{capture_AUC})
#'   \item Predictive abundance log-likelihood
#'         (\code{lambda_log_score})
#' }
#'
#' @param phyloseq A phyloseq object.
#' @param models Named list of model-fitting functions.
#' @param site_col Character. Site column name.
#' @param sample_col Character. Sample column name.
#' @param replicate_col Character. Replicate column name.
#' @param otu_col Character. OTU column name.
#' @param count_col Character. Count column name.
#' @param cv_level Character. One of:
#' \code{"site"}, \code{"sample"},
#' \code{"OTU"}, or \code{"replicate"}.
#' @param K Integer. Number of folds.
#' @param seed Integer random seed.
#'
#' @return A named list containing fold-by-fold
#' cross-validation metrics for each model.
#'
#' 
#'
#' @examples
#' \dontrun{
#'
#' cv_res <- compare_models_cv_gllvm(
#'   phyloseq = ps,
#'   models = models,
#'   site_col = "site_month",
#'   sample_col = "Name",
#'   replicate_col = "Replicate",
#'   otu_col = "OTU",
#'   count_col = "y",
#'   cv_level = "sample",
#'   K = 3
#' )
#'
#' }
#' @export
compare_models_cv_gllvm <- function(
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
  
  # ==========================================================
  # Metadata
  # ==========================================================
  
  meta <- as.data.frame(
    phyloseq::sample_data(phyloseq)
  )
  
  meta$.__sample_id__ <- rownames(meta)
  
  if (!(sample_col %in% names(meta))) {
    meta[[sample_col]] <- meta$.__sample_id__
  }
  
  # ==========================================================
  # Define CV groups
  # ==========================================================
  
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
    
  } else {
    
    groups <- phyloseq::taxa_names(
      phyloseq
    )
  }
  
  groups <- groups[!is.na(groups)]
  
  if (length(groups) < 2) {
    stop("Not enough groups for CV.")
  }
  
  K_eff <- min(K, length(groups))
  
  fold_ids <- sample(
    rep(seq_len(K_eff),
        length.out = length(groups))
  )
  
  folds <- split(groups, fold_ids)
  
  results <- list()
  
  # ==========================================================
  # LOOP OVER MODELS
  # ==========================================================
  
  for (model_name in names(models)) {
    
    cat("\n==============================\n")
    cat("Model:", model_name, "\n")
    cat("CV level:", cv_level, "\n")
    cat("==============================\n")
    
    model_fun <- models[[model_name]]
    
    fold_metrics <- vector(
      "list",
      length(folds)
    )
    
    # ========================================================
    # LOOP OVER FOLDS
    # ========================================================
    
    for (k in seq_along(folds)) {
      
      cat("\n--- Fold", k, "/", length(folds), "---\n")
      
      test_groups <- folds[[k]]
      
      # ------------------------------------------------------
      # Split data
      # ------------------------------------------------------
      
      if (cv_level == "OTU") {
        
        all_otus <- phyloseq::taxa_names(
          phyloseq
        )
        
        test_otus <- intersect(
          test_groups,
          all_otus
        )
        
        train_otus <- setdiff(
          all_otus,
          test_otus
        )
        
        if (length(test_otus) == 0 ||
            length(train_otus) == 0) {
          
          fold_metrics[[k]] <- data.frame(
            psi_AUC = NA_real_,
            capture_AUC = NA_real_,
            lambda_log_score = NA_real_
          )
          
          next
        }
        
        ps_train <- phyloseq::prune_taxa(
          train_otus,
          phyloseq
        )
        
        ps_test <- phyloseq::prune_taxa(
          test_otus,
          phyloseq
        )
        
      } else {
        
        if (cv_level == "site") {
          
          test_samples <- meta$.__sample_id__[
            meta[[site_col]] %in% test_groups
          ]
          
        } else if (cv_level == "sample") {
          
          test_samples <- meta$.__sample_id__[
            meta[[sample_col]] %in% test_groups
          ]
          
        } else {
          
          test_samples <- meta$.__sample_id__[
            meta[[replicate_col]] %in% test_groups
          ]
        }
        
        train_samples <- setdiff(
          meta$.__sample_id__,
          test_samples
        )
        
        if (length(test_samples) == 0 ||
            length(train_samples) == 0) {
          
          fold_metrics[[k]] <- data.frame(
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
      
      # ------------------------------------------------------
      # Fit model
      # ------------------------------------------------------
      
      fit <- tryCatch(
        model_fun(ps_train),
        error = function(e) {
          
          message(
            "Model failed at fold ",
            k,
            ": ",
            e$message
          )
          
          NULL
        }
      )
      
      if (is.null(fit)) {
        
        fold_metrics[[k]] <- data.frame(
          psi_AUC = NA_real_,
          capture_AUC = NA_real_,
          lambda_log_score = NA_real_
        )
        
        next
      }
      
      # ======================================================
      # Prepare test data
      # ======================================================
      
      test_long <- tryCatch(
        
        prepare_long_data(
          physeq_obj = ps_test,
          site_col = site_col,
          nested_cols = unique(
            c(sample_col, replicate_col)
          )
        )$long_df,
        
        error = function(e) NULL
      )
      
      if (is.null(test_long) ||
          nrow(test_long) == 0) {
        
        fold_metrics[[k]] <- data.frame(
          psi_AUC = NA_real_,
          capture_AUC = NA_real_,
          lambda_log_score = NA_real_
        )
        
        next
      }
      
      # ======================================================
      # total_reads
      # ======================================================
      
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
      
      # ======================================================
      # TRUE OCCUPANCY
      # ======================================================
      
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
      
      # ======================================================
      # TRUE CAPTURE
      # ======================================================
      
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
      
      # ======================================================
      # OCCUPANCY PREDICTIONS (GLLVM)
      # ======================================================
      
      psi_pred <- fit$summary |>
        
        dplyr::select(
          dplyr::all_of(site_col),
          dplyr::all_of(otu_col),
          psi_mean
        )
      
      site_eval <- dplyr::left_join(
        site_eval,
        psi_pred,
        by = c(site_col, otu_col)
      )
      
      site_eval$psi <- site_eval$psi_mean
      
      # ======================================================
      # CAPTURE PREDICTIONS
      # ======================================================
      
      capture_pred <- fit$capture |>
        
        dplyr::select(
          dplyr::all_of(site_col),
          dplyr::all_of(sample_col),
          dplyr::all_of(otu_col),
          capture_mean
        )
      
      sample_eval <- dplyr::left_join(
        sample_eval,
        capture_pred,
        by = c(
          site_col,
          sample_col,
          otu_col
        )
      )
      
      sample_eval$capture <- sample_eval$capture_mean
      
      # ======================================================
      # OCCUPANCY AUC
      # ======================================================
      
      site_eval_auc <- site_eval |>
        
        dplyr::filter(
          !is.na(z),
          !is.na(psi)
        )
      
      n_pos_z <- sum(
        site_eval_auc$z == 1
      )
      
      n_neg_z <- sum(
        site_eval_auc$z == 0
      )
      
      psi_auc <- NA_real_
      
      if (n_pos_z > 0 &&
          n_neg_z > 0) {
        
        psi_auc <- tryCatch({
          
          as.numeric(
            pROC::roc(
              response = site_eval_auc$z,
              predictor = site_eval_auc$psi,
              quiet = TRUE
            )$auc
          )
          
        }, error = function(e) {
          
          message(
            "psi_AUC failed: ",
            e$message
          )
          
          NA_real_
        })
        
      } else {
        
        message(
          "psi_AUC skipped: only one class"
        )
      }
      
      # ======================================================
      # CAPTURE AUC
      # ======================================================
      
      sample_eval_auc <- sample_eval |>
        
        dplyr::filter(
          !is.na(a),
          !is.na(capture)
        )
      
      n_pos_a <- sum(
        sample_eval_auc$a == 1
      )
      
      n_neg_a <- sum(
        sample_eval_auc$a == 0
      )
      
      capture_auc <- NA_real_
      
      if (n_pos_a > 0 &&
          n_neg_a > 0) {
        
        capture_auc <- tryCatch({
          
          as.numeric(
            pROC::roc(
              response = sample_eval_auc$a,
              predictor = sample_eval_auc$capture,
              quiet = TRUE
            )$auc
          )
          
        }, error = function(e) {
          
          message(
            "capture_AUC failed: ",
            e$message
          )
          
          NA_real_
        })
        
      } else {
        
        message(
          "capture_AUC skipped: only one class"
        )
      }
      
      # ======================================================
      # ABUNDANCE MODEL
      # ======================================================
      
      abund_fit <- tail(
        fit$abundance_models,
        1
      )[[1]]
      
      # ======================================================
      # Predict lambda
      # ======================================================
      
      test_long$lambda <- suppressWarnings(
        
        as.numeric(
          
          predict(
            abund_fit,
            newdata = test_long,
            type = "response",
            allow.new.levels = TRUE
          )
        )
      )
      
      # ======================================================
      # Predictive log-likelihood
      # ======================================================
      
      y_obs <- test_long[[count_col]]
      
      mu <- pmax(
        test_long$lambda,
        1e-8
      )
      
      family_name <- abund_fit$modelInfo$family$family
      
      # ------------------------------------------------------
      # Poisson
      # ------------------------------------------------------
      
      if (family_name == "poisson") {
        
        loglik_vec <- stats::dpois(
          y_obs,
          lambda = mu,
          log = TRUE
        )
        
      # ------------------------------------------------------
      # Negative Binomial
      # ------------------------------------------------------
        
      } else if (family_name == "nbinom2") {
        
        theta <- sigma(abund_fit)
        
        loglik_vec <- stats::dnbinom(
          y_obs,
          mu = mu,
          size = theta,
          log = TRUE
        )
        
      # ------------------------------------------------------
      # ZIP / ZINB
      # ------------------------------------------------------
        
      } else {
        
        theta <- tryCatch(
          sigma(abund_fit),
          error = function(e) NULL
        )
        
        zi_prob <- tryCatch(
          
          predict(
            abund_fit,
            newdata = test_long,
            type = "zprob",
            allow.new.levels = TRUE
          ),
          
          error = function(e)
            rep(0, nrow(test_long))
        )
        
        if (is.null(theta)) {
          
          # ZIP
          
          loglik_vec <- ifelse(
            
            y_obs == 0,
            
            log(
              zi_prob +
                (1 - zi_prob) *
                exp(-mu)
            ),
            
            log(1 - zi_prob) +
              stats::dpois(
                y_obs,
                lambda = mu,
                log = TRUE
              )
          )
          
        } else {
          
          # ZINB
          
          p0_nb <- stats::dnbinom(
            0,
            mu = mu,
            size = theta
          )
          
          loglik_vec <- ifelse(
            
            y_obs == 0,
            
            log(
              zi_prob +
                (1 - zi_prob) *
                p0_nb
            ),
            
            log(1 - zi_prob) +
              stats::dnbinom(
                y_obs,
                mu = mu,
                size = theta,
                log = TRUE
              )
          )
        }
      }
      
      # ======================================================
      # FINAL LOG SCORE
      # ======================================================
      
      lambda_log_score <- mean(
        loglik_vec,
        na.rm = TRUE
      )
      
      # ======================================================
      # STORE METRICS
      # ======================================================
      
      fold_metrics[[k]] <- data.frame(
        psi_AUC = psi_auc,
        capture_AUC = capture_auc,
        lambda_log_score = lambda_log_score
      )
      
      cat(
        "Fold metrics:",
        "psi_AUC =", round(psi_auc, 3),
        "| capture_AUC =", round(capture_auc, 3),
        "| lambda_log_score =", round(lambda_log_score, 3),
        "\n"
      )
    }
    
    results[[model_name]] <- dplyr::bind_rows(
      fold_metrics
    )
  }
  
  return(results)
}
