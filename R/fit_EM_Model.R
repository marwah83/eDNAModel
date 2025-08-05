library(lme4)

fit_occ_siteZIP_repl_glmm <- function(df,
                                      iters = 400, burn_in = 100,
                                      a.formula = ~ Site + (1|SampleID) + (1|ReplicateID),
                                      o.formula = ~ Site,
                                      z_tol_stop = TRUE, param_tol = 1e-5) {
  stopifnot(all(c("i","s","j","k","x_site","w","y") %in% names(df)))
  
  # --- IDs / factors ---
  df$OTU         <- factor(df$i)
  df$Site        <- factor(df$s)
  df$SampleID    <- interaction(df$s, df$j, drop = TRUE)
  df$ReplicateID <- interaction(df$s, df$j, df$k, drop = TRUE)
  
  # OTU × Site group
  df$g <- interaction(df$OTU, df$Site, drop = TRUE)
  levels_g <- levels(df$g)
  
  # --- Safe init: occupied if any positive in group ---
  any_pos0 <- tapply(df$y > 0, df$g, any)
  any_pos0 <- any_pos0[levels_g]; any_pos0[is.na(any_pos0)] <- FALSE
  if (!any(any_pos0)) stop("All OTU×Site groups have zero counts; cannot initialize.")
  Z_hat_g <- as.integer(any_pos0); names(Z_hat_g) <- levels_g
  
  # Unique OTU×Site frame for occupancy model
  site_frame <- unique(df[c("g","OTU","Site","x_site")])
  rownames(site_frame) <- as.character(site_frame$g)
  
  # --- Helpers ---
  fmt_rhs <- function(f) {
    # turn a one-sided formula into a clean RHS string; fallback to "1" if empty
    s <- paste(deparse(f), collapse = " ")
    s <- sub("^\\s*~\\s*", "", s)
    s <- trimws(s)
    if (nchar(s) == 0) "1" else s
  }
  get_fixef <- function(fit) {
    cf <- summary(fit)$coef
    list(est = cf[, "Estimate"], se = cf[, "Std. Error"])
  }
  
  # Build column layout lazily on first successful fits
  trace_est <- trace_se <- NULL
  
  # Fallback lambda to allow first E-step if Poisson fit fails
  lambda_hat <- rep(0.1, nrow(df))
  
  for (t in seq_len(iters)) {
    # ------- M-step: Poisson GLMM on rows with Z_hat=1 -------
    df$Z_hat <- as.integer(Z_hat_g[as.character(df$g)] == 1L)
    Q <- df[df$Z_hat == 1L, , drop = FALSE]
    
    fit_pois <- NULL
    fe_p <- list(est = numeric(0), se = numeric(0))
    if (nrow(Q) > 0) {
      rhs_a <- fmt_rhs(a.formula)             # e.g., "Site + (1|SampleID) + (1|ReplicateID)"
      rhs_full <- if (identical(rhs_a, "1")) "w" else paste("w +", rhs_a)
      f_pois <- as.formula(paste("y ~", rhs_full))
      fit_pois <- tryCatch(
        glmer(f_pois, data = Q, family = poisson(link = "log"),
              control = glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)),
        error = function(e) NULL
      )
      if (!is.null(fit_pois)) {
        fe_p <- get_fixef(fit_pois)
        # Marginal mean (exclude RE) for E-step
        lambda_hat <- pmax(predict(fit_pois, newdata = df, type = "response", re.form = NA), 1e-12)
      }
    }
    
    # ------- M-step: Occupancy (GLM with fixed Site) on OTU×Site -------
    site_frame$Z_hat <- Z_hat_g[as.character(site_frame$g)]
    rhs_o <- fmt_rhs(o.formula)               # e.g., "Site" (or "1")
    f_logit <- as.formula(paste("Z_hat ~", rhs_o))
    fit_logit <- tryCatch(
      glm(f_logit, data = site_frame, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    
    fe_b <- list(est = numeric(0), se = numeric(0))
    if (!is.null(fit_logit)) {
      fe_b <- get_fixef(fit_logit)
      psi_hat_g <- pmin(pmax(predict(fit_logit, newdata = site_frame, type = "response"), 1e-12), 1 - 1e-12)
      names(psi_hat_g) <- rownames(site_frame)
    } else {
      psi_hat_g <- setNames(rep(mean(Z_hat_g), nrow(site_frame)), rownames(site_frame))
    }
    
    # ------- Collect fixed effects into stable columns -------
    pois_fix_names <- names(fe_p$est)
    bin_fix_names  <- names(fe_b$est)
    
    if (is.null(trace_est)) {
      pois_cols <- c("(Intercept)", "w", grep("^Site", pois_fix_names, value = TRUE))
      bin_cols  <- c("(Intercept)", grep("^Site", bin_fix_names, value = TRUE))
      all_cols  <- c(paste0("Pois:", pois_cols), paste0("Occ:", bin_cols))
      trace_est <- matrix(NA_real_, nrow = iters, ncol = length(all_cols),
                          dimnames = list(NULL, all_cols))
      trace_se  <- matrix(NA_real_, nrow = iters, ncol = length(all_cols),
                          dimnames = list(NULL, all_cols))
    } else {
      pois_cols <- sub("^Pois:", "", grep("^Pois:", colnames(trace_est), value = TRUE))
      bin_cols  <- sub("^Occ:",  "", grep("^Occ:",  colnames(trace_est), value = TRUE))
    }
    
    # Fill Poisson FE
    for (nm in pois_cols) {
      colnm <- paste0("Pois:", nm)
      trace_est[t, colnm] <- if (nm %in% pois_fix_names) fe_p$est[nm] else NA_real_
      trace_se [t, colnm] <- if (nm %in% pois_fix_names) fe_p$se [nm] else NA_real_
    }
    # Fill Occupancy FE
    for (nm in bin_cols) {
      colnm <- paste0("Occ:", nm)
      trace_est[t, colnm] <- if (nm %in% bin_fix_names) fe_b$est[nm] else NA_real_
      trace_se [t, colnm] <- if (nm %in% bin_fix_names) fe_b$se [nm] else NA_real_
    }
    
    # ------- E-step: update Z_hat_g -------
    any_pos <- tapply(df$y > 0, df$g, any)
    any_pos <- any_pos[levels_g]; any_pos[is.na(any_pos)] <- FALSE
    Z_new <- as.integer(any_pos); names(Z_new) <- levels_g
    
    zero_groups <- levels_g[!any_pos]
    if (length(zero_groups)) {
      idx_zero <- df$g %in% zero_groups
      Lambda_sum <- tapply(lambda_hat[idx_zero], droplevels(df$g[idx_zero]), sum)
      Lambda_sum <- Lambda_sum[zero_groups]; Lambda_sum[is.na(Lambda_sum)] <- 0
      psi_vec <- psi_hat_g[zero_groups]; psi_vec[is.na(psi_vec)] <- mean(Z_hat_g)
      q  <- exp(-Lambda_sum)
      p1 <- (psi_vec * q) / (psi_vec * q + (1 - psi_vec))
      p1 <- pmin(pmax(p1, 1e-8), 1 - 1e-8)
      # deterministic update (or use Bernoulli draws for stochastic EM)
      Z_new[zero_groups] <- as.integer(p1 > 0.5)
    }
    
    # Convergence check
    z_diff <- sum(Z_new != Z_hat_g)
    param_diff <- if (t > 1) max(abs(trace_est[t, ] - trace_est[t-1, ]), na.rm = TRUE) else Inf
    Z_hat_g <- Z_new
    if (z_tol_stop && z_diff == 0 && is.finite(param_diff) && param_diff < param_tol) {
      message(sprintf("Converged at iter %d: Z unchanged, max |Δparam| = %.3g", t, param_diff))
      trace_est <- trace_est[1:t, , drop = FALSE]
      trace_se  <- trace_se [1:t, , drop = FALSE]
      break
    }
  }
  
  # Burn-in
  keep_from <- min(burn_in + 1L, nrow(trace_est))
  post_est  <- trace_est[keep_from:nrow(trace_est), , drop = FALSE]
  post_se   <- trace_se [keep_from:nrow(trace_se ), , drop = FALSE]
  
  list(
    trace_est = trace_est,
    trace_se  = trace_se,
    estimates = post_est,
    ses       = post_se,
    fe_means  = colMeans(post_est, na.rm = TRUE),
    last_pois_re_sd = if (exists("fit_pois") && !is.null(fit_pois)) lme4::VarCorr(fit_pois) else NULL
  )
}

sim <- simulate_occ_siteZIP_repl(I = 60, S = 12, J = 10, K = 4,
                                 beta_0 = 0.5, beta_1 = -0.05,
                                 gamma_0 = 0.5, gamma_1 = 0.3,
                                 seed = 42)

fit <- fit_occ_siteZIP_repl_glmm(
  sim$data_long,
  iters = 600,
  burn_in = 200,
  a.formula = ~ Site + (1|SampleID) + (1|ReplicateID),
  o.formula = ~ Site
)
