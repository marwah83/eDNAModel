#' Fit a Two-Part Occupancy-Abundance Model for Microbiome Data
#'
#' Implements a hierarchical Bayesian-style framework using repeated GLMM fitting
#' to estimate microbial occupancy (presence/absence) and abundance conditional
#' on presence, for OTU-level microbiome data stored in a `phyloseq` object.
#'
#' @param phyloseq A `phyloseq` object containing OTU table, sample data, and taxonomy.
#' @param site_col Character. Name of the column in `sample_data` representing sampling sites.
#' @param abundance_rhs Right-hand side of the formula for abundance model (e.g., `Treatment * OTU`).
#' @param occupancy_rhs Right-hand side of the formula for occupancy model (e.g., `Treatment + OTU`).
#' @param min_species_sum Integer. Minimum total count for an OTU to be retained in analysis. Default: 50.
#' @param abundance_threshold Integer. Minimum count threshold to consider OTU as present. Default: 1.
#' @param n_iter Integer. Number of iterations for Gibbs-like sampling. Default: 50.
#' @param burn_in Integer. Number of iterations to discard as burn-in. Default: 10.
#' @param abundance_family Character. Family for abundance model. One of `"poisson"`, `"nbinom"`, `"zip"`, or `"zinb"`. Default: `"poisson"`.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{Posterior summaries of occupancy (psi), abundance (lambda), and detection probability (p_detect)}
#'   \item{psi_list}{List of occupancy estimates per iteration}
#'   \item{lambda_list}{List of abundance estimates per iteration}
#'   \item{p_detect_list}{List of detection probabilities per iteration}
#'   \item{occupancy_models}{Fitted GLLVM occupancy models}
#'   \item{abundance_models}{Fitted GLMM abundance models}
#'   \item{reduced_data}{Processed input data}
#'   }
#'
#' @details
#' This function decomposes species occurrence into two processes:
#' 1. **Occupancy** (presence/absence): modeled with a binomial GLMM (`glmmTMB`).
#' 2. **Abundance** (counts given presence): modeled with Poisson, Negative Binomial, or Zero-Inflated Poisson/NegBin.
#'
#' At each iteration:
#' - Simulate detection given occupancy and abundance.
#' - Fit occupancy and abundance models.
#' - Update simulated presence (`z_sim`) using current parameter estimates.
#'
#' Posterior summaries are computed using inverse-variance weighting across iterations (after burn-in).
#'
#' @section Model Features:
#' - Supports nested or interaction terms with `OTU` (e.g., `Treatment * OTU`, `Group / OTU`).
#' - Automatically detects treatment or nested variables used in the formula.
#' - Applies reweighting based on prediction uncertainty to generate robust estimates.
#' - Automatically relevels `OTU` to the most abundant taxon for model stability.
#'
#' @examples
#' \dontrun{
#' # Example usage of FitModel with a phyloseq object named `physeq_one`
#'
#' out <- FitModel(
#'  phyloseq = ps_obj,
#' site_col = "Site",
#' abundance_rhs = Treatment * OTU,
#' occupancy_rhs = Treatment + OTU,
#' abundance_family = "nbinom",
#' n_iter = 100,
#' burn_in = 20
#')
#' # Check the output summary
#' head(out$summary)
#'
#' }
#'
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats reformulate rnorm plogis predict
#' @import dplyr
#' @import tidyr
#'
#' @export
FitModel <- function(phyloseq,
                     site_col,
                     abundance_rhs,
                     occupancy_rhs,
                     min_species_sum = 50,
                     abundance_threshold = 1,
                     n_iter = 50,
                     burn_in = 10,
                     abundance_family = "poisson") {

  abundance_rhs <- substitute(abundance_rhs)
  occupancy_rhs <- substitute(occupancy_rhs)

  # Validate family input
  valid_families <- c("poisson", "nbinom", "zip")
  if (!abundance_family %in% valid_families) {
    stop(" 'abundance_family' must be one of: ", paste(valid_families, collapse = ", "))
  }

  # Detect columns from abundance_rhs
  rhs_text <- paste(deparse(abundance_rhs), collapse = " ")
  rhs_text <- gsub("`", "", rhs_text)
  rhs_terms <- unique(all.vars(abundance_rhs))

  # Detect treatment (interaction with OTU)
  treatment_candidates <- rhs_terms[sapply(rhs_terms, function(x)
    grepl(paste0("\\b", x, "\\b\\s*\\*\\s*OTU"), rhs_text) |
      grepl(paste0("OTU\\s*\\*\\s*\\b", x, "\\b"), rhs_text))]

  # Detect nested variables (/ OTU)
  nested_candidates <- rhs_terms[sapply(rhs_terms, function(x)
    grepl(paste0("\\b", x, "\\b\\s*/\\s*OTU"), rhs_text) |
      grepl(paste0("OTU\\s*/\\s*\\b", x, "\\b"), rhs_text))]

  treatment_col <- if (length(treatment_candidates) > 0) treatment_candidates[1] else NULL
  nested_cols   <- if (length(nested_candidates) > 0) nested_candidates else NULL

  message(" Detected columns:")
  if (!is.null(treatment_col)) message("   Treatment: ", treatment_col)
  if (!is.null(nested_cols))   message("   Nested: ", paste(nested_cols, collapse = ", "))
  message(" Using abundance model family: ", abundance_family)

  # Prepare long data
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    min_species_sum = min_species_sum,
    site_col = site_col,
    nested_cols = nested_cols
  )

  long_df <- prep$long_df

  # Set reference OTU
  message(" Calculating most abundant OTU...")
  otu_abundances <- taxa_sums(prep$physeq_filtered)
  top_otu <- names(sort(otu_abundances, decreasing = TRUE))[1]
  long_df$OTU <- factor(long_df$OTU, levels = unique(long_df$OTU))
  long_df$OTU <- relevel(long_df$OTU, ref = top_otu)

  # Initialize storage
  psi_list <- list()
  lambda_list <- list()
  p_detect_list <- list()
  occupancy_models <- list()
  abundance_models <- list()

  if (!is.null(treatment_col)) {
    treatment_levels <- levels(factor(long_df[[treatment_col]]))
    message(" Treatment levels: ", paste(treatment_levels, collapse = ", "))
  }

  # Reduce data: site × OTU
  group_vars <- c(site_col, "OTU")
  reduced_data <- long_df %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::summarise(z_obs = as.integer(sum(y > 0) > abundance_threshold),
                     across(-y, dplyr::first), .groups = "drop") %>%
    dplyr::mutate(z_sim = z_obs)

  if (!is.null(treatment_col)) {
    reduced_data[[treatment_col]] <- factor(reduced_data[[treatment_col]], levels = treatment_levels)
    if (nlevels(reduced_data[[treatment_col]]) < 2) {
      stop("treatment_col must have greater than or equal 2 levels.")
    }
  }

  # Gibbs sampling loop
  for (i in 1:n_iter) {
    message(" Iteration ", i)

    # Occupancy model
    model_occupancy <- glmmTMB::glmmTMB(
      formula = reformulate(deparse(occupancy_rhs), response = "z_sim"),
      data = reduced_data,
      family = binomial
    )

    pred_link <- predict(model_occupancy, type = "link", se.fit = TRUE, newdata = reduced_data)
    logit_draw <- rnorm(nrow(reduced_data), mean = pred_link$fit, sd = pred_link$se.fit)
    psi_pred <- plogis(logit_draw)

    psi_cols <- c(site_col, "OTU")
    if (!is.null(treatment_col)) psi_cols <- c(psi_cols, treatment_col)

    psi_list[[i]] <- setNames(
      data.frame(reduced_data[psi_cols], eta = logit_draw),
      c(psi_cols, "eta")
    )
    occupancy_models[[i]] <- model_occupancy

    # Abundance model
    abundance_data <- long_df %>%
      left_join(reduced_data %>% select(all_of(c(site_col, "OTU", "z_sim"))),
                by = c(site_col, "OTU")) %>%
      filter(z_sim == 1)

    if (!is.null(treatment_col)) {
      abundance_data[[treatment_col]] <- droplevels(factor(abundance_data[[treatment_col]]))
    }

    #abundance_glmm_family <- switch(abundance_family,
     # poisson = poisson(),
      #nbinom  = nbinom2(),
      #zip     = glmmTMB::ziPoisson()
    #)

    # Set abundance GLMM family and zero-inflation formula
if (abundance_family == "poisson") {
  abundance_glmm_family <- poisson
  zi_formula <- ~0
} else if (abundance_family == "nbinom") {
  abundance_glmm_family <- nbinom2
  zi_formula <- ~0
} else if (abundance_family == "zip") {
  abundance_glmm_family <- poisson
  zi_formula <- ~1
} else if (abundance_family == "zinb") {
  abundance_glmm_family <- nbinom2
  zi_formula <- ~1
} else {
  stop(" Invalid abundance_family. Use 'poisson', 'nbinom', 'zip', or 'zinb'.")
}

# Fit abundance model (using your variable name and structure)
model_abundance <- glmmTMB::glmmTMB(
  formula = reformulate(deparse(abundance_rhs), response = "y"),
  data = abundance_data,
  family = abundance_glmm_family,
  ziformula = zi_formula
)


    lambda_pred <- predict(model_abundance, type = "link", se.fit = TRUE, newdata = abundance_data)
    lambda <- exp(lambda_pred$fit)

    lambda_list[[i]] <- setNames(
      data.frame(abundance_data[psi_cols], eta = lambda_pred$fit),
      c(psi_cols, "eta")
    )

    p_detect_list[[i]] <- setNames(
      data.frame(abundance_data[psi_cols], eta = 1 - exp(-lambda)),
      c(psi_cols, "eta")
    )

    abundance_models[[i]] <- model_abundance

    lambda_total <- abundance_data %>%
      mutate(lambda_pred = lambda) %>%
      group_by(across(all_of(c(site_col, "OTU")))) %>%
      summarise(lambda_prod = 1 - prod(1 - exp(-lambda_pred)), .groups = "drop")

    reduced_data <- reduced_data %>%
      select(-any_of("lambda_prod")) %>%
      left_join(lambda_total, by = c(site_col, "OTU")) %>%
      mutate(lambda_prod = tidyr::replace_na(lambda_prod, 0))

    zero_indices <- which(reduced_data$z_obs == 0)
    adjusted_prob <- psi_pred[zero_indices] * reduced_data$lambda_prod[zero_indices]
    reduced_data$z_sim[zero_indices] <- rbinom(length(zero_indices), 1,
                                               pmin(pmax(adjusted_prob, 0.001), 0.999))
  }

  # Posterior summaries
  bind_summary_link <- function(lst, link_type = c("logit", "log", "cloglog")) {
    link_type <- match.arg(link_type)
    df <- bind_rows(lst)

    df_summary <- df %>%
      group_by(across(all_of(psi_cols))) %>%
      summarise(eta_mean = mean(eta), eta_var = var(eta), .groups = "drop")

    df <- df %>%
      left_join(df_summary, by = psi_cols) %>%
      mutate(weight = 1 / ifelse(is.na(eta_var) | eta_var == 0, 1e-6, eta_var))

    df %>%
      group_by(across(all_of(psi_cols))) %>%
      summarise(
        eta_mean = sum(eta * weight) / sum(weight),
        eta_var = sum(weight * (eta - eta_mean)^2) /
          ((sum(weight)^2 - sum(weight^2)) / sum(weight)),
        se_eta = sqrt(eta_var),
        lwr_eta = eta_mean - 1.96 * se_eta,
        upr_eta = eta_mean + 1.96 * se_eta,
        mean = switch(link_type,
                      logit = plogis(eta_mean),
                      log = exp(eta_mean),
                      cloglog = 1 - exp(-exp(eta_mean))),
        lwr = switch(link_type,
                     logit = plogis(lwr_eta),
                     log = exp(lwr_eta),
                     cloglog = 1 - exp(-exp(upr_eta))),
        upr = switch(link_type,
                     logit = plogis(upr_eta),
                     log = exp(upr_eta),
                     cloglog = 1 - exp(-exp(lwr_eta))),
        se = se_eta,
        .groups = "drop"
      ) %>%
      select(all_of(psi_cols), mean, se, lwr, upr)
  }

  psi_summary      <- bind_summary_link(psi_list[-seq_len(burn_in)], "logit")
  lambda_summary   <- bind_summary_link(lambda_list[-seq_len(burn_in)], "log")
  p_detect_summary <- bind_summary_link(p_detect_list[-seq_len(burn_in)], "cloglog")

  final_summary <- psi_summary %>%
    rename_with(~paste0("psi_", .), -all_of(psi_cols)) %>%
    left_join(rename_with(lambda_summary, ~paste0("lambda_", .), -all_of(psi_cols)),
              by = psi_cols) %>%
    left_join(rename_with(p_detect_summary, ~paste0("p_detect_", .), -all_of(psi_cols)),
              by = psi_cols)

  return(list(
    summary = final_summary,
    psi_list = psi_list[-seq_len(burn_in)],
    lambda_list = lambda_list[-seq_len(burn_in)],
    p_detect_list = p_detect_list[-seq_len(burn_in)],
    occupancy_models = occupancy_models,
    abundance_models = abundance_models,
    reduced_data = reduced_data
  ))
}
