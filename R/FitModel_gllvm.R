
#' Fit a Multispecies Three-Level Occupancy–Detection–Abundance Model
#'
#' @description
#' Fits an approximate hierarchical occupancy–detection–abundance model
#' for multispecies eDNA or microbiome count data stored in a
#' \code{phyloseq} object.
#'
#' The model distinguishes three hierarchical levels:
#' \itemize{
#'   \item \strong{Site level (Z):} latent taxon occupancy.
#'   \item \strong{Biological sample level (A):} latent taxon capture.
#'   \item \strong{PCR replicate level (Y):} observed sequence counts.
#' }
#'
#' Site-level occupancy is modeled using a binomial generalized linear
#' latent variable model (GLLVM). Capture and abundance are modeled
#' using generalized linear mixed models fitted with \code{glmmTMB}.
#'
#' Latent occupancy and capture states are updated iteratively using
#' a stochastic data-augmentation procedure.
#'
#' @param phyloseq A \code{phyloseq} object containing taxon counts
#'   and sample metadata.
#' @param site_col Character string identifying the site column.
#' @param sample_col Character string identifying biological samples.
#'   Defaults to \code{"Name"}.
#' @param replicate_col Character string identifying PCR replicates,
#'   or \code{NULL} when no explicit replicate identifier is used.
#' @param otu_col Character string identifying the taxon column.
#'   Defaults to \code{"OTU"}.
#' @param count_col Character string identifying the observed count
#'   variable. Defaults to \code{"y"}.
#' @param abundance_rhs Right-hand side of the abundance model, or a
#'   complete formula with response matching \code{count_col}.
#'   For example, \code{(1 | OTU)} or
#'   \code{y ~ offset(log(total_reads)) + (1 | OTU)}.
#' @param capture_formula Formula for the capture model, with response
#'   \code{a_sim}. By default, a binomial random-intercept model
#'   with a taxon-specific random effect is used.
#' @param occupancy_covars Optional character vector of site-level
#'   occupancy covariates. Covariate values must be constant within
#'   each site.
#' @param min_species_sum Minimum total count required to retain
#'   a taxon. Defaults to 50.
#' @param min_detection_replicates Minimum number of observations
#'   exceeding \code{abundance_threshold} required to retain a taxon.
#'   Defaults to 1.
#' @param abundance_threshold Count threshold used to define an
#'   observed detection. Detection occurs when the count is strictly
#'   greater than this threshold. Defaults to 1.
#' @param n_iter Number of stochastic model-fitting iterations.
#'   Defaults to 50.
#' @param burn_in Number of initial iterations excluded from
#'   summary calculations. Must be smaller than \code{n_iter}.
#'   Defaults to 10.
#' @param abundance_family Abundance distribution. One of
#'   \code{"poisson"}, \code{"nbinom"}, \code{"zip"}, or
#'   \code{"zinb"}. Defaults to \code{"poisson"}.
#' @param num_lv_c Number of constrained latent variables.
#'   Currently must equal 0 because constrained latent-variable
#'   estimation is not implemented in this wrapper.
#' @param num_lv Number of unconstrained latent variables in the
#'   GLLVM occupancy model. Defaults to 1. Larger values increase
#'   model complexity and may lead to numerical instability.
#' @param gllvm_n_init Number of GLLVM initialization attempts.
#'   Defaults to 10.
#' @param gllvm_gradient_tol Positive numerical tolerance for the
#'   maximum absolute gradient. Defaults to \code{1e-3}.
#' @param gllvm_strict Logical. If \code{TRUE}, skip iterations
#'   failing the implemented numerical acceptance criteria.
#'   If \code{FALSE}, record diagnostic failures but allow
#'   subsequent stochastic updates. Defaults to \code{FALSE}.
#' @param gllvm_condition_tol Hessian condition-number threshold
#'   used to identify ill-conditioned fits. Defaults to \code{1e8}.
#' @param gllvm_eigen_tol Positive relative eigenvalue threshold
#'   used to classify near-zero or negative eigenvalues.
#'   Defaults to \code{1e-8}.
#' @param gllvm_warn_separation Logical. If \code{TRUE}, warn when
#'   an OTU is universally present or absent in the current
#'   augmented occupancy matrix. Defaults to \code{FALSE}.
#' @param gllvm_auto_exclude Logical. If \code{TRUE}, exclude OTUs
#'   with insufficient variation in their augmented occupancy
#'   states from the GLLVM fit. Excluded OTUs remain in the
#'   downstream capture and abundance components.
#'   Defaults to \code{TRUE}.
#' @param gllvm_min_occupied_sites Non-negative integer specifying
#'   the occupancy screening threshold. An OTU is excluded from
#'   the GLLVM when its number of occupied sites or unoccupied
#'   sites is less than or equal to this threshold.
#'   A value of 0 excludes constant occupancy columns.
#'   Defaults to 0.
#' @param gllvm_excluded_prior Numeric vector of length two
#'   containing positive Beta smoothing parameters \eqn{(a,b)}
#'   for OTUs excluded from the GLLVM. Defaults to
#'   \code{c(1, 1)}.
#' @param verbose Logical. If \code{TRUE}, print progress and
#'   diagnostic information. Defaults to \code{TRUE}.
#'
#' @details
#' \strong{Hierarchical model}
#'
#' Let \eqn{i} index sites, \eqn{j} biological samples,
#' \eqn{k} PCR replicates, and \eqn{m} taxa.
#'
#' The latent site occupancy state is
#' \deqn{Z_{im} \sim \mathrm{Bernoulli}(\psi_{im}).}
#'
#' Conditional on site occupancy, the biological sample capture
#' state is
#' \deqn{
#' A_{ijm} \mid Z_{im}
#' \sim \mathrm{Bernoulli}(Z_{im}p_{ijm}).
#' }
#'
#' Conditional on capture, observed counts follow the selected
#' abundance distribution. In the Poisson case,
#' \deqn{
#' Y_{ijkm} \mid A_{ijm}=1
#' \sim \mathrm{Poisson}(\lambda_{ijkm}).
#' }
#'
#' Structural zeros arise when \eqn{A_{ijm}=0}.
#'
#' \strong{Observed detection}
#'
#' The observed site-level detection indicator is
#' \deqn{
#' z^{obs}_{im}
#' = I\left(\max_{j,k}Y_{ijkm}>c\right),
#' }
#' where \eqn{c} is \code{abundance_threshold}.
#'
#' \strong{Occupancy model}
#'
#' The GLLVM models taxon-specific occupancy probabilities using
#' site-level covariates and latent variables:
#' \deqn{
#' \mathrm{logit}(\psi_{im})
#' = \alpha_m + X_i^\top\beta_m
#' + u_i^\top\theta_m.
#' }
#'
#' Here, \eqn{u_i} denotes site latent variables and
#' \eqn{\theta_m} denotes taxon-specific loadings.
#' When no occupancy covariates are supplied, the model uses
#' taxon intercepts and latent variables.
#'
#' \strong{Automatic OTU screening}
#'
#' At each iteration, the augmented occupancy matrix is examined
#' for taxa with insufficient variation.
#'
#' When \code{gllvm_auto_exclude = TRUE}, such taxa are excluded
#' from the GLLVM fitting matrix, but are not removed from the
#' hierarchical dataset.
#'
#' Their occupancy probabilities are approximated using
#' Beta-Binomial smoothing:
#' \deqn{
#' \widehat{\psi}_m
#' = \frac{n_{occupied,m}+a}{n_{sites,m}+a+b}.
#' }
#'
#' These probabilities are conditional approximations, not
#' joint GLLVM estimates.
#'
#' \strong{Capture model}
#'
#' Capture probabilities are modeled using a binomial GLMM:
#' \deqn{
#' \mathrm{logit}(p_{ijm})
#' = W_{ijm}^{\top}\gamma + b_m,
#' }
#' where the random-effects structure is determined by
#' \code{capture_formula}.
#'
#' \strong{Abundance model}
#'
#' Conditional on capture, the abundance model is fitted using
#' \code{glmmTMB}, supporting Poisson, negative-binomial,
#' zero-inflated Poisson, and zero-inflated negative-binomial
#' distributions.
#'
#' Optional offsets, such as
#' \code{offset(log(total_reads))}, can account for differences
#' in sequencing depth.
#'
#' \strong{Detection probabilities}
#'
#' For an abundance distribution with zero-count probability
#' \eqn{P_0}, the probability of at least one read is
#' \deqn{p_{detect}=1-P_0.}
#'
#' For a Poisson distribution,
#' \deqn{P_0=\exp(-\lambda).}
#'
#' For a negative-binomial distribution with mean \eqn{\lambda}
#' and size parameter \eqn{\theta},
#' \deqn{
#' P_0=\left(\frac{\theta}{\theta+\lambda}\right)^\theta.
#' }
#'
#' For zero-inflated distributions with structural-zero
#' probability \eqn{\pi},
#' \deqn{P_0=\pi+(1-\pi)P_{0,count}.}
#'
#' When zero probabilities are combined over PCR replicates,
#' the sample-level probability of no reads is calculated
#' from the product of replicate-level zero probabilities,
#' conditional on the fitted model.
#'
#' The complementary log-log transformation may be used for
#' the probability of at least one detection:
#' \deqn{\eta=\log(-\log(P_0)).}
#'
#' Detection probabilities are derived from fitted abundance
#' parameters; they are not fitted as an independent GLMM.
#'
#' @section Stochastic estimation algorithm:
#' The function uses an iterative, EM-like data-augmentation
#' procedure rather than full joint Bayesian inference.
#'
#' \enumerate{
#'   \item Construct the current augmented site-by-taxon
#'     occupancy matrix.
#'   \item Screen OTUs and fit the binomial GLLVM.
#'   \item Obtain occupancy probabilities for all retained OTUs.
#'   \item Fit the capture GLMM conditional on occupied sites.
#'   \item Fit the abundance GLMM conditional on captured samples.
#'   \item Derive zero-count and detection probabilities.
#'   \item Update latent capture and occupancy states using
#'     stochastic conditional sampling.
#'   \item Repeat and summarize retained iterations after burn-in.
#' }
#'
#' @section Numerical convergence diagnostics:
#' Each GLLVM fit is assessed using optimizer status,
#' maximum absolute gradient, Hessian eigenvalues, covariance
#' eigenvalues, and Hessian condition number.
#'
#' The returned \code{gllvm_diagnostics} data frame includes:
#' \itemize{
#'   \item \code{optimizer_converged}: optimizer convergence flag.
#'   \item \code{max_abs_gradient}: largest absolute gradient.
#'   \item \code{hessian_min_eigen}: minimum Hessian eigenvalue.
#'   \item \code{hessian_negative_eigen}: number of materially
#'     negative Hessian eigenvalues.
#'   \item \code{hessian_near_zero_eigen}: number of near-zero
#'     Hessian eigenvalues.
#'   \item \code{covariance_min_eigen}: minimum covariance eigenvalue.
#'   \item \code{covariance_negative_eigen}: number of materially
#'     negative covariance eigenvalues.
#'   \item \code{covariance_near_zero_eigen}: number of near-zero
#'     covariance eigenvalues.
#'   \item \code{hessian_condition_number}: estimated Hessian
#'     condition number.
#'   \item \code{hessian_ill_conditioned}: indicator that the
#'     condition number exceeds \code{gllvm_condition_tol}.
#'   \item \code{hessian_positive_definite}: whether the
#'     Hessian has strictly positive eigenvalues.
#'   \item \code{numerical_ok}: overall numerical acceptance flag.
#' }
#'
#' Optimizer convergence alone does not guarantee a well-conditioned
#' Hessian or reliable standard errors.
#'
#' For example, increasing the number of latent variables may
#' cause weak identifiability in datasets with few sites or taxa.
#'
#' @section Iteration retention:
#' The returned \code{iteration_success} vector identifies iterations
#' that completed the required estimation steps.
#'
#' The \code{retained_iterations} vector identifies successful
#' iterations retained after burn-in.
#'
#' When \code{gllvm_strict = FALSE}, an iteration may be retained
#' even when \code{numerical_ok = FALSE}.
#'
#' Users should inspect both computational success and numerical
#' diagnostics before interpreting the fitted model.
#'
#' @section Interpretation and limitations:
#' \itemize{
#'   \item The procedure provides approximate hierarchical inference
#'     rather than a full joint likelihood or posterior distribution.
#'   \item Iteration-based intervals should not automatically be
#'     interpreted as calibrated Bayesian credible intervals.
#'   \item Occupancy estimates for automatically excluded OTUs
#'     are based on Beta-Binomial smoothing.
#'   \item Detection probabilities depend on the fitted abundance
#'     distribution and its parameters.
#'   \item Numerical instability may arise from sparse taxa,
#'     separation, overparameterization, or weak identifiability.
#'   \item Latent-variable signs and orientations may differ across
#'     iterations; averaging latent coordinates requires appropriate
#'     alignment.
#'   \item Large abundance means can imply near-certain detection.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{Site-by-taxon summaries for occupancy
#'     (\code{psi_*}), capture (\code{capture_*}), abundance
#'     (\code{lambda_*}), and detection (\code{p_detect_*}).}
#'   \item{capture}{Biological sample-level capture summaries.}
#'   \item{capture_site}{Site-level aggregated capture summaries.}
#'   \item{psi_list}{Per-iteration occupancy linear predictors.}
#'   \item{capture_list}{Per-iteration capture linear predictors.}
#'   \item{lambda_list}{Per-iteration abundance linear predictors.}
#'   \item{p_detect_list}{Per-iteration detection linear predictors.}
#'   \item{occupancy_models}{Fitted GLLVM occupancy models.}
#'   \item{capture_models}{Fitted capture GLMMs.}
#'   \item{abundance_models}{Fitted abundance GLMMs.}
#'   \item{reduced_data}{Site-by-taxon augmented occupancy data.}
#'   \item{sample_data}{Biological sample-by-taxon capture data.}
#'   \item{long_df}{Processed long-format count data.}
#'   \item{lv_sites}{Site latent-variable estimates.}
#'   \item{lv_species}{Taxon latent-variable loadings.}
#'   \item{mean_lv_sites}{Mean site latent-variable coordinates.}
#'   \item{mean_lv_species}{Mean taxon latent-variable loadings.}
#'   \item{filter_summary}{Taxon filtering information.}
#'   \item{diagnostic_AIC}{Per-iteration model fit statistics.}
#'   \item{gllvm_diagnostics}{Per-iteration GLLVM numerical diagnostics.}
#'   \item{gllvm_separation}{OTU separation diagnostics.}
#'   \item{gllvm_excluded_otus}{OTUs excluded from the GLLVM
#'     component at each iteration.}
#'   \item{occupancy_model_iterations}{Iteration identifiers
#'     corresponding to fitted occupancy models.}
#'   \item{gllvm_settings}{GLLVM fitting and diagnostic settings.}
#'   \item{iteration_success}{Logical vector identifying
#'     computationally successful iterations.}
#'   \item{retained_iterations}{Iteration identifiers retained
#'     after burn-in.}
#'   \item{abundance_formula}{Formula used for the abundance model.}
#'   \item{capture_formula}{Formula used for the capture model.}
#'   \item{note}{Additional interpretation notes.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- FitModel_gllvm(
#'   phyloseq = ps,
#'   site_col = "Site",
#'   sample_col = "Name",
#'   replicate_col = NULL,
#'   otu_col = "OTU",
#'   count_col = "y",
#'   abundance_rhs = (1 | OTU),
#'   capture_formula = a_sim ~ 1 + (1 | OTU),
#'   abundance_family = "poisson",
#'   n_iter = 10,
#'   burn_in = 2,
#'   num_lv = 1,
#'   gllvm_n_init = 10,
#'   gllvm_gradient_tol = 1e-3,
#'   gllvm_strict = FALSE,
#'   gllvm_auto_exclude = TRUE,
#'   gllvm_min_occupied_sites = 1
#' )
#'
#' head(fit$summary)
#' fit$gllvm_diagnostics
#' fit$retained_iterations
#' }
#'
#' @import dplyr
#' @import tidyr
#' @importFrom gllvm gllvm
#' @importFrom glmmTMB glmmTMB nbinom2
#' @importFrom reshape2 acast
#' @export
FitModel_gllvm <- function(
    phyloseq,
    site_col,
    sample_col = "Name",
    replicate_col = NULL,
    otu_col = "OTU",
    count_col = "y",
    abundance_rhs,
    capture_formula = NULL,
    occupancy_covars = NULL,
    min_species_sum = 50,
    min_detection_replicates = 1,
    abundance_threshold = 1,
    n_iter = 50,
    burn_in = 10,
    abundance_family = c("poisson", "nbinom", "zip", "zinb"),
    num_lv_c = 0,
    num_lv = 1,
    gllvm_n_init = 10,
    gllvm_gradient_tol = 1e-3,
    gllvm_strict = FALSE,
    gllvm_condition_tol = 1e8,
    gllvm_eigen_tol = 1e-8,
    gllvm_warn_separation = FALSE,
    gllvm_auto_exclude = TRUE,
    gllvm_min_occupied_sites = 0,
    gllvm_excluded_prior = c(1, 1),
    verbose = TRUE
) {
  
  # ============================================================
  # 0. Basic argument checks
  # ============================================================
  
  abundance_family <- match.arg(abundance_family)
  
  if (!is.numeric(n_iter) || length(n_iter) != 1L ||
      n_iter < 1 || n_iter != as.integer(n_iter)) {
    stop("n_iter must be a positive integer.")
  }
  
  if (!is.numeric(burn_in) || length(burn_in) != 1L ||
      burn_in < 0 || burn_in != as.integer(burn_in)) {
    stop("burn_in must be a non-negative integer.")
  }
  
  if (burn_in >= n_iter) {
    stop("burn_in must be smaller than n_iter.")
  }
  
  if (!is.numeric(num_lv_c) || length(num_lv_c) != 1L ||
      num_lv_c < 0 || num_lv_c != as.integer(num_lv_c)) {
    stop("num_lv_c must be a non-negative integer.")
  }
  
  if (!is.numeric(num_lv) || length(num_lv) != 1L ||
      is.na(num_lv) || !is.finite(num_lv) ||
      num_lv < 0 || num_lv != as.integer(num_lv)) {
    stop("num_lv must be a non-negative integer.")
  }
  
  # This implementation fits an unconstrained GLLVM with optional
  # site-level covariates X_cov. Constrained latent variables require
  # a separate specification and are not implemented here.
  if (num_lv_c != 0L) {
    stop("num_lv_c must be 0 in this implementation. Use num_lv = 1 or 2 for unconstrained latent variables.")
  }
  
  if (length(gllvm_n_init) != 1L || !is.finite(gllvm_n_init) ||
      gllvm_n_init < 1 || gllvm_n_init != as.integer(gllvm_n_init)) {
    stop("gllvm_n_init must be a positive integer.")
  }
  if (length(gllvm_gradient_tol) != 1L || !is.finite(gllvm_gradient_tol) ||
      gllvm_gradient_tol <= 0) {
    stop("gllvm_gradient_tol must be positive and finite.")
  }
  if (!is.logical(gllvm_strict) || length(gllvm_strict) != 1L ||
      is.na(gllvm_strict)) {
    stop("gllvm_strict must be TRUE or FALSE.")
  }
  
  if (!is.numeric(gllvm_condition_tol) || length(gllvm_condition_tol) != 1L ||
      is.na(gllvm_condition_tol) || gllvm_condition_tol <= 1) {
    stop("gllvm_condition_tol must be numeric and greater than 1.")
  }
  if (!is.numeric(gllvm_eigen_tol) || length(gllvm_eigen_tol) != 1L ||
      !is.finite(gllvm_eigen_tol) || gllvm_eigen_tol <= 0) {
    stop("gllvm_eigen_tol must be positive and finite.")
  }
  if (!is.logical(gllvm_warn_separation) || length(gllvm_warn_separation) != 1L ||
      is.na(gllvm_warn_separation)) {
    stop("gllvm_warn_separation must be TRUE or FALSE.")
  }
  
  if (!is.logical(gllvm_auto_exclude) || length(gllvm_auto_exclude) != 1L ||
      is.na(gllvm_auto_exclude)) stop("gllvm_auto_exclude must be TRUE or FALSE.")
  if (length(gllvm_min_occupied_sites) != 1L ||
      !is.finite(gllvm_min_occupied_sites) ||
      gllvm_min_occupied_sites < 0 ||
      gllvm_min_occupied_sites != as.integer(gllvm_min_occupied_sites)) {
    stop("gllvm_min_occupied_sites must be a nonnegative integer.")
  }
  if (!is.numeric(gllvm_excluded_prior) ||
      length(gllvm_excluded_prior) != 2L ||
      any(!is.finite(gllvm_excluded_prior)) ||
      any(gllvm_excluded_prior <= 0)) {
    stop("gllvm_excluded_prior must contain two positive Beta prior parameters.")
  }
  
  # ============================================================
  # Helper functions
  # ============================================================
  
  bt <- function(x) {
    paste0("`", x, "`")
  }
  
  # Force structural identifiers to factors
  to_factor_cols <- function(df, cols) {
    
    cols <- unique(stats::na.omit(cols))
    cols <- intersect(cols, names(df))
    
    for (col in cols) {
      df[[col]] <- as.factor(as.character(df[[col]]))
    }
    
    df
  }
  
  # Convert character predictors to factors, while preserving
  # genuinely numeric covariates as numeric
  character_to_factor <- function(df, cols) {
    
    cols <- unique(stats::na.omit(cols))
    cols <- intersect(cols, names(df))
    
    for (col in cols) {
      
      if (is.character(df[[col]])) {
        df[[col]] <- factor(df[[col]])
      }
    }
    
    df
  }
  
  get_formula_vars <- function(formula, response) {
    setdiff(all.vars(formula), response)
  }
  
  # ============================================================
  # 1. SAFE abundance formula handling
  #
  # Supports BOTH:
  #
  # abundance_rhs =
  #   (1 | OTU) + (1 | Samplingmonth / OTU)
  #
  # AND:
  #
  # abundance_rhs =
  #   y ~ (1 | OTU) + (1 | Samplingmonth / OTU)
  #
  # IMPORTANT:
  # substitute() MUST be called before abundance_rhs is evaluated.
  # ============================================================
  
  abundance_rhs_expr <- substitute(abundance_rhs)
  
  abundance_rhs_txt <- paste(
    deparse(abundance_rhs_expr),
    collapse = " "
  )
  
  is_full_formula <-
    is.call(abundance_rhs_expr) &&
    identical(abundance_rhs_expr[[1L]], as.name("~"))
  
  if (is_full_formula) {
    
    abundance_formula <- stats::as.formula(
      abundance_rhs_txt,
      env = parent.frame()
    )
    
  } else {
    
    abundance_formula <- stats::as.formula(
      paste(
        bt(count_col),
        "~",
        abundance_rhs_txt
      ),
      env = parent.frame()
    )
  }
  
  # ============================================================
  # 2. Default capture formula
  # ============================================================
  
  if (is.null(capture_formula)) {
    
    capture_formula <- stats::as.formula(
      paste(
        "a_sim ~ 1 + (1 |",
        bt(otu_col),
        ")"
      ),
      env = parent.frame()
    )
  }
  
  if (!inherits(capture_formula, "formula")) {
    stop("capture_formula must be a formula.")
  }
  
  # ============================================================
  # 3. Prepare long data
  # ============================================================
  
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col,
    nested_cols = unique(
      stats::na.omit(
        c(sample_col, replicate_col)
      )
    )
  )
  
  long_df <- prep$long_df
  
  # ============================================================
  # 4. Required-column checks
  # ============================================================
  
  needed_cols <- unique(
    c(
      site_col,
      sample_col,
      otu_col,
      count_col
    )
  )
  
  missing_cols <- setdiff(
    needed_cols,
    names(long_df)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in long_df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.null(replicate_col) &&
      !(replicate_col %in% names(long_df))) {
    
    stop(
      "replicate_col '",
      replicate_col,
      "' not found in long_df."
    )
  }
  
  if (!is.null(occupancy_covars) &&
      length(occupancy_covars) > 0) {
    
    missing_occ <- setdiff(
      occupancy_covars,
      names(long_df)
    )
    
    if (length(missing_occ) > 0) {
      stop(
        "Missing occupancy covariates in long_df: ",
        paste(missing_occ, collapse = ", ")
      )
    }
  }
  
  # ============================================================
  # 5. Type handling
  # ============================================================
  
  long_df[[site_col]]   <- as.character(long_df[[site_col]])
  long_df[[sample_col]] <- as.character(long_df[[sample_col]])
  long_df[[otu_col]]    <- as.character(long_df[[otu_col]])
  long_df[[count_col]]  <- as.numeric(long_df[[count_col]])
  
  if (!is.null(replicate_col)) {
    long_df[[replicate_col]] <-
      as.character(long_df[[replicate_col]])
  }
  
  long_df <- long_df[
    !is.na(long_df[[count_col]]),
    ,
    drop = FALSE
  ]
  
  if (nrow(long_df) == 0) {
    stop("No non-missing count observations remain.")
  }
  
  # ============================================================
  # 6. Compute sample-level sequencing depth
  #
  # Use site + sample rather than sample alone so the function
  # remains valid even when sample labels repeat among sites.
  # ============================================================
  
  sample_id_keys <- unique(
    c(site_col, sample_col)
  )
  
  total_reads_df <- long_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(sample_id_keys)
      )
    ) |>
    dplyr::summarise(
      total_reads = sum(
        .data[[count_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  long_df <- dplyr::left_join(
    long_df,
    total_reads_df,
    by = sample_id_keys
  )
  
  # ============================================================
  # 7. Offset handling
  # ============================================================
  
  uses_offset <- any(
    grepl(
      "offset\\s*\\(",
      deparse(abundance_formula)
    )
  )
  
  if (uses_offset) {
    
    zero_samples <- total_reads_df |>
      dplyr::filter(
        is.na(.data$total_reads) |
          .data$total_reads <= 0
      )
    
    if (nrow(zero_samples) > 0) {
      
      if (verbose) {
        message(
          "Offset detected -> removing ",
          nrow(zero_samples),
          " zero-read sample(s)."
        )
      }
      
      zero_samples$.drop_sample <- TRUE
      
      long_df <- long_df |>
        dplyr::left_join(
          zero_samples |>
            dplyr::select(
              dplyr::all_of(sample_id_keys),
              .data$.drop_sample
            ),
          by = sample_id_keys
        ) |>
        dplyr::filter(
          is.na(.data$.drop_sample)
        ) |>
        dplyr::select(
          -dplyr::any_of(".drop_sample")
        )
    }
    
    if (any(
      long_df$total_reads <= 0,
      na.rm = TRUE
    )) {
      
      stop(
        "total_reads must be positive when using ",
        "offset(log(total_reads))."
      )
    }
    
  } else {
    
    if (verbose) {
      message("No abundance offset detected.")
    }
  }
  
  # ============================================================
  # 8. OTU filtering
  # ============================================================
  
  otu_stats <- long_df |>
    dplyr::group_by(
      .data[[otu_col]]
    ) |>
    dplyr::summarise(
      total_count = sum(
        .data[[count_col]],
        na.rm = TRUE
      ),
      detected_replicates = sum(
        .data[[count_col]] > abundance_threshold,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  keep_otus <- otu_stats |>
    dplyr::filter(
      .data$total_count >= min_species_sum,
      .data$detected_replicates >=
        min_detection_replicates
    ) |>
    dplyr::pull(
      .data[[otu_col]]
    )
  
  long_df <- long_df |>
    dplyr::filter(
      .data[[otu_col]] %in% keep_otus
    )
  
  if (nrow(long_df) == 0) {
    stop("No OTUs remain after filtering.")
  }
  
  if (length(keep_otus) < 2) {
    stop(
      "At least two OTUs are required for the GLLVM occupancy model."
    )
  }
  
  if (verbose) {
    
    message(
      "OTUs before filtering: ",
      dplyr::n_distinct(
        otu_stats[[otu_col]]
      )
    )
    
    message(
      "OTUs after filtering: ",
      dplyr::n_distinct(
        long_df[[otu_col]]
      )
    )
  }
  
  # ============================================================
  # 9. Factor / covariate handling
  # ============================================================
  
  id_factor_cols <- unique(
    stats::na.omit(
      c(
        site_col,
        sample_col,
        replicate_col,
        otu_col
      )
    )
  )
  
  formula_predictors <- unique(
    c(
      occupancy_covars,
      all.vars(capture_formula),
      all.vars(abundance_formula)
    )
  )
  
  formula_predictors <- setdiff(
    formula_predictors,
    c(
      count_col,
      "total_reads",
      "z_sim",
      "a_sim"
    )
  )
  
  # IDs always factors
  long_df <- to_factor_cols(
    long_df,
    id_factor_cols
  )
  
  # Character covariates -> factors,
  # numeric covariates stay numeric
  long_df <- character_to_factor(
    long_df,
    formula_predictors
  )
  
  # ============================================================
  # 10. Set reference OTU
  # ============================================================
  
  otu_abundances <- long_df |>
    dplyr::group_by(
      .data[[otu_col]]
    ) |>
    dplyr::summarise(
      total_count = sum(
        .data[[count_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  top_otu <- otu_abundances |>
    dplyr::arrange(
      dplyr::desc(.data$total_count)
    ) |>
    dplyr::slice(1) |>
    dplyr::pull(
      .data[[otu_col]]
    )
  
  long_df[[otu_col]] <-
    stats::relevel(
      factor(long_df[[otu_col]]),
      ref = as.character(top_otu)
    )
  
  # ============================================================
  # 11. Define hierarchical keys
  # ============================================================
  
  site_keys <- unique(
    c(site_col, otu_col)
  )
  
  sample_keys <- unique(
    c(site_col, sample_col, otu_col)
  )
  
  pcr_keys <- unique(
    c(
      site_col,
      sample_col,
      replicate_col,
      otu_col
    )
  )
  
  pcr_keys <- pcr_keys[
    !is.na(pcr_keys)
  ]
  
  # ============================================================
  # 12. Build site-level latent occupancy data Z
  # ============================================================
  
  site_keep_vars <- if (
    !is.null(occupancy_covars)
  ) {
    
    occupancy_covars
    
  } else {
    
    character(0)
  }
  
  reduced_data <- long_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(site_keys)
      )
    ) |>
    dplyr::summarise(
      z_obs = as.integer(
        any(
          .data[[count_col]] >
            abundance_threshold
        )
      ),
      dplyr::across(
        dplyr::all_of(site_keep_vars),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      z_sim = .data$z_obs
    )
  
  reduced_data <- to_factor_cols(
    reduced_data,
    intersect(
      id_factor_cols,
      names(reduced_data)
    )
  )
  
  reduced_data <- character_to_factor(
    reduced_data,
    occupancy_covars
  )
  
  # ============================================================
  # 13. Build biological-sample latent capture data A
  # ============================================================
  
  cap_vars <- get_formula_vars(
    capture_formula,
    "a_sim"
  )
  
  sample_keep_vars <- setdiff(
    cap_vars,
    c(
      sample_keys,
      "a_sim"
    )
  )
  
  sample_keep_vars <- intersect(
    sample_keep_vars,
    names(long_df)
  )
  
  sample_data <- long_df |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(sample_keys)
      )
    ) |>
    dplyr::summarise(
      a_obs = as.integer(
        any(
          .data[[count_col]] >
            abundance_threshold
        )
      ),
      dplyr::across(
        dplyr::all_of(sample_keep_vars),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      a_sim = .data$a_obs
    )
  
  sample_data <- to_factor_cols(
    sample_data,
    intersect(
      id_factor_cols,
      names(sample_data)
    )
  )
  
  sample_data <- character_to_factor(
    sample_data,
    cap_vars
  )
  
  # ============================================================
  # 14. Formula validation
  # ============================================================
  
  validate_formula <- function(
    formula,
    data,
    expected_response,
    model_name
  ) {
    
    response <- all.vars(formula)[1]
    
    if (!identical(
      response,
      expected_response
    )) {
      
      stop(
        model_name,
        " formula response must be '",
        expected_response,
        "', but got '",
        response,
        "'."
      )
    }
    
    missing_vars <- setdiff(
      all.vars(formula),
      names(data)
    )
    
    if (length(missing_vars) > 0) {
      
      stop(
        "Missing variables in ",
        model_name,
        " formula: ",
        paste(
          missing_vars,
          collapse = ", "
        )
      )
    }
    
    invisible(TRUE)
  }
  
  validate_formula(
    capture_formula,
    sample_data,
    "a_sim",
    "capture"
  )
  
  validate_formula(
    abundance_formula,
    long_df,
    count_col,
    "abundance"
  )
  
  # ============================================================
  # 15. Abundance family setup
  # ============================================================
  
  abundance_glmm_family <- switch(
    abundance_family,
    
    poisson =
      stats::poisson(),
    
    nbinom =
      glmmTMB::nbinom2(),
    
    zip =
      stats::poisson(),
    
    zinb =
      glmmTMB::nbinom2()
  )
  
  zi_formula <- if (
    abundance_family %in% c(
      "zip",
      "zinb"
    )
  ) {
    
    ~1
    
  } else {
    
    ~0
  }
  
  # ============================================================
  # 16. Zero probability from abundance model
  # ============================================================
  
  get_p0_abundance <- function(
    fit,
    newdata,
    family_name
  ) {
    
    # For zero-inflated models we need the conditional count mean,
    # not the zero-inflation-adjusted response expectation.
    pred_type <- if (
      family_name %in% c(
        "zip",
        "zinb"
      )
    ) {
      
      "conditional"
      
    } else {
      
      "response"
    }
    
    mu <- as.numeric(
      stats::predict(
        fit,
        type = pred_type,
        newdata = newdata,
        allow.new.levels = TRUE
      )
    )
    
    mu <- pmax(
      mu,
      1e-12
    )
    
    zi <- rep(
      0,
      length(mu)
    )
    
    if (
      family_name %in% c(
        "zip",
        "zinb"
      )
    ) {
      
      zi <- tryCatch(
        as.numeric(
          stats::predict(
            fit,
            type = "zprob",
            newdata = newdata,
            allow.new.levels = TRUE
          )
        ),
        error = function(e) {
          rep(
            0,
            length(mu)
          )
        }
      )
      
      zi <- pmin(
        pmax(
          zi,
          0
        ),
        1
      )
    }
    
    p0_cond <- switch(
      family_name,
      
      poisson =
        stats::dpois(
          0,
          lambda = mu
        ),
      
      zip =
        stats::dpois(
          0,
          lambda = mu
        ),
      
      nbinom = {
        
        theta <- tryCatch(
          as.numeric(
            stats::sigma(fit)
          ),
          error = function(e) NA_real_
        )
        
        if (
          !is.finite(theta) ||
          theta <= 0
        ) {
          
          stop(
            "Invalid negative-binomial dispersion parameter."
          )
        }
        
        stats::dnbinom(
          0,
          mu = mu,
          size = theta
        )
      },
      
      zinb = {
        
        theta <- tryCatch(
          as.numeric(
            stats::sigma(fit)
          ),
          error = function(e) NA_real_
        )
        
        if (
          !is.finite(theta) ||
          theta <= 0
        ) {
          
          stop(
            "Invalid negative-binomial dispersion parameter."
          )
        }
        
        stats::dnbinom(
          0,
          mu = mu,
          size = theta
        )
      }
    )
    
    p0 <- zi +
      (1 - zi) *
      p0_cond
    
    pmin(
      pmax(
        p0,
        1e-12
      ),
      1
    )
  }
  
  # ============================================================
  # 17. Precompute site x OTU template
  # ============================================================
  
  z_template <- reshape2::acast(
    reduced_data,
    stats::as.formula(
      paste(
        bt(site_col),
        "~",
        bt(otu_col)
      )
    ),
    value.var = "z_sim",
    fill = 0
  )
  
  z_sites <- rownames(z_template)
  z_otus  <- colnames(z_template)
  
  if (length(z_sites) < 2) {
    stop(
      "At least two sites are required for the GLLVM occupancy model."
    )
  }
  
  # ============================================================
  # 18. Occupancy design matrix
  # ============================================================
  
  X_cov <- NULL
  
  if (
    !is.null(occupancy_covars) &&
    length(occupancy_covars) > 0
  ) {
    
    cov_cols <- unique(
      c(
        site_col,
        occupancy_covars
      )
    )
    
    cov_df <- reduced_data |>
      dplyr::select(
        dplyr::all_of(cov_cols)
      ) |>
      dplyr::distinct()
    
    cov_df[[site_col]] <-
      as.character(
        cov_df[[site_col]]
      )
    
    # Covariates must be constant within site
    cov_count <- cov_df |>
      dplyr::count(
        .data[[site_col]],
        name = ".n_cov_rows"
      )
    
    bad_sites <- cov_count |>
      dplyr::filter(
        .data$.n_cov_rows > 1
      )
    
    if (nrow(bad_sites) > 0) {
      
      stop(
        "occupancy_covars must have one unique combination per site. ",
        "Covariates vary within ",
        nrow(bad_sites),
        " site(s)."
      )
    }
    
    cov_df <- character_to_factor(
      cov_df,
      occupancy_covars
    )
    
    cov_df <- cov_df[
      match(
        z_sites,
        cov_df[[site_col]]
      ),
      ,
      drop = FALSE
    ]
    
    if (anyNA(cov_df[[site_col]])) {
      
      stop(
        "Could not align occupancy covariates with all GLLVM sites."
      )
    }
    
    X_full <- stats::model.matrix(
      ~ .,
      data = cov_df[
        ,
        occupancy_covars,
        drop = FALSE
      ]
    )
    
    # Remove intercept column only if present
    if (
      ncol(X_full) > 0 &&
      "(Intercept)" %in% colnames(X_full)
    ) {
      
      X_cov <- X_full[
        ,
        colnames(X_full) != "(Intercept)",
        drop = FALSE
      ]
      
    } else {
      
      X_cov <- X_full
    }
    
    if (ncol(X_cov) == 0) {
      X_cov <- NULL
    } else {
      # model.matrix() retains row positions but may assign numeric row names.
      # The GLLVM requires matching site identifiers in y and X.
      rownames(X_cov) <- z_sites
      if (anyNA(X_cov) || any(!is.finite(X_cov))) {
        stop("Occupancy design matrix contains missing/non-finite values.")
      }
    }
  }
  
  # ============================================================
  # 19. Storage
  # ============================================================
  
  psi_list <- vector(
    "list",
    n_iter
  )
  
  capture_list <- vector(
    "list",
    n_iter
  )
  
  lambda_list <- vector(
    "list",
    n_iter
  )
  
  p_detect_list <- vector(
    "list",
    n_iter
  )
  
  occupancy_models <- vector(
    "list",
    n_iter
  )
  
  capture_models <- vector(
    "list",
    n_iter
  )
  
  abundance_models <- vector(
    "list",
    n_iter
  )
  
  lv_sites_list <- vector(
    "list",
    n_iter
  )
  
  lv_species_list <- vector(
    "list",
    n_iter
  )
  
  iteration_success <- rep(
    FALSE,
    n_iter
  )
  
  diagnostic_AIC <- data.frame(
    iteration = integer(),
    occupancy_AIC = numeric(),
    capture_AIC = numeric(),
    abundance_AIC = numeric(),
    abundance_logLik = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Numerical diagnostics for each attempted occupancy fit.
  # These are diagnostics, not a substitute for examining chain stability.
  gllvm_diagnostics <- data.frame(
    iteration = seq_len(n_iter),
    optimizer_converged = rep(NA, n_iter),
    max_abs_gradient = rep(NA_real_, n_iter),
    hessian_min_eigen = rep(NA_real_, n_iter),
    hessian_negative_eigen = rep(NA_integer_, n_iter),
    hessian_near_zero_eigen = rep(NA_integer_, n_iter),
    covariance_min_eigen = rep(NA_real_, n_iter),
    covariance_negative_eigen = rep(NA_integer_, n_iter),
    covariance_near_zero_eigen = rep(NA_integer_, n_iter),
    numerical_ok = rep(NA, n_iter),
    n_all_present = rep(NA_integer_, n_iter),
    n_all_absent = rep(NA_integer_, n_iter),
    n_gllvm_otus = rep(NA_integer_, n_iter),
    n_excluded_otus = rep(NA_integer_, n_iter),
    hessian_condition_number = rep(NA_real_, n_iter),
    hessian_ill_conditioned = rep(NA, n_iter),
    hessian_positive_definite = rep(NA, n_iter),
    stringsAsFactors = FALSE
  )
  # Preserve the exact OTU identities for each augmented occupancy matrix.
  gllvm_separation <- vector("list", n_iter)
  names(gllvm_separation) <- as.character(seq_len(n_iter))
  gllvm_excluded_otus <- vector("list", n_iter)
  names(gllvm_excluded_otus) <- as.character(seq_len(n_iter))
  gllvm_fit_warnings <- vector("list", n_iter)
  names(gllvm_fit_warnings) <- as.character(seq_len(n_iter))
  
  
  matrix_eigen_diagnostic <- function(M, tolerance = gllvm_eigen_tol) {
    empty <- c(minimum = NA_real_, negative = NA_real_, near_zero = NA_real_)
    if (is.null(M) || !is.matrix(M) || nrow(M) != ncol(M) ||
        any(!is.finite(M))) return(empty)
    M <- (M + t(M)) / 2
    ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NULL)
    if (is.null(ev) || !length(ev)) return(empty)
    threshold <- tolerance * max(abs(ev))
    c(minimum = min(ev),
      negative = sum(ev < -threshold),
      near_zero = sum(abs(ev) <= threshold))
  }
  
  # ============================================================
  # 20. EM-like stochastic iterations
  # ============================================================
  
  for (i in seq_len(n_iter)) {
    
    if (verbose) {
      message(
        "Iteration ",
        i
      )
    }
    
    t_iter <- Sys.time()
    
    # ----------------------------------------------------------
    # 20.1 Update occupancy matrix
    # ----------------------------------------------------------
    
    z_matrix <- z_template
    z_matrix[,] <- 0
    
    idx <- cbind(
      match(
        as.character(
          reduced_data[[site_col]]
        ),
        z_sites
      ),
      match(
        as.character(
          reduced_data[[otu_col]]
        ),
        z_otus
      )
    )
    
    if (anyNA(idx)) {
      
      warning(
        "Could not align site/OTU indices at iteration ",
        i,
        ". Iteration skipped."
      )
      
      next
    }
    
    z_matrix[idx] <-
      reduced_data$z_sim
    
    # Constant columns in the AUGMENTED occupancy matrix may cause
    # complete separation of species-specific binomial intercepts.
    # Keep these OTUs; record the issue rather than changing the model.
    occupied <- colSums(z_matrix == 1, na.rm = TRUE)
    n_valid <- colSums(!is.na(z_matrix))
    all_present <- colnames(z_matrix)[n_valid > 0 & occupied == n_valid]
    all_absent <- colnames(z_matrix)[n_valid > 0 & occupied == 0]
    gllvm_separation[[i]] <- list(
      iteration = i,
      all_present = all_present,
      all_absent = all_absent,
      occupied_sites = stats::setNames(as.integer(occupied), colnames(z_matrix)),
      n_sites = nrow(z_matrix)
    )
    gllvm_diagnostics$n_all_present[i] <- length(all_present)
    gllvm_diagnostics$n_all_absent[i] <- length(all_absent)
    if (gllvm_warn_separation && (length(all_present) || length(all_absent))) {
      warning(sprintf(
        "GLLVM iteration %d: augmented occupancy has %d universally present OTUs (%s) and %d universally absent OTUs (%s).",
        i, length(all_present), paste(all_present, collapse = ", "),
        length(all_absent), paste(all_absent, collapse = ", ")
      ), call. = FALSE)
    }
    
    # Automatic, iteration-specific screening of constant augmented
    # occupancy columns. Optional min_occupied_sites > 0 also removes
    # rare / near-universal columns (sensitivity analysis only).
    # The full z_matrix is NEVER modified.
    excluded <- character()
    if (gllvm_auto_exclude) {
      excluded <- colnames(z_matrix)[
        n_valid == 0L |
          occupied <= gllvm_min_occupied_sites |
          (n_valid - occupied) <= gllvm_min_occupied_sites
      ]
    }
    included <- colnames(z_matrix)[
      !colnames(z_matrix) %in% excluded
    ]
    gllvm_excluded_otus[[i]] <- excluded
    gllvm_diagnostics$n_excluded_otus[i] <- length(excluded)
    gllvm_diagnostics$n_gllvm_otus[i] <- length(included)
    
    if (length(included) < max(3L, num_lv + 2L)) {
      warning(sprintf(
        "Iteration %d: only %d OTUs remain after screening; skipping.",
        i, length(included)
      ), call. = FALSE)
      next
    }
    z_gllvm <- z_matrix[, included, drop = FALSE]
    # Verify covariate identity rather than silently relabeling rows.
    if (!is.null(X_cov)) {
      stopifnot(nrow(X_cov) == nrow(z_gllvm),
                identical(rownames(X_cov), rownames(z_gllvm)))
    }
    if (verbose && length(excluded)) {
      message(sprintf(
        "GLLVM iteration %d: excluding %d OTUs from GLLVM only: %s",
        i, length(excluded), paste(excluded, collapse = ", ")
      ))
    }
    
    # ----------------------------------------------------------
    # 20.2 GLLVM occupancy model
    # ----------------------------------------------------------
    
    fit_warnings <- character()
    model_occupancy <- tryCatch(
      withCallingHandlers(
        gllvm::gllvm(
          y = z_gllvm,
          X = X_cov,
          family = "binomial",
          num.lv = num_lv,
          method = "VA",
          control.start = list(
            n.init = gllvm_n_init,
            starting.val = "res"
          )
        ),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        
        message(
          "GLLVM failed at iteration ",
          i,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    gllvm_fit_warnings[[i]] <- unique(fit_warnings)
    if (is.null(model_occupancy)) {
      next
    }
    
    # gllvm 2.0.x reports convergence as logical TRUE/FALSE.
    # Some other optimizer interfaces use numeric code 0 for success.
    conv <- model_occupancy$convergence
    optimizer_ok <- if (is.logical(conv) && length(conv) == 1L) {
      isTRUE(conv)
    } else if (is.numeric(conv) && length(conv) == 1L && !is.na(conv)) {
      conv == 0
    } else {
      NA
    }
    
    # Use the parameter vector explicitly, rather than relying on the
    # objective closure's current environment state.
    gradient <- tryCatch({
      fn <- model_occupancy$TMBfn
      if (is.null(fn) || !is.function(fn$gr)) stop("No gradient function")
      pars <- fn$par
      if (is.null(pars)) fn$gr() else fn$gr(pars)
    }, error = function(e) NA_real_)
    max_gradient <- if (length(gradient) && all(is.finite(gradient))) {
      max(abs(gradient))
    } else NA_real_
    
    hess_diag <- matrix_eigen_diagnostic(model_occupancy$Hess$Hess.full)
    cov_diag <- matrix_eigen_diagnostic(model_occupancy$Hess$cov.mat.mod)
    
    # Condition number: Inf when Hessian is non-positive-definite.
    H <- model_occupancy$Hess$Hess.full
    hessian_condition <- NA_real_
    hessian_pd <- NA
    if (is.matrix(H) && nrow(H) == ncol(H) && all(is.finite(H))) {
      hessian_ev <- tryCatch(
        eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values,
        error = function(e) NULL
      )
      if (!is.null(hessian_ev) && length(hessian_ev)) {
        hessian_pd <- min(hessian_ev) > 0
        hessian_condition <- if (hessian_pd) {
          max(hessian_ev) / min(hessian_ev)
        } else Inf
      }
    }
    gllvm_diagnostics$hessian_condition_number[i] <- hessian_condition
    gllvm_diagnostics$hessian_positive_definite[i] <- hessian_pd
    gllvm_diagnostics$hessian_ill_conditioned[i] <- if (is.na(hessian_condition)) {
      NA
    } else hessian_condition > gllvm_condition_tol
    
    # numerical_ok retains the original optimizer/gradient criterion.
    # Ill-conditioning is separately recorded, not used to silently reject
    # stochastic updates or alter the estimator.
    # Missing Hessian/covariance information is not evidence of a bad fit.
    # Conversely, an observed negative eigenvalue or ill-conditioned Hessian
    # is a genuine reason to flag the iteration.
    hessian_bad <- (is.finite(hess_diag["negative"]) &&
                      hess_diag["negative"] > 0) ||
      (is.finite(hess_diag["near_zero"]) &&
         hess_diag["near_zero"] > 0) ||
      identical(hessian_pd, FALSE) ||
      (is.finite(hessian_condition) &&
         hessian_condition > gllvm_condition_tol)
    covariance_bad <- is.finite(cov_diag["negative"]) &&
      cov_diag["negative"] > 0
    se_warning <- any(grepl(
      "negative variance|not converged|not positive definite|non-positive definite",
      fit_warnings, ignore.case = TRUE
    ))
    numerical_ok <- isTRUE(optimizer_ok) && is.finite(max_gradient) &&
      max_gradient <= gllvm_gradient_tol &&
      !hessian_bad && !covariance_bad && !se_warning
    
    gllvm_diagnostics$optimizer_converged[i] <- optimizer_ok
    gllvm_diagnostics$max_abs_gradient[i] <- max_gradient
    gllvm_diagnostics$hessian_min_eigen[i] <- hess_diag["minimum"]
    gllvm_diagnostics$hessian_negative_eigen[i] <- as.integer(hess_diag["negative"])
    gllvm_diagnostics$hessian_near_zero_eigen[i] <- as.integer(hess_diag["near_zero"])
    gllvm_diagnostics$covariance_min_eigen[i] <- cov_diag["minimum"]
    gllvm_diagnostics$covariance_negative_eigen[i] <- as.integer(cov_diag["negative"])
    gllvm_diagnostics$covariance_near_zero_eigen[i] <- as.integer(cov_diag["near_zero"])
    gllvm_diagnostics$numerical_ok[i] <- numerical_ok
    
    if (verbose) {
      message(sprintf(
        "GLLVM iteration %d: optimizer=%s; max gradient=%s; Hessian min=%s; covariance min=%s",
        i, as.character(optimizer_ok),
        format(max_gradient, digits = 4),
        format(hess_diag["minimum"], digits = 4),
        format(cov_diag["minimum"], digits = 4)
      ))
    }
    if (!numerical_ok) {
      if (verbose) message(sprintf(
        "GLLVM iteration %d flagged: optimizer=%s; max gradient=%s; Hessian issue=%s; covariance issue=%s; SE warning=%s",
        i, as.character(optimizer_ok), format(max_gradient, digits = 4),
        hessian_bad, covariance_bad, se_warning
      ))
      # Strict mode rejects the iteration. Non-strict mode keeps the
      # stochastic update but retains an explicit diagnostic flag.
      if (gllvm_strict) next
    }
    
    occupancy_models[[i]] <-
      model_occupancy
    
    # ----------------------------------------------------------
    # 20.3 GLLVM latent site variables
    # ----------------------------------------------------------
    
    if (
      !is.null(model_occupancy$lvs) &&
      length(model_occupancy$lvs) > 0 &&
      ncol(
        as.data.frame(
          model_occupancy$lvs
        )
      ) > 0
    ) {
      
      lv_sites <- as.data.frame(
        model_occupancy$lvs
      )
      
      colnames(lv_sites) <-
        paste0(
          "LV",
          seq_len(
            ncol(lv_sites)
          )
        )
      
      site_names <- rownames(
        model_occupancy$lvs
      )
      
      if (
        is.null(site_names) ||
        length(site_names) != nrow(lv_sites)
      ) {
        
        site_names <- z_sites
      }
      
      lv_sites[[site_col]] <-
        as.character(site_names)
      
      lv_sites$Iteration <- i
      
    } else {
      
      lv_sites <- data.frame()
    }
    
    # ----------------------------------------------------------
    # 20.4 Taxon latent loadings
    # ----------------------------------------------------------
    
    theta_obj <- tryCatch(
      model_occupancy$params$theta,
      error = function(e) NULL
    )
    
    if (
      !is.null(theta_obj) &&
      length(theta_obj) > 0
    ) {
      
      theta_mat <- as.matrix(
        theta_obj
      )
      
      if (
        nrow(theta_mat) > 0 &&
        ncol(theta_mat) > 0
      ) {
        
        lv_species <- as.data.frame(
          theta_mat
        )
        
        colnames(lv_species) <-
          paste0(
            "LV",
            seq_len(
              ncol(lv_species)
            )
          )
        
        otu_names <- rownames(
          theta_mat
        )
        
        if (
          is.null(otu_names) ||
          length(otu_names) != nrow(lv_species)
        ) {
          
          otu_names <- included[
            seq_len(
              min(
                length(z_otus),
                nrow(lv_species)
              )
            )
          ]
        }
        
        if (
          length(otu_names) ==
          nrow(lv_species)
        ) {
          
          lv_species[[otu_col]] <-
            as.character(otu_names)
          
        } else {
          
          lv_species[[otu_col]] <-
            as.character(
              seq_len(
                nrow(lv_species)
              )
            )
        }
        
        lv_species$Iteration <- i
        
      } else {
        
        lv_species <- data.frame()
      }
      
    } else {
      
      lv_species <- data.frame()
    }
    
    lv_sites_list[[i]] <-
      lv_sites
    
    lv_species_list[[i]] <-
      lv_species
    
    # ----------------------------------------------------------
    # 20.5 Occupancy probability
    # ----------------------------------------------------------
    
    psi_prob <- tryCatch(
      
      stats::predict(
        model_occupancy,
        type = "response"
      ),
      
      error = function(e) {
        
        message(
          "GLLVM prediction failed at iteration ",
          i,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (is.null(psi_prob)) {
      next
    }
    
    psi_reduced <- as.matrix(psi_prob)
    
    if (!identical(dim(psi_reduced), dim(z_gllvm))) {
      warning(sprintf(
        "Iteration %d: reduced GLLVM prediction dimensions mismatch; skipping.",
        i
      ), call. = FALSE)
      next
    }
    
    # Map GLLVM probabilities back to all OTUs in their original order.
    psi_prob <- matrix(
      NA_real_, nrow = nrow(z_matrix), ncol = ncol(z_matrix),
      dimnames = dimnames(z_matrix)
    )
    psi_prob[, included] <- psi_reduced
    
    # Excluded OTUs receive smoothed intercept-only probabilities.
    # This is a Beta-Binomial empirical approximation conditional on
    # the CURRENT augmented Z, not a joint hierarchical estimate.
    # It avoids assigning 0 or 1 and preserves all OTUs downstream.
    if (length(excluded)) {
      a <- gllvm_excluded_prior[1L]
      b <- gllvm_excluded_prior[2L]
      for (otu in excluded) {
        n <- n_valid[[otu]]
        k <- occupied[[otu]]
        psi_prob[, otu] <- (k + a) / (n + a + b)
      }
    }
    if (anyNA(psi_prob) || any(!is.finite(psi_prob))) {
      warning(sprintf("Iteration %d: invalid occupancy probabilities; skipping.", i),
              call. = FALSE)
      next
    }
    
    
    # ----------------------------------------------------------
    # 20.5 Convert occupancy probabilities to long format
    # ----------------------------------------------------------
    # expand.grid() uses the first column as the fastest-varying
    # index, matching the column-major ordering of as.vector().
    stopifnot(
      is.matrix(psi_prob),
      !is.null(rownames(psi_prob)),
      !is.null(colnames(psi_prob)),
      !anyNA(rownames(psi_prob)),
      !anyNA(colnames(psi_prob)),
      !anyDuplicated(rownames(psi_prob)),
      !anyDuplicated(colnames(psi_prob)),
      all(is.finite(psi_prob))
    )
    
    psi_long <- expand.grid(
      rownames(psi_prob),
      colnames(psi_prob),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    names(psi_long) <- c(site_col, otu_col)
    psi_long$psi_prob <- as.vector(psi_prob)
    
    # Use a local character-key copy so reduced_data retains its
    # original factor columns for subsequent model operations.
    reduced_keys <- reduced_data[, site_keys, drop = FALSE]
    for (key in site_keys) {
      reduced_keys[[key]] <- as.character(reduced_keys[[key]])
      psi_long[[key]] <- as.character(psi_long[[key]])
    }
    
    if (anyDuplicated(psi_long[, site_keys, drop = FALSE])) {
      stop(sprintf("Iteration %d: duplicate occupancy prediction keys.", i))
    }
    
    missing_keys <- dplyr::anti_join(
      dplyr::distinct(reduced_keys),
      dplyr::distinct(psi_long[, site_keys, drop = FALSE]),
      by = site_keys
    )
    if (nrow(missing_keys) > 0L) {
      warning(sprintf(
        "Iteration %d: %d site-OTU keys lack occupancy predictions; skipping.",
        i, nrow(missing_keys)
      ), call. = FALSE)
      next
    }
    
    psi_long$psi_prob <- pmin(
      pmax(as.numeric(psi_long$psi_prob), 1e-6),
      1 - 1e-6
    )
    psi_long$eta <- stats::qlogis(psi_long$psi_prob)
    psi_list[[i]] <- psi_long[, c(site_col, otu_col, "eta"), drop = FALSE]
    
    # ----------------------------------------------------------
    # 20.6 Add current Z state to sample data
    # ----------------------------------------------------------
    
    sample_data <- sample_data |>
      dplyr::select(
        -dplyr::any_of(
          c(
            "z_sim",
            "p0_sample",
            "capture_prob"
          )
        )
      ) |>
      dplyr::left_join(
        reduced_data |>
          dplyr::select(
            dplyr::all_of(
              c(
                site_keys,
                "z_sim"
              )
            )
          ),
        by = site_keys
      )
    
    capture_fit_data <- sample_data |>
      dplyr::filter(
        .data$z_sim == 1
      )
    
    if (nrow(capture_fit_data) == 0) {
      
      message(
        "No rows available for capture model at iteration ",
        i,
        ". Iteration skipped."
      )
      
      next
    }
    
    # ----------------------------------------------------------
    # 20.7 Capture GLMM
    # ----------------------------------------------------------
    
    cap_fit <- tryCatch(
      
      glmmTMB::glmmTMB(
        formula = capture_formula,
        data = capture_fit_data,
        family = stats::binomial()
      ),
      
      error = function(e) {
        
        message(
          "Capture model failed at iteration ",
          i,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (is.null(cap_fit)) {
      next
    }
    
    capture_models[[i]] <-
      cap_fit
    
    eta_capture <- tryCatch(
      
      as.numeric(
        stats::predict(
          cap_fit,
          type = "link",
          newdata = sample_data,
          allow.new.levels = TRUE
        )
      ),
      
      error = function(e) {
        
        message(
          "Capture prediction failed at iteration ",
          i,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (
      is.null(eta_capture) ||
      length(eta_capture) != nrow(sample_data)
    ) {
      
      next
    }
    
    capture_prob <-
      stats::plogis(
        eta_capture
      )
    
    capture_list[[i]] <-
      data.frame(
        sample_data[
          sample_keys
        ],
        eta = eta_capture,
        check.names = FALSE
      )
    
    # ----------------------------------------------------------
    # 20.8 Abundance data conditional on A = 1
    # ----------------------------------------------------------
    
    abundance_data <- long_df |>
      dplyr::left_join(
        sample_data |>
          dplyr::select(
            dplyr::all_of(
              c(
                sample_keys,
                "a_sim"
              )
            )
          ),
        by = sample_keys
      ) |>
      dplyr::filter(
        .data$a_sim == 1
      )
    
    if (nrow(abundance_data) == 0) {
      
      message(
        "No rows available for abundance model at iteration ",
        i,
        ". Iteration skipped."
      )
      
      next
    }
    
    abundance_data <- to_factor_cols(
      abundance_data,
      intersect(
        id_factor_cols,
        names(abundance_data)
      )
    )
    
    abundance_data <- character_to_factor(
      abundance_data,
      formula_predictors
    )
    
    # ----------------------------------------------------------
    # 20.9 Abundance GLMM
    # ----------------------------------------------------------
    
    model_abundance <- tryCatch(
      
      glmmTMB::glmmTMB(
        formula = abundance_formula,
        data = abundance_data,
        family = abundance_glmm_family,
        ziformula = zi_formula
      ),
      
      error = function(e) {
        
        message(
          "Abundance model failed at iteration ",
          i,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (is.null(model_abundance)) {
      next
    }
    
    abundance_models[[i]] <-
      model_abundance
    
    # ----------------------------------------------------------
    # 20.10 Abundance predictions
    # ----------------------------------------------------------
    
    lambda_eta <- tryCatch(
      
      as.numeric(
        stats::predict(
          model_abundance,
          type = "link",
          newdata = long_df,
          allow.new.levels = TRUE
        )
      ),
      
      error = function(e) {
        
        message(
          "Abundance prediction failed at iteration ",
          i,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (
      is.null(lambda_eta) ||
      length(lambda_eta) != nrow(long_df)
    ) {
      
      next
    }
    
    lambda_list[[i]] <-
      data.frame(
        long_df[
          site_keys
        ],
        eta = lambda_eta,
        check.names = FALSE
      ) |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            site_keys
          )
        )
      ) |>
      dplyr::summarise(
        eta = mean(
          .data$eta,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    # ----------------------------------------------------------
    # 20.11 PCR/read-level zero probabilities
    # ----------------------------------------------------------
    
    pred_data_all <- long_df |>
      dplyr::left_join(
        sample_data |>
          dplyr::select(
            dplyr::all_of(
              c(
                sample_keys,
                "a_sim"
              )
            )
          ),
        by = sample_keys
      )
    
    pred_data_all <- to_factor_cols(
      pred_data_all,
      intersect(
        id_factor_cols,
        names(pred_data_all)
      )
    )
    
    pred_data_all <- character_to_factor(
      pred_data_all,
      formula_predictors
    )
    
    p0_pcr <- tryCatch(
      
      get_p0_abundance(
        fit = model_abundance,
        newdata = pred_data_all,
        family_name = abundance_family
      ),
      
      error = function(e) {
        
        message(
          "Zero-probability calculation failed at iteration ",
          i,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    if (
      is.null(p0_pcr) ||
      length(p0_pcr) != nrow(pred_data_all)
    ) {
      
      next
    }
    
    # ----------------------------------------------------------
    # 20.12 Collapse PCR zeros to sample level
    # ----------------------------------------------------------
    
    sample_p0 <- pred_data_all |>
      dplyr::mutate(
        p0_pcr = p0_pcr
      ) |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            sample_keys
          )
        )
      ) |>
      dplyr::summarise(
        p0_sample = prod(
          .data$p0_pcr,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    sample_data <- sample_data |>
      dplyr::select(
        -dplyr::any_of(
          c(
            "p0_sample",
            "capture_prob"
          )
        )
      ) |>
      dplyr::mutate(
        capture_prob =
          capture_prob
      ) |>
      dplyr::left_join(
        sample_p0,
        by = sample_keys
      )
    
    sample_data$p0_sample[
      is.na(sample_data$p0_sample)
    ] <- 1
    
    # ----------------------------------------------------------
    # 20.13 Update latent A
    # ----------------------------------------------------------
    
    zero_sample <- which(
      sample_data$a_obs == 0 &
        sample_data$z_sim == 1
    )
    
    if (length(zero_sample) > 0) {
      
      numerator_a <-
        sample_data$capture_prob[
          zero_sample
        ] *
        sample_data$p0_sample[
          zero_sample
        ]
      
      denominator_a <-
        (
          1 -
            sample_data$capture_prob[
              zero_sample
            ]
        ) +
        numerator_a
      
      posterior_a <-
        numerator_a /
        pmax(
          denominator_a,
          1e-12
        )
      
      posterior_a <-
        pmin(
          pmax(
            posterior_a,
            1e-4
          ),
          1 - 1e-4
        )
      
      sample_data$a_sim[
        zero_sample
      ] <- stats::rbinom(
        n = length(
          zero_sample
        ),
        size = 1,
        prob = posterior_a
      )
    }
    
    # Positive observations remain positive
    sample_data$a_sim[
      sample_data$a_obs == 1
    ] <- 1
    
    # Cannot capture if site absent
    sample_data$a_sim[
      sample_data$z_sim == 0
    ] <- 0
    
    # ----------------------------------------------------------
    # 20.14 Collapse sample probabilities to site level
    # ----------------------------------------------------------
    
    site_p0 <- sample_data |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            site_keys
          )
        )
      ) |>
      dplyr::summarise(
        p0_site = prod(
          (
            1 -
              .data$capture_prob
          ) +
            .data$capture_prob *
            .data$p0_sample,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    site_p0$p0_site <-
      pmin(
        pmax(
          site_p0$p0_site,
          0
        ),
        1
      )
    
    # ----------------------------------------------------------
    # 20.15 Detection probability on cloglog scale
    # ----------------------------------------------------------
    
    p_detect_list[[i]] <-
      site_p0 |>
      dplyr::mutate(
        .p0_safe = pmin(
          pmax(
            .data$p0_site,
            1e-12
          ),
          1 - 1e-12
        ),
        eta = log(
          -log(
            .data$.p0_safe
          )
        )
      ) |>
      dplyr::select(
        dplyr::all_of(
          site_keys
        ),
        "eta"
      )
    
    # ----------------------------------------------------------
    # 20.16 Update site-level p0
    # ----------------------------------------------------------
    
    reduced_data <- reduced_data |>
      dplyr::select(
        -dplyr::any_of(
          "p0_site"
        )
      ) |>
      dplyr::left_join(
        site_p0,
        by = site_keys
      ) |>
      dplyr::mutate(
        p0_site =
          tidyr::replace_na(
            .data$p0_site,
            1
          )
      )
    
    # ----------------------------------------------------------
    # 20.17 Merge GLLVM psi
    # ----------------------------------------------------------
    
    
    # Match by character labels, without modifying the original
    # factor-valued identifiers in reduced_data.
    z_merge <- reduced_data
    for (key in site_keys) {
      z_merge[[key]] <- as.character(z_merge[[key]])
    }
    z_merge <- dplyr::left_join(
      z_merge,
      psi_long[, c(site_keys, "psi_prob"), drop = FALSE],
      by = site_keys
    )
    if (nrow(z_merge) != nrow(reduced_data) ||
        anyNA(z_merge$psi_prob) || any(!is.finite(z_merge$psi_prob))) {
      warning(sprintf(
        "Iteration %d: missing, non-finite, or duplicated occupancy predictions; skipping.", i
      ), call. = FALSE)
      next
    }
    
    # ----------------------------------------------------------
    # 20.18 Update latent Z
    # ----------------------------------------------------------
    
    zero_indices <- which(
      z_merge$z_obs == 0
    )
    
    if (length(zero_indices) > 0) {
      
      numerator_z <-
        z_merge$psi_prob[
          zero_indices
        ] *
        z_merge$p0_site[
          zero_indices
        ]
      
      denominator_z <-
        (
          1 -
            z_merge$psi_prob[
              zero_indices
            ]
        ) +
        numerator_z
      
      posterior_z <-
        numerator_z /
        pmax(
          denominator_z,
          1e-12
        )
      
      posterior_z <-
        pmin(
          pmax(
            posterior_z,
            1e-4
          ),
          1 - 1e-4
        )
      
      reduced_data$z_sim[
        zero_indices
      ] <- stats::rbinom(
        n = length(
          zero_indices
        ),
        size = 1,
        prob = posterior_z
      )
    }
    
    # Positive site detections remain occupied
    reduced_data$z_sim[
      reduced_data$z_obs == 1
    ] <- 1
    
    # ----------------------------------------------------------
    # 20.19 Diagnostics
    # ----------------------------------------------------------
    
    diagnostic_AIC <- rbind(
      diagnostic_AIC,
      data.frame(
        iteration = i,
        
        occupancy_AIC =
          tryCatch(
            stats::AIC(
              model_occupancy
            ),
            error = function(e) {
              NA_real_
            }
          ),
        
        capture_AIC =
          tryCatch(
            stats::AIC(
              cap_fit
            ),
            error = function(e) {
              NA_real_
            }
          ),
        
        abundance_AIC =
          tryCatch(
            stats::AIC(
              model_abundance
            ),
            error = function(e) {
              NA_real_
            }
          ),
        
        abundance_logLik =
          tryCatch(
            as.numeric(
              stats::logLik(
                model_abundance
              )
            ),
            error = function(e) {
              NA_real_
            }
          )
      )
    )
    
    iteration_success[i] <- TRUE
    
    if (verbose) {
      
      message(
        "Iteration ",
        i,
        " finished in ",
        round(
          difftime(
            Sys.time(),
            t_iter,
            units = "secs"
          ),
          2
        ),
        " seconds"
      )
    }
  }
  
  # ============================================================
  # 21. Retained successful iterations
  # ============================================================
  
  keep_requested <- seq.int(
    burn_in + 1,
    n_iter
  )
  
  keep <- keep_requested[
    iteration_success[
      keep_requested
    ]
  ]
  
  if (length(keep) == 0) {
    
    stop(
      "No successful iterations remained after burn-in."
    )
  }
  
  if (
    verbose &&
    length(keep) <
    length(keep_requested)
  ) {
    
    message(
      length(keep),
      " of ",
      length(keep_requested),
      " post-burn-in iterations were fully successful."
    )
  }
  
  flagged <- which(!is.na(gllvm_diagnostics$numerical_ok) &
                     !gllvm_diagnostics$numerical_ok)
  if (length(flagged)) {
    warning(sprintf(
      "GLLVM numerical diagnostics flagged %d of %d attempted fits; inspect $gllvm_diagnostics and $gllvm_fit_warnings (iterations: %s).",
      length(flagged), n_iter, paste(flagged, collapse = ", ")
    ), call. = FALSE)
  }
  
  # ============================================================
  # 22. Generic link-scale summarizer
  # ============================================================
  
  summarise_link <- function(
    lst,
    link_type = c(
      "logit",
      "log",
      "cloglog"
    )
  ) {
    
    link_type <- match.arg(
      link_type
    )
    
    df <- dplyr::bind_rows(
      lst
    )
    
    if (nrow(df) == 0) {
      return(
        data.frame()
      )
    }
    
    inv_link <- switch(
      link_type,
      
      logit =
        stats::plogis,
      
      log =
        exp,
      
      cloglog =
        function(x) {
          1 - exp(-exp(x))
        }
    )
    
    df |>
      dplyr::mutate(
        value =
          inv_link(
            .data$eta
          )
      ) |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            site_keys
          )
        )
      ) |>
      dplyr::summarise(
        mean = mean(
          .data$value,
          na.rm = TRUE
        ),
        median = stats::median(
          .data$value,
          na.rm = TRUE
        ),
        lwr = stats::quantile(
          .data$value,
          0.025,
          na.rm = TRUE,
          names = FALSE
        ),
        upr = stats::quantile(
          .data$value,
          0.975,
          na.rm = TRUE,
          names = FALSE
        ),
        .groups = "drop"
      )
  }
  
  # ============================================================
  # 23. Capture summaries
  # ============================================================
  
  summarise_capture <- function(
    lst
  ) {
    
    df <- dplyr::bind_rows(
      lst
    )
    
    if (nrow(df) == 0) {
      return(
        data.frame()
      )
    }
    
    df |>
      dplyr::mutate(
        capture =
          stats::plogis(
            .data$eta
          )
      ) |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            sample_keys
          )
        )
      ) |>
      dplyr::summarise(
        capture_mean =
          mean(
            .data$capture,
            na.rm = TRUE
          ),
        
        capture_median =
          stats::median(
            .data$capture,
            na.rm = TRUE
          ),
        
        capture_lwr =
          stats::quantile(
            .data$capture,
            0.025,
            na.rm = TRUE,
            names = FALSE
          ),
        
        capture_upr =
          stats::quantile(
            .data$capture,
            0.975,
            na.rm = TRUE,
            names = FALSE
          ),
        
        .groups = "drop"
      )
  }
  
  # ============================================================
  # 24. Site-level capture summaries
  # ============================================================
  
  summarise_capture_site <- function(
    capture_summary
  ) {
    
    if (nrow(capture_summary) == 0) {
      return(
        data.frame()
      )
    }
    
    capture_summary |>
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(
            site_keys
          )
        )
      ) |>
      dplyr::summarise(
        capture_mean =
          mean(
            .data$capture_mean,
            na.rm = TRUE
          ),
        
        capture_median =
          stats::median(
            .data$capture_median,
            na.rm = TRUE
          ),
        
        capture_lwr =
          stats::quantile(
            .data$capture_lwr,
            0.025,
            na.rm = TRUE,
            names = FALSE
          ),
        
        capture_upr =
          stats::quantile(
            .data$capture_upr,
            0.975,
            na.rm = TRUE,
            names = FALSE
          ),
        
        .groups = "drop"
      )
  }
  
  # ============================================================
  # 25. Occupancy summaries
  # ============================================================
  
  psi_summary <-
    summarise_link(
      psi_list[keep],
      "logit"
    ) |>
    dplyr::rename_with(
      ~ paste0(
        "psi_",
        .x
      ),
      -dplyr::all_of(
        site_keys
      )
    )
  
  # ============================================================
  # 26. Abundance summaries
  # ============================================================
  
  lambda_summary <-
    summarise_link(
      lambda_list[keep],
      "log"
    ) |>
    dplyr::rename_with(
      ~ paste0(
        "lambda_",
        .x
      ),
      -dplyr::all_of(
        site_keys
      )
    )
  
  # ============================================================
  # 27. Detection summaries
  # ============================================================
  
  p_detect_summary <-
    summarise_link(
      p_detect_list[keep],
      "cloglog"
    ) |>
    dplyr::rename_with(
      ~ paste0(
        "p_detect_",
        .x
      ),
      -dplyr::all_of(
        site_keys
      )
    )
  
  # ============================================================
  # 28. Capture summaries
  # ============================================================
  
  capture_summary <-
    summarise_capture(
      capture_list[keep]
    )
  
  capture_site_summary <-
    summarise_capture_site(
      capture_summary
    )
  
  # ============================================================
  # 29. Final site x OTU summary
  # ============================================================
  
  final_summary <- psi_summary
  
  if (nrow(capture_site_summary) > 0) {
    
    final_summary <-
      dplyr::left_join(
        final_summary,
        capture_site_summary,
        by = site_keys
      )
  }
  
  if (nrow(lambda_summary) > 0) {
    
    final_summary <-
      dplyr::left_join(
        final_summary,
        lambda_summary,
        by = site_keys
      )
  }
  
  if (nrow(p_detect_summary) > 0) {
    
    final_summary <-
      dplyr::left_join(
        final_summary,
        p_detect_summary,
        by = site_keys
      )
  }
  
  # ============================================================
  # 30. Combine GLLVM latent variables
  # ============================================================
  
  lv_sites_combined <-
    dplyr::bind_rows(
      lv_sites_list[keep]
    )
  
  lv_species_combined <-
    dplyr::bind_rows(
      lv_species_list[keep]
    )
  
  # ============================================================
  # 31. Mean site latent scores
  # ============================================================
  
  mean_lv_sites <- if (
    nrow(lv_sites_combined) > 0 &&
    site_col %in% names(lv_sites_combined)
  ) {
    
    lv_sites_combined |>
      dplyr::group_by(
        .data[[site_col]]
      ) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with(
            "LV"
          ),
          ~ mean(
            .x,
            na.rm = TRUE
          )
        ),
        .groups = "drop"
      )
    
  } else {
    
    data.frame()
  }
  
  # ============================================================
  # 32. Mean OTU latent loadings
  # ============================================================
  
  mean_lv_species <- if (
    nrow(lv_species_combined) > 0 &&
    otu_col %in% names(lv_species_combined)
  ) {
    
    lv_species_combined |>
      dplyr::group_by(
        .data[[otu_col]]
      ) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with(
            "LV"
          ),
          ~ mean(
            .x,
            na.rm = TRUE
          )
        ),
        .groups = "drop"
      )
    
  } else {
    
    data.frame()
  }
  
  # ============================================================
  # 33. Return results
  # ============================================================
  
  list(
    
    summary = final_summary,
    
    capture = capture_summary,
    
    capture_site =
      capture_site_summary,
    
    psi_list =
      psi_list[keep],
    
    capture_list =
      capture_list[keep],
    
    lambda_list =
      lambda_list[keep],
    
    p_detect_list =
      p_detect_list[keep],
    
    occupancy_models =
      occupancy_models[keep],
    
    capture_models =
      capture_models[keep],
    
    abundance_models =
      abundance_models[keep],
    
    reduced_data =
      reduced_data,
    
    sample_data =
      sample_data,
    
    long_df =
      long_df,
    
    lv_sites =
      lv_sites_combined,
    
    lv_species =
      lv_species_combined,
    
    mean_lv_sites =
      mean_lv_sites,
    
    mean_lv_species =
      mean_lv_species,
    
    filter_summary = list(
      otu_stats =
        otu_stats,
      
      kept_otus =
        keep_otus,
      
      min_species_sum =
        min_species_sum,
      
      min_detection_replicates =
        min_detection_replicates,
      
      abundance_threshold =
        abundance_threshold
    ),
    
    diagnostic_AIC =
      diagnostic_AIC,
    
    gllvm_diagnostics =
      gllvm_diagnostics,
    
    gllvm_separation =
      gllvm_separation,
    
    gllvm_excluded_otus =
      gllvm_excluded_otus,
    
    gllvm_fit_warnings =
      gllvm_fit_warnings,
    
    # Maps retained list positions back to original iteration numbers.
    occupancy_model_iterations =
      keep,
    
    
    gllvm_settings = list(
      num_lv_c = num_lv_c,
      num_lv = num_lv,
      n_init = gllvm_n_init,
      gradient_tol = gllvm_gradient_tol,
      strict = gllvm_strict,
      condition_tol = gllvm_condition_tol,
      eigen_tol = gllvm_eigen_tol,
      warn_separation = gllvm_warn_separation,
      auto_exclude = gllvm_auto_exclude,
      min_occupied_sites = gllvm_min_occupied_sites,
      excluded_prior = gllvm_excluded_prior
    ),
    
    iteration_success =
      iteration_success,
    
    retained_iterations =
      keep,
    
    abundance_formula =
      abundance_formula,
    
    capture_formula =
      capture_formula,
    
    note = paste(
      "Three-level GLLVM-GLMM model.",
      "Hierarchy: site occupancy Z estimated jointly across taxa using gllvm,",
      "biological-sample capture A estimated using a binomial GLMM,",
      "and replicate-level abundance Y estimated using a count GLMM."
    )
  )
}
                         
    
   
