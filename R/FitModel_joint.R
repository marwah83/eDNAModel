#' Fit a Joint Hierarchical eDNA Occupancy--Capture--Abundance Model
#'
#' @description
#' Fits a joint hierarchical model for environmental DNA (eDNA) metabarcoding
#' count data using Template Model Builder (TMB). The model represents the
#' observation process through linked occupancy, capture/detection, and
#' sequencing-abundance components.
#'
#' The model can be fitted with Poisson, negative-binomial, zero-inflated
#' Poisson (ZIP), or zero-inflated negative-binomial (ZINB) abundance
#' distributions. OTU-level random effects can be included separately in the
#' occupancy, capture, abundance, and zero-inflation components. Additional
#' sample and sample-by-OTU random effects can also be included.
#'
#' Latent occupancy and biological-sample capture states are analytically
#' marginalized in the observed-data likelihood. Enabled Gaussian random
#' effects are integrated using the Laplace approximation implemented in TMB.
#'
#' @details
#' The hierarchical model contains three principal ecological/observation
#' components.
#'
#' **Occupancy model**
#'
#' For site \eqn{i} and OTU \eqn{k}, occupancy probability is modeled as
#'
#' \deqn{
#' \mathrm{logit}(\psi_{ik})
#' =
#' X^{(\psi)}_{ik}\beta_{\psi}
#' +
#' b^{(\psi)}_k,
#' }
#'
#' where \eqn{X^{(\psi)}} is the occupancy design matrix,
#' \eqn{\beta_{\psi}} contains fixed effects, and
#' \eqn{b^{(\psi)}_k} is an optional OTU-level random intercept.
#'
#' **Capture/detection model**
#'
#' Conditional on occupancy, the probability that OTU \eqn{k} is captured in
#' biological sample \eqn{j} is modeled as
#'
#' \deqn{
#' \mathrm{logit}(p_{jk})
#' =
#' X^{(p)}_{jk}\beta_p
#' +
#' b^{(p)}_k.
#' }
#'
#' **Abundance model**
#'
#' Conditional on capture, expected read abundance is modeled as
#'
#' \deqn{
#' \log(\lambda_{rjk})
#' =
#' X^{(\lambda)}_{rjk}\beta_{\lambda}
#' +
#' b^{(\lambda)}_k
#' +
#' b^{(s)}_j
#' +
#' b^{(s \times k)}_{jk}
#' +
#' o_{rjk},
#' }
#'
#' where \eqn{o_{rjk}} is an optional log-scale abundance offset.
#'
#' For ZIP and ZINB models, additional structural zeros are represented using
#' a zero-inflation probability
#'
#' \deqn{
#' \mathrm{logit}(\pi_k)
#' =
#' \alpha_{\pi}
#' +
#' b^{(\pi)}_k,
#' }
#'
#' where \eqn{\alpha_{\pi}} is the population-level zero-inflation intercept
#' and \eqn{b^{(\pi)}_k} is an optional OTU-level random effect. This allows
#' structural-zero probabilities to vary among OTUs rather than assuming a
#' common dropout probability for all taxa.
#'
#'
#' @section Numerical convergence:
#'
#' Numerical convergence is assessed using several complementary diagnostics
#' rather than relying solely on the optimizer return code.
#'
#' A fit is classified as strictly numerically converged only when all of the
#' following conditions are satisfied:
#'
#' \itemize{
#'   \item \code{nlminb()} returns convergence code 0;
#'   \item all components of the final outer TMB gradient are finite;
#'   \item the maximum absolute outer gradient is no greater than
#'         \code{gradient_tol};
#'   \item \code{TMB::sdreport()} completes successfully;
#'   \item the Hessian is positive definite according to
#'         \code{sdreport()$pdHess}; and
#'   \item all standard errors returned for the fixed parameter vector are
#'         finite.
#' }
#'
#' The default strict gradient tolerance is \eqn{10^{-3}}. The gradient
#' criterion is
#'
#' \deqn{
#' \max_j
#' \left|
#' \frac{\partial \tilde{\ell}(\theta)}
#'      {\partial \theta_j}
#' \right|
#' \leq
#' \mathrm{gradient\_tol},
#' }
#'
#' where \eqn{\tilde{\ell}(\theta)} denotes the TMB objective after Laplace
#' approximation over enabled random effects.
#'
#' Because the outer TMB gradient can exhibit small numerical fluctuations near
#' an otherwise stable optimum, the function distinguishes strict convergence
#' from a user-facing marginal-gradient region. With the defaults,
#' \code{gradient_tol = 1e-3} and \code{gradient_marginal_factor = 5},
#' gradient diagnostics are classified as:
#'
#' \itemize{
#'   \item \code{"PASS"} when
#'         \eqn{\max |\nabla\tilde{\ell}| \le 0.001};
#'   \item \code{"MARGINAL"} when
#'         \eqn{0.001 < \max |\nabla\tilde{\ell}| \le 0.005}; and
#'   \item \code{"FAIL"} when
#'         \eqn{\max |\nabla\tilde{\ell}| > 0.005}.
#' }
#'
#' The marginal classification does not redefine the strict mathematical
#' gradient criterion. Instead, it distinguishes small residual gradients near
#' an otherwise numerically stable optimum from more substantial optimization
#' problems. A fit with optimizer code 0, successful \code{sdreport()},
#' positive-definite Hessian, finite fixed-parameter standard errors, and a
#' marginal gradient is therefore reported as \code{"MARGINAL"} rather than
#' automatically being labelled a numerical failure.
#'
#'
#' @section Optimization restarts and solution selection:
#'
#' If the strict optimizer/gradient criteria are not met after the first
#' optimization pass, the function can restart \code{nlminb()} from the
#' current solution up to \code{max_restarts} times.
#'
#' Restarted optimizations can return essentially identical likelihood values
#' while their evaluated outer gradients fluctuate slightly because of
#' numerical noise in the Laplace-approximated objective. Consequently, the
#' function does not automatically retain the final optimization pass.
#'
#' After all required passes have been completed, the minimum objective value
#' is identified. Optimization passes whose objective values differ from this
#' minimum by no more than
#'
#' \deqn{
#' \mathrm{objective\_rel\_tol}
#' \times
#' \max(1, |\tilde{\ell}_a|, |\tilde{\ell}_b|)
#' }
#'
#' are treated as effectively equivalent likelihood solutions. Among these
#' equivalent solutions, the pass with the smallest finite maximum absolute
#' outer gradient is retained.
#'
#' This prevents a later restart with an effectively identical likelihood but
#' a slightly larger numerical gradient from replacing an otherwise better
#' solution.
#'
#' Restarting may terminate early when consecutive optimization passes have
#' effectively equivalent objective values and the improvement in maximum
#' absolute gradient is no greater than \code{gradient_improvement_tol}.
#' This provides a practical stall criterion for cases in which repeated
#' optimization reaches the same likelihood while the gradient fluctuates
#' numerically.
#'
#'
#' @section Hessian diagnostic:
#'
#' Positive definiteness of the Hessian is assessed using the
#' \code{pdHess} diagnostic returned by \code{TMB::sdreport()}. A
#' non-positive-definite Hessian may indicate a flat likelihood direction,
#' weak or non-identification of one or more parameters, a saddle point, or
#' other numerical difficulties.
#'
#' This diagnostic evaluates positive definiteness rather than merely positive
#' semidefiniteness. It is considered jointly with the optimizer status and
#' gradient diagnostic; \code{pdHess = TRUE} alone is not interpreted as
#' evidence of convergence.
#'
#'
#' @section Standard-error diagnostic:
#'
#' Fixed-parameter standard errors are obtained from
#' \code{summary(sdreport_object, "fixed")}. They pass the finite-standard-error
#' diagnostic only when every returned standard error satisfies
#' \code{is.finite()}. Consequently, \code{NA}, \code{NaN}, \code{Inf}, and
#' \code{-Inf} standard errors cause this diagnostic to fail.
#'
#' This is a numerical sanity check and does not imply that finite standard
#' errors are necessarily small or scientifically informative. Large but
#' finite standard errors are handled separately as heuristic parameter
#' warnings.
#'
#'
#' @section Heuristic parameter diagnostics:
#'
#' In addition to the formal numerical convergence criteria, the function can
#' flag potentially problematic parameter estimates. These diagnostics are
#' intended to assist interpretation and do not by themselves cause a fit to
#' be classified as non-converged.
#'
#' Two heuristic checks are currently used:
#'
#' \itemize{
#'   \item estimated log-standard-deviation parameters below
#'         \code{log_sd_warning_threshold}, which may indicate a variance
#'         component close to zero; and
#'   \item unusually large standard errors relative to the absolute parameter
#'         estimate, based on \code{se_estimate_ratio_threshold}.
#' }
#'
#' These thresholds are diagnostic guidelines rather than universal
#' statistical criteria and should be interpreted in the context of the fitted
#' model and data.
#'
#'
#' @section Memory use and sdreport:
#'
#' For large eDNA OTU matrices, covariance and joint-precision calculations
#' performed by \code{TMB::sdreport()} can require substantial memory.
#' Therefore, both \code{get_report_covariance = FALSE} and
#' \code{get_joint_precision = FALSE} are memory-efficient defaults.
#'
#' Users requiring the covariance matrix of reported quantities can explicitly
#' set \code{get_report_covariance = TRUE}. Users requiring the joint precision
#' matrix of fixed and random effects can set
#' \code{get_joint_precision = TRUE}. Either option may substantially increase
#' memory requirements and computation time for large models.
#'
#'
#' @section Verbose output:
#'
#' Package-level progress messages and the internal TMB optimizer trace are
#' controlled separately.
#'
#' Setting \code{verbose = TRUE} prints useful eDNAModel progress information,
#' including optimization-pass summaries, selected-pass information,
#' \code{sdreport()} progress, and final numerical diagnostics.
#'
#' The lower-level TMB optimization trace (for example, repeated
#' \code{iter:} and \code{mgc:} messages) is controlled by
#' \code{tmb_verbose}. Its default is \code{FALSE}, preventing potentially
#' thousands of TMB inner-optimizer messages from being printed during normal
#' use. Advanced users can set \code{tmb_verbose = TRUE} when detailed TMB
#' tracing is required.
#'
#'
#' @param phyloseq A \code{phyloseq} object containing the OTU count table and
#'   associated sample metadata.
#'
#' @param site_col Character string giving the sample-data column identifying
#'   sampling sites.
#'
#' @param sample_col Character string giving the column identifying biological
#'   samples. Default is \code{"Name"}.
#'
#' @param replicate_col Optional character string identifying technical or PCR
#'   replicates. Default is \code{NULL}.
#'
#' @param otu_col Character string identifying the OTU column in the long-format
#'   data. Default is \code{"OTU"}.
#'
#' @param count_col Character string identifying the observed read-count
#'   column. Default is \code{"y"}.
#'
#' @param occupancy_formula A model formula specifying fixed effects for the
#'   occupancy component. Default is \code{~ 1}.
#'
#' @param capture_formula A model formula specifying fixed effects for the
#'   capture/detection component. Default is \code{~ 1}.
#'
#' @param abundance_formula A model formula specifying fixed effects for the
#'   abundance component. Default is \code{~ 1}.
#'
#' @param abundance_offset Optional character string identifying a positive
#'   abundance-exposure variable. The logarithm of this variable is included
#'   as an offset. Default is \code{NULL}.
#'
#' @param abundance_family Character string specifying the abundance
#'   distribution. One of \code{"poisson"}, \code{"nbinom"}, \code{"zip"}, or
#'   \code{"zinb"}.
#'
#' @param min_species_sum Minimum total read count required for an OTU to be
#'   retained. Default is 10.
#'
#' @param min_detection_replicates Minimum number of positive observations
#'   required for an OTU to be retained. Default is 1.
#'
#' @param random_occ_otu Logical. Include an OTU-level random intercept in the
#'   occupancy model. Default is \code{TRUE}.
#'
#' @param random_capture_otu Logical. Include an OTU-level random intercept in
#'   the capture model. Default is \code{TRUE}.
#'
#' @param random_abund_otu Logical. Include an OTU-level random intercept in the
#'   abundance model. Default is \code{TRUE}.
#'
#' @param random_sample Logical. Include a biological-sample random intercept
#'   in the abundance model. Default is \code{TRUE}.
#'
#' @param random_sample_otu Logical. Include a sample-by-OTU random intercept
#'   in the abundance model. Default is \code{FALSE}.
#'
#' @param random_zi_otu Logical. Include an OTU-level random intercept in the
#'   zero-inflation model. This option is used only for ZIP and ZINB models.
#'   Default is \code{TRUE}.
#'
#' @param gradient_tol Positive numeric value specifying the maximum acceptable
#'   absolute component of the final outer TMB gradient for strict numerical
#'   convergence. Default is \code{1e-3}.
#'
#' @param gradient_marginal_factor Numeric value greater than 1 defining the
#'   user-facing marginal-gradient region as
#'   \code{gradient_marginal_factor * gradient_tol}. The default is 5, giving
#'   a marginal upper limit of 0.005 when \code{gradient_tol = 0.001}.
#'   This diagnostic classification does not alter the strict gradient
#'   convergence criterion.
#'
#' @param max_restarts Non-negative integer giving the maximum number of
#'   automatic optimizer restarts after the initial optimization pass.
#'   Default is 2.
#'
#' @param iter_max Maximum number of \code{nlminb()} iterations per
#'   optimization pass. Default is 5000.
#'
#' @param eval_max Maximum number of objective evaluations per optimization
#'   pass. Default is 10000.
#'
#' @param rel_tol Relative convergence tolerance supplied to \code{nlminb()}.
#'   Default is \code{1e-10}.
#'
#' @param objective_rel_tol Relative tolerance used to determine whether
#'   objective values from different optimization passes are effectively
#'   equivalent. Two objective values \eqn{f_a} and \eqn{f_b} are treated as
#'   equivalent when
#'   \code{abs(f_a - f_b) <= objective_rel_tol *
#'   max(1, abs(f_a), abs(f_b))}. Default is \code{1e-8}.
#'
#' @param gradient_improvement_tol Non-negative numeric value specifying the
#'   minimum reduction in maximum absolute gradient considered a meaningful
#'   improvement between effectively equivalent optimization passes.
#'   Default is \code{1e-4}.
#'
#' @param n_gradient_report Number of parameters with the largest absolute
#'   gradients to report when the strict gradient criterion is not satisfied.
#'   Default is 2.
#'
#' @param get_report_covariance Logical. Passed to
#'   \code{TMB::sdreport(getReportCovariance = ...)}. The default is
#'   \code{FALSE} to reduce memory use for large OTU datasets.
#'
#' @param get_joint_precision Logical. Request the joint precision matrix from
#'   \code{TMB::sdreport()} when random effects are present. The default is
#'   \code{FALSE} to reduce memory use and computation time. Set to
#'   \code{TRUE} only when the joint precision matrix is required.
#'
#' @param log_sd_warning_threshold Numeric threshold used to flag very small
#'   estimated random-effect standard deviations on the log-SD scale. This is
#'   a heuristic diagnostic only. Default is -5.
#'
#' @param se_estimate_ratio_threshold Numeric threshold used to flag fixed
#'   parameters whose standard error is large relative to the absolute
#'   estimate. This is a heuristic diagnostic and does not determine
#'   convergence. Default is 10.
#'
#' @param estimate_zero_tol Positive numeric tolerance below which an estimate
#'   is treated as effectively zero when computing SE-to-estimate ratios.
#'   Default is \code{1e-6}.
#'
#' @param DLL Character string giving the name of the compiled TMB dynamic
#'   library. Default is \code{"eDNAModel"}.
#'
#' @param verbose Logical. If \code{TRUE}, print package-level fitting progress,
#'   optimization-pass summaries, \code{sdreport()} progress, and numerical
#'   diagnostics. Default is \code{TRUE}.
#'
#' @param tmb_verbose Logical. If \code{TRUE}, allow the internal TMB
#'   optimization trace to be printed. Default is \code{FALSE}. This is
#'   separate from \code{verbose} so normal package progress can be displayed
#'   without printing thousands of TMB inner-optimization messages.
#'
#'
#' @return
#' A list containing the fitted optimizer object, TMB objective,
#' \code{sdreport} results, parameter summaries, processed model data, and
#' numerical diagnostics. Important components include:
#'
#' \describe{
#'
#'   \item{\code{fit}}{
#'     The selected \code{nlminb()} optimization result.
#'   }
#'
#'   \item{\code{tmb_object}}{
#'     The fitted TMB objective object.
#'   }
#'
#'   \item{\code{sdreport}}{
#'     The \code{TMB::sdreport()} object, when successful.
#'   }
#'
#'   \item{\code{fixed_effects}}{
#'     Estimates and standard errors for the fixed parameter vector.
#'   }
#'
#'   \item{\code{derived}}{
#'     Reported or derived parameter summaries.
#'   }
#'
#'   \item{\code{convergence}}{
#'     Detailed optimizer, gradient, Hessian, standard-error, restart,
#'     selected-pass, and overall numerical-convergence diagnostics.
#'   }
#'
#'   \item{\code{diagnostics}}{
#'     User-facing diagnostic table and heuristic parameter warnings.
#'   }
#'
#'   \item{\code{sdreport_settings}}{
#'     Settings used when calling \code{TMB::sdreport()}.
#'   }
#' }
#'
#'
#' @note
#' The convergence diagnostics assess numerical optimization and local
#' curvature of the Laplace-approximated objective. They do not, by themselves,
#' establish the accuracy of the Laplace approximation or the biological
#' adequacy of the fitted model.
#'
#' A \code{"MARGINAL"} gradient classification should be interpreted together
#' with objective stability across optimization passes, optimizer status,
#' Hessian positive definiteness, and standard-error diagnostics. It is
#' intended to distinguish small residual numerical gradients from clear
#' optimization failures rather than to redefine the strict gradient
#' convergence criterion.
#'
#' Simulation-based parameter-recovery and coverage studies are recommended
#' when evaluating approximation accuracy, particularly for sparse and highly
#' discrete eDNA datasets.
#'
#'
#' @examples
#' \dontrun{
#'
#' fit_zinb <- FitModel_joint(
#'
#'     phyloseq = ps,
#'
#'     site_col = "Sampling.area.Name",
#'     sample_col = "Name",
#'     replicate_col = "Replicate",
#'
#'     otu_col = "OTU",
#'     count_col = "y",
#'
#'     occupancy_formula = ~ 1,
#'     capture_formula = ~ 1,
#'     abundance_formula = ~ 1,
#'
#'     abundance_family = "zinb",
#'
#'     random_occ_otu = TRUE,
#'     random_capture_otu = TRUE,
#'     random_abund_otu = TRUE,
#'     random_zi_otu = TRUE,
#'
#'     random_sample = TRUE,
#'     random_sample_otu = FALSE,
#'
#'     gradient_tol = 1e-3,
#'     gradient_marginal_factor = 5,
#'
#'     max_restarts = 2,
#'
#'     objective_rel_tol = 1e-8,
#'     gradient_improvement_tol = 1e-4,
#'
#'     get_report_covariance = FALSE,
#'     get_joint_precision = FALSE,
#'
#'     verbose = TRUE,
#'     tmb_verbose = FALSE
#' )
#'
#'
#' # Parameter estimates
#' fit_zinb$fixed_effects
#'
#' # Overall convergence classification
#' fit_zinb$convergence$overall_status
#'
#' # Strict convergence
#' fit_zinb$convergence$strict_convergence
#'
#' # Acceptable PASS/MARGINAL convergence
#' fit_zinb$convergence$acceptable_convergence
#'
#' # Maximum absolute outer gradient
#' fit_zinb$convergence$max_abs_gradient
#'
#' # Gradient classification
#' fit_zinb$convergence$gradient_status
#'
#' # Hessian diagnostic
#' fit_zinb$convergence$pd_hessian
#'
#' # Standard-error diagnostic
#' fit_zinb$convergence$finite_standard_errors
#'
#' # Optimization passes and selected pass
#' fit_zinb$convergence$n_optimization_passes
#' fit_zinb$convergence$selected_pass
#'
#' # Complete user-facing diagnostic table
#' fit_zinb$diagnostics$table
#' }
#'
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats model.matrix nlminb qlogis
#' @importFrom dplyr group_by summarise filter select pull first across
#'   all_of left_join
#'
#' @export
FitModel_joint <- function(
    phyloseq,
    site_col,
    sample_col = "Name",
    replicate_col = NULL,
    otu_col = "OTU",
    count_col = "y",
    
    occupancy_formula = ~ 1,
    capture_formula = ~ 1,
    abundance_formula = ~ 1,
    abundance_offset = NULL,
    
    abundance_family = c(
      "poisson",
      "nbinom",
      "zip",
      "zinb"
    ),
    
    min_species_sum = 10,
    min_detection_replicates = 1,
    
    # ==========================================================
    # BUILT-IN RANDOM EFFECTS
    # ==========================================================
    
    random_occ_otu = TRUE,
    random_capture_otu = TRUE,
    random_abund_otu = TRUE,
    
    random_sample = TRUE,
    random_sample_otu = FALSE,
    
    random_zi_otu = TRUE,
    
    # ==========================================================
    # NUMERICAL CONVERGENCE
    # ==========================================================
    
    # Strict gradient criterion
    gradient_tol = 1e-3,
    
    # PASS:     gradient <= 0.001
    # MARGINAL: gradient <= 0.005
    # FAIL:     gradient > 0.005
    gradient_marginal_factor = 5,
    
    max_restarts = 2L,
    
    iter_max = 5000L,
    eval_max = 10000L,
    rel_tol = 1e-10,
    
    # Two optimization passes are regarded as having effectively
    # equivalent objectives when their difference is smaller than
    # this relative tolerance.
    objective_rel_tol = 1e-8,
    
    # Improvement in maximum gradient smaller than this is not
    # considered practically meaningful.
    gradient_improvement_tol = 1e-4,
    
    n_gradient_report = 2L,
    
    # ==========================================================
    # SDREPORT MEMORY CONTROLS
    # ==========================================================
    
    get_report_covariance = FALSE,
    
    # FALSE is safer for large OTU matrices.
    get_joint_precision = FALSE,
    
    # ==========================================================
    # HEURISTIC DIAGNOSTICS
    # ==========================================================
    
    log_sd_warning_threshold = -5,
    se_estimate_ratio_threshold = 10,
    estimate_zero_tol = 1e-6,
    
    DLL = "eDNAModel",
    
    # Package-level messages
    verbose = TRUE,
    
    # TMB inner optimiser output
    tmb_verbose = FALSE
) {
  
  # ==========================================================
  # 1. MATCH FAMILY
  # ==========================================================
  
  abundance_family <- match.arg(
    abundance_family
  )
  
  
  # ==========================================================
  # 2. VALIDATE NUMERICAL CONTROLS
  # ==========================================================
  
  if (
    length(gradient_tol) != 1L ||
    !is.finite(gradient_tol) ||
    gradient_tol <= 0
  ) {
    stop(
      "gradient_tol must be a single positive finite number."
    )
  }
  
  if (
    length(gradient_marginal_factor) != 1L ||
    !is.finite(gradient_marginal_factor) ||
    gradient_marginal_factor <= 1
  ) {
    stop(
      "gradient_marginal_factor must be greater than 1."
    )
  }
  
  if (
    length(objective_rel_tol) != 1L ||
    !is.finite(objective_rel_tol) ||
    objective_rel_tol <= 0
  ) {
    stop(
      "objective_rel_tol must be a positive finite number."
    )
  }
  
  if (
    length(gradient_improvement_tol) != 1L ||
    !is.finite(gradient_improvement_tol) ||
    gradient_improvement_tol < 0
  ) {
    stop(
      "gradient_improvement_tol must be a non-negative finite number."
    )
  }
  
  max_restarts <- as.integer(max_restarts)
  iter_max <- as.integer(iter_max)
  eval_max <- as.integer(eval_max)
  n_gradient_report <- as.integer(n_gradient_report)
  
  if (
    length(max_restarts) != 1L ||
    is.na(max_restarts) ||
    max_restarts < 0L
  ) {
    stop(
      "max_restarts must be a non-negative integer."
    )
  }
  
  if (
    length(iter_max) != 1L ||
    is.na(iter_max) ||
    iter_max < 1L
  ) {
    stop(
      "iter_max must be a positive integer."
    )
  }
  
  if (
    length(eval_max) != 1L ||
    is.na(eval_max) ||
    eval_max < 1L
  ) {
    stop(
      "eval_max must be a positive integer."
    )
  }
  
  if (
    length(n_gradient_report) != 1L ||
    is.na(n_gradient_report) ||
    n_gradient_report < 1L
  ) {
    stop(
      "n_gradient_report must be a positive integer."
    )
  }
  
  
  # ==========================================================
  # 3. Convert the phyloseq object to long format
  #
  # prepare_long_data() is assumed to be an internal
  # eDNAModel helper already defined in the package.
  # ==========================================================
  
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col,
    nested_cols = unique(
      stats::na.omit(
        c(sample_col, replicate_col)
      )
    )
  )
  
  dat <- prep$long_df
  
  
  # ==========================================================
  # 4. Check required variables
  # ==========================================================
  
  required <- unique(
    c(
      site_col,
      sample_col,
      otu_col,
      count_col
    )
  )
  
  if (!is.null(replicate_col)) {
    required <- unique(
      c(required, replicate_col)
    )
  }
  
  missing_cols <- setdiff(
    required,
    names(dat)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  
  # ==========================================================
  # 5. Standardize basic variable types
  # ==========================================================
  
  dat[[site_col]] <- as.character(dat[[site_col]])
  dat[[sample_col]] <- as.character(dat[[sample_col]])
  dat[[otu_col]] <- as.character(dat[[otu_col]])
  dat[[count_col]] <- as.numeric(dat[[count_col]])
  
  if (!is.null(replicate_col)) {
    dat[[replicate_col]] <-
      as.character(dat[[replicate_col]])
  }
  
  
  # ==========================================================
  # 6. Remove observations with invalid counts
  # ==========================================================
  
  dat <- dat[
    !is.na(dat[[count_col]]) &
      is.finite(dat[[count_col]]),
    ,
    drop = FALSE
  ]
  
  if (any(dat[[count_col]] < 0)) {
    stop("Counts must be non-negative.")
  }
  
  
  # ==========================================================
  # 7. Filter rare OTUs
  # ==========================================================
  
  otu_stats <- dat |>
    dplyr::group_by(.data[[otu_col]]) |>
    dplyr::summarise(
      total_count =
        sum(.data[[count_col]], na.rm = TRUE),
      
      detected_replicates =
        sum(.data[[count_col]] > 0, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  retained_otus <- otu_stats |>
    dplyr::filter(
      .data$total_count >= min_species_sum,
      .data$detected_replicates >= min_detection_replicates
    ) |>
    dplyr::pull(.data[[otu_col]])
  
  dat <- dat |>
    dplyr::filter(
      .data[[otu_col]] %in% retained_otus
    )
  
  if (nrow(dat) == 0L) {
    stop("No OTUs remain after filtering.")
  }
  
  
  # ==========================================================
  # 8. Construct hierarchical grouping factors
  # ==========================================================
  
  dat$.site <- factor(dat[[site_col]])
  dat$.otu <- factor(dat[[otu_col]])
  dat$.sample <- factor(dat[[sample_col]])
  
  # One latent occupancy state for each site x OTU.
  dat$.site_otu <- interaction(
    dat$.site,
    dat$.otu,
    drop = TRUE,
    lex.order = TRUE
  )
  
  # One latent capture state for each biological sample x OTU.
  dat$.sample_otu <- interaction(
    dat$.site,
    dat$.sample,
    dat$.otu,
    drop = TRUE,
    lex.order = TRUE
  )
  
  
  # ==========================================================
  # 9. Construct sample x OTU data
  # ==========================================================
  
  sample_df <- dat |>
    dplyr::group_by(.data$.sample_otu) |>
    dplyr::summarise(
      .site_otu = dplyr::first(.data$.site_otu),
      .otu = dplyr::first(.data$.otu),
      .sample = dplyr::first(.data$.sample),
      
      sample_positive =
        as.integer(any(.data[[count_col]] > 0)),
      
      .groups = "drop"
    )
  
  
  # ==========================================================
  # 10. Add capture covariates
  # ==========================================================
  
  cap_vars <- intersect(
    all.vars(capture_formula),
    names(dat)
  )
  
  if (length(cap_vars) > 0L) {
    
    cap_covars <- dat |>
      dplyr::select(
        .data$.sample_otu,
        dplyr::all_of(cap_vars)
      ) |>
      dplyr::group_by(.data$.sample_otu) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(cap_vars),
          ~ dplyr::first(.x)
        ),
        .groups = "drop"
      )
    
    sample_df <- dplyr::left_join(
      sample_df,
      cap_covars,
      by = ".sample_otu"
    )
  }
  
  
  # ==========================================================
  # 11. Construct site x OTU data
  # ==========================================================
  
  site_df <- dat |>
    dplyr::group_by(.data$.site_otu) |>
    dplyr::summarise(
      .otu = dplyr::first(.data$.otu),
      
      site_positive =
        as.integer(any(.data[[count_col]] > 0)),
      
      .groups = "drop"
    )
  
  
  # ==========================================================
  # 12. Add occupancy covariates
  # ==========================================================
  
  occ_vars <- intersect(
    all.vars(occupancy_formula),
    names(dat)
  )
  
  if (length(occ_vars) > 0L) {
    
    occ_covars <- dat |>
      dplyr::select(
        .data$.site_otu,
        dplyr::all_of(occ_vars)
      ) |>
      dplyr::group_by(.data$.site_otu) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(occ_vars),
          ~ dplyr::first(.x)
        ),
        .groups = "drop"
      )
    
    site_df <- dplyr::left_join(
      site_df,
      occ_covars,
      by = ".site_otu"
    )
  }
  
  
  # ==========================================================
  # 13. Construct fixed-effect design matrices
  # ==========================================================
  
  X_occ <- stats::model.matrix(
    occupancy_formula,
    data = site_df
  )
  
  X_cap <- stats::model.matrix(
    capture_formula,
    data = sample_df
  )
  
  X_abund <- stats::model.matrix(
    abundance_formula,
    data = dat
  )
  
  
  # ==========================================================
  # 14. Construct abundance offset
  #
  # The user supplies the exposure itself. We put it on
  # the log scale before passing it to the C++ model.
  # ==========================================================
  
  if (is.null(abundance_offset)) {
    
    offset_abund <- rep(0, nrow(dat))
    
  } else {
    
    if (!(abundance_offset %in% names(dat))) {
      stop(
        "abundance_offset '",
        abundance_offset,
        "' was not found."
      )
    }
    
    off <- as.numeric(
      dat[[abundance_offset]]
    )
    
    if (any(!is.finite(off) | off <= 0)) {
      stop(
        "abundance_offset must contain positive finite values."
      )
    }
    
    offset_abund <- log(off)
  }
  
  
  # ==========================================================
  # 15. Create mappings between hierarchy levels
  #
  # C++/TMB indices are zero based.
  # ==========================================================
  
  site_group_levels <- levels(dat$.site_otu)
  sample_group_levels <- levels(dat$.sample_otu)
  
  row_sample_group <- match(
    as.character(dat$.sample_otu),
    sample_group_levels
  ) - 1L
  
  sample_site_group <- match(
    as.character(sample_df$.site_otu),
    site_group_levels
  ) - 1L
  
  occ_otu <- as.integer(site_df$.otu) - 1L
  sample_otu <- as.integer(sample_df$.otu) - 1L
  row_otu <- as.integer(dat$.otu) - 1L
  row_sample_id <- as.integer(dat$.sample) - 1L
  
  # The sample x OTU abundance random effect uses the
  # same grouping as the sample x OTU capture state.
  row_sample_otu_re <- row_sample_group
  
  if (anyNA(row_sample_group)) {
    stop("Invalid row_sample_group mapping.")
  }
  
  if (anyNA(sample_site_group)) {
    stop("Invalid sample_site_group mapping.")
  }
  
  
  # ==========================================================
  # 16. Encode abundance distribution for C++
  # ==========================================================
  
  family_code <- switch(
    abundance_family,
    poisson = 0L,
    nbinom = 1L,
    zip = 2L,
    zinb = 3L
  )
  
  # OTU-specific zero inflation exists only for ZIP/ZINB.
  use_zi_otu <-
    isTRUE(random_zi_otu) &&
    abundance_family %in% c("zip", "zinb")
  
  
  # ==========================================================
  # 17. Construct TMB data list
  # ==========================================================
  
  data_tmb <- list(
    y = dat[[count_col]],
    
    X_occ = X_occ,
    X_cap = X_cap,
    X_abund = X_abund,
    
    offset_abund = offset_abund,
    
    occ_otu = as.integer(occ_otu),
    sample_otu = as.integer(sample_otu),
    
    row_sample_group =
      as.integer(row_sample_group),
    
    sample_site_group =
      as.integer(sample_site_group),
    
    row_otu =
      as.integer(row_otu),
    
    row_sample_id =
      as.integer(row_sample_id),
    
    row_sample_otu_re =
      as.integer(row_sample_otu_re),
    
    sample_positive =
      as.integer(sample_df$sample_positive),
    
    site_positive =
      as.integer(site_df$site_positive),
    
    n_site_groups =
      nrow(site_df),
    
    n_sample_groups =
      nrow(sample_df),
    
    family_code =
      family_code,
    
    use_occ_otu =
      as.integer(random_occ_otu),
    
    use_cap_otu =
      as.integer(random_capture_otu),
    
    use_abund_otu =
      as.integer(random_abund_otu),
    
    use_sample_re =
      as.integer(random_sample),
    
    use_sample_otu_re =
      as.integer(random_sample_otu),
    
    use_zi_otu =
      as.integer(use_zi_otu)
  )
  
  
  # ==========================================================
  # 18. Dimensions
  # ==========================================================
  
  n_otu <- nlevels(dat$.otu)
  n_sample <- nlevels(dat$.sample)
  n_sample_otu <- nrow(sample_df)
  
  
  # ==========================================================
  # 19. Generate sensible starting values
  # ==========================================================
  
  # Keep initial probabilities away from exactly 0 and 1.
  site_positive_rate <- mean(
    site_df$site_positive
  )
  
  site_positive_rate <- pmin(
    pmax(site_positive_rate, 0.05),
    0.95
  )
  
  sample_positive_rate <- mean(
    sample_df$sample_positive
  )
  
  sample_positive_rate <- pmin(
    pmax(sample_positive_rate, 0.05),
    0.95
  )
  
  positive_counts <- dat[[count_col]][
    dat[[count_col]] > 0
  ]
  
  mean_count_start <- if (length(positive_counts) > 0L) {
    mean(positive_counts)
  } else {
    1
  }
  
  beta_occ_start <- rep(
    0,
    ncol(X_occ)
  )
  
  beta_occ_start[1] <- stats::qlogis(
    site_positive_rate
  )
  
  beta_cap_start <- rep(
    0,
    ncol(X_cap)
  )
  
  beta_cap_start[1] <- stats::qlogis(
    sample_positive_rate
  )
  
  beta_abund_start <- rep(
    0,
    ncol(X_abund)
  )
  
  beta_abund_start[1] <- log(
    pmax(mean_count_start, 1e-4)
  )
  
  
  # ==========================================================
  # 20. Define TMB parameter list
  # ==========================================================
  
  parameters <- list(
    
    # Fixed effects
    beta_occ = beta_occ_start,
    beta_cap = beta_cap_start,
    beta_abund = beta_abund_start,
    
    # OTU random effects
    b_occ_otu = rep(0, n_otu),
    b_cap_otu = rep(0, n_otu),
    b_abund_otu = rep(0, n_otu),
    
    # Sample random effects
    b_sample = rep(0, n_sample),
    b_sample_otu = rep(0, n_sample_otu),
    
    # OTU-specific zero-inflation deviations
    b_zi_otu = rep(0, n_otu),
    
    # Random-effect log standard deviations
    log_sd_occ_otu = log(0.5),
    log_sd_cap_otu = log(0.5),
    log_sd_abund_otu = log(0.5),
    log_sd_sample = log(0.5),
    log_sd_sample_otu = log(0.5),
    log_sd_zi_otu = log(0.5),
    
    # Negative-binomial dispersion
    log_theta = log(10),
    
    # Population-level zero-inflation intercept
    zi_intercept = stats::qlogis(0.05)
  )
  
  
  # ==========================================================
  # 21. Parameter mapping
  #
  # Parameters that are not used by a particular model are
  # mapped to NA so that TMB does not estimate them.
  # ==========================================================
  
  map <- list()
  
  if (!random_occ_otu) {
    map$b_occ_otu <- factor(rep(NA, n_otu))
    map$log_sd_occ_otu <- factor(NA)
  }
  
  if (!random_capture_otu) {
    map$b_cap_otu <- factor(rep(NA, n_otu))
    map$log_sd_cap_otu <- factor(NA)
  }
  
  if (!random_abund_otu) {
    map$b_abund_otu <- factor(rep(NA, n_otu))
    map$log_sd_abund_otu <- factor(NA)
  }
  
  if (!random_sample) {
    map$b_sample <- factor(rep(NA, n_sample))
    map$log_sd_sample <- factor(NA)
  }
  
  if (!random_sample_otu) {
    map$b_sample_otu <- factor(rep(NA, n_sample_otu))
    map$log_sd_sample_otu <- factor(NA)
  }
  
  if (!use_zi_otu) {
    map$b_zi_otu <- factor(rep(NA, n_otu))
    map$log_sd_zi_otu <- factor(NA)
  }
  
  # Distribution-specific parameters.
  if (abundance_family == "poisson") {
    map$log_theta <- factor(NA)
    map$zi_intercept <- factor(NA)
  }
  
  if (abundance_family == "nbinom") {
    map$zi_intercept <- factor(NA)
  }
  
  if (abundance_family == "zip") {
    map$log_theta <- factor(NA)
  }
  
  
  # ==========================================================
  # 22. Tell TMB which parameters are random
  # ==========================================================
  
  random_effects <- character(0)
  
  if (random_occ_otu) {
    random_effects <- c(
      random_effects,
      "b_occ_otu"
    )
  }
  
  if (random_capture_otu) {
    random_effects <- c(
      random_effects,
      "b_cap_otu"
    )
  }
  
  if (random_abund_otu) {
    random_effects <- c(
      random_effects,
      "b_abund_otu"
    )
  }
  
  if (random_sample) {
    random_effects <- c(
      random_effects,
      "b_sample"
    )
  }
  
  if (random_sample_otu) {
    random_effects <- c(
      random_effects,
      "b_sample_otu"
    )
  }
  
  if (use_zi_otu) {
    random_effects <- c(
      random_effects,
      "b_zi_otu"
    )
  }
  
  
  # ==========================================================
  # 23. Print basic model information
  # ==========================================================
  
  if (verbose) {
    
    message("Sites x OTUs: ", nrow(site_df))
    message("Samples x OTUs: ", nrow(sample_df))
    message("Read observations: ", nrow(dat))
    message("OTUs: ", n_otu)
    message("Biological samples: ", n_sample)
    
    message(
      "Random effects: ",
      if (length(random_effects) == 0L) {
        "none"
      } else {
        paste(random_effects, collapse = ", ")
      }
    )
    
    if (
      random_zi_otu &&
      !(abundance_family %in% c("zip", "zinb"))
    ) {
      message(
        "Note: random_zi_otu = TRUE is ignored because ",
        "abundance_family = '",
        abundance_family,
        "' has no zero-inflation component."
      )
    }
    
    message(
      "sdreport report covariance: ",
      if (get_report_covariance) {
        "enabled"
      } else {
        "disabled (memory-efficient default)"
      }
    )
  }
  
  # ==========================================================
  
  
  # ==========================================================
  # 24. BUILD TMB OBJECTIVE
  # ==========================================================
  
  if (verbose) {
    
    message("")
    message("Creating TMB objective...")
    
    message(
      "TMB internal trace: ",
      if (tmb_verbose) "enabled" else "disabled"
    )
  }
  
  obj <- TMB::MakeADFun(
    
    data = data_tmb,
    
    parameters = parameters,
    
    random =
      if (length(random_effects) == 0L) {
        NULL
      } else {
        random_effects
      },
    
    map =
      if (length(map) == 0L) {
        NULL
      } else {
        map
      },
    
    DLL = DLL,
    
    # IMPORTANT:
    #
    # verbose controls eDNAModel messages.
    # tmb_verbose controls TMB's iter/mgc trace.
    silent = !tmb_verbose
  )
  
  
  # ==========================================================
  # 25. INITIAL OBJECTIVE
  # ==========================================================
  
  initial_nll <- obj$fn(
    obj$par
  )
  
  if (!is.finite(initial_nll)) {
    
    stop(
      "Initial joint negative log-likelihood is not finite."
    )
  }
  
  if (verbose) {
    
    message(
      "Initial negative log-likelihood: ",
      format(
        initial_nll,
        digits = 12
      )
    )
    
    message(
      "Optimizing joint likelihood..."
    )
  }
  
  
  # ==========================================================
  # 26. GRADIENT EVALUATION
  # ==========================================================
  
  evaluate_gradient <- function(par) {
    
    # Update inner random-effect mode at supplied fixed
    # parameter vector.
    invisible(
      obj$fn(par)
    )
    
    g <- tryCatch(
      
      obj$gr(par),
      
      error = function(e) {
        
        rep(
          NA_real_,
          length(par)
        )
      }
    )
    
    if (
      length(g) == length(par)
    ) {
      
      names(g) <- names(par)
    }
    
    gradient_finite <-
      length(g) > 0L &&
      all(is.finite(g))
    
    max_gradient <-
      if (gradient_finite) {
        
        max(
          abs(g)
        )
        
      } else {
        
        Inf
      }
    
    list(
      
      gradient = g,
      
      finite =
        gradient_finite,
      
      max_abs_gradient =
        max_gradient,
      
      gradient_ok =
        gradient_finite &&
        max_gradient <= gradient_tol
    )
  }
  
  
  # ==========================================================
  # 27. HELPER: OBJECTIVE EQUIVALENCE
  # ==========================================================
  
  objective_tolerance <- function(x, y) {
    
    objective_rel_tol *
      max(
        1,
        abs(x),
        abs(y)
      )
  }
  
  
  objectives_equivalent <- function(x, y) {
    
    if (
      !is.finite(x) ||
      !is.finite(y)
    ) {
      return(FALSE)
    }
    
    abs(x - y) <=
      objective_tolerance(x, y)
  }
  
  
  # ==========================================================
  # 28. OPTIMIZATION + RESTARTS
  # ==========================================================
  
  optimization_history <- list()
  
  current_start <- obj$par
  
  max_passes <-
    max_restarts + 1L
  
  stalled <- FALSE
  
  
  for (
    pass in seq_len(max_passes)
  ) {
    
    if (verbose) {
      
      message("")
      message(
        "Optimization pass ",
        pass,
        " of ",
        max_passes,
        "..."
      )
    }
    
    
    opt_pass <- stats::nlminb(
      
      start =
        current_start,
      
      objective =
        obj$fn,
      
      gradient =
        obj$gr,
      
      control = list(
        
        iter.max =
          iter_max,
        
        eval.max =
          eval_max,
        
        rel.tol =
          rel_tol
      )
    )
    
    
    grad_pass <- evaluate_gradient(
      opt_pass$par
    )
    
    
    optimizer_ok_pass <-
      identical(
        as.integer(
          opt_pass$convergence
        ),
        0L
      )
    
    
    # ------------------------------------------------------
    # Store COMPLETE optimization pass
    # ------------------------------------------------------
    
    optimization_history[[pass]] <- list(
      
      pass = pass,
      
      opt = opt_pass,
      
      par = opt_pass$par,
      
      objective =
        opt_pass$objective,
      
      convergence_code =
        opt_pass$convergence,
      
      optimizer_message =
        opt_pass$message,
      
      optimizer_ok =
        optimizer_ok_pass,
      
      gradient =
        grad_pass$gradient,
      
      max_abs_gradient =
        grad_pass$max_abs_gradient,
      
      gradient_ok =
        grad_pass$gradient_ok
    )
    
    
    # ------------------------------------------------------
    # Print pass summary
    # ------------------------------------------------------
    
    if (verbose) {
      
      message(
        "  NLL: ",
        format(
          opt_pass$objective,
          digits = 13
        )
      )
      
      message(
        "  Optimizer code: ",
        opt_pass$convergence
      )
      
      message(
        "  Optimizer message: ",
        opt_pass$message
      )
      
      message(
        "  Maximum absolute gradient: ",
        signif(
          grad_pass$max_abs_gradient,
          9
        )
      )
      
      message(
        "  Strict gradient target met: ",
        grad_pass$gradient_ok
      )
    }
    
    
    # ------------------------------------------------------
    # Strict optimizer/gradient criteria reached
    # ------------------------------------------------------
    
    if (
      optimizer_ok_pass &&
      grad_pass$gradient_ok
    ) {
      
      if (verbose) {
        
        message(
          "  Strict optimizer/gradient criteria satisfied."
        )
      }
      
      break
    }
    
    
    # ------------------------------------------------------
    # Check restart behavior
    # ------------------------------------------------------
    
    if (pass > 1L) {
      
      previous <-
        optimization_history[[pass - 1L]]
      
      if (
        is.finite(
          previous$objective
        ) &&
        is.finite(
          opt_pass$objective
        ) &&
        is.finite(
          previous$max_abs_gradient
        ) &&
        is.finite(
          grad_pass$max_abs_gradient
        )
      ) {
        
        obj_change <-
          abs(
            opt_pass$objective -
              previous$objective
          )
        
        
        obj_tol <-
          objective_tolerance(
            opt_pass$objective,
            previous$objective
          )
        
        
        objective_stable <-
          obj_change <= obj_tol
        
        
        gradient_improvement <-
          previous$max_abs_gradient -
          grad_pass$max_abs_gradient
        
        
        gradient_meaningfully_improved <-
          gradient_improvement >
          gradient_improvement_tol
        
        
        if (verbose) {
          
          message(
            "  Objective change: ",
            signif(
              obj_change,
              8
            )
          )
          
          message(
            "  Objective equivalence tolerance: ",
            signif(
              obj_tol,
              8
            )
          )
          
          message(
            "  Gradient improvement: ",
            signif(
              gradient_improvement,
              8
            )
          )
        }
        
        
        # Stop if likelihood is effectively unchanged
        # and the gradient has not improved enough.
        stalled <-
          objective_stable &&
          !gradient_meaningfully_improved
        
        
        if (stalled) {
          
          if (verbose) {
            
            message(
              "  Restart sequence stalled."
            )
            
            message(
              "  Objective is effectively unchanged ",
              "and gradient improvement is below ",
              "gradient_improvement_tol."
            )
          }
          
          break
        }
      }
    }
    
    
    if (pass >= max_passes) {
      break
    }
    
    
    if (verbose) {
      
      message(
        "  Restarting optimization from current solution..."
      )
    }
    
    
    current_start <-
      opt_pass$par
  }
  
  
  # ==========================================================
  # 29. SELECT BEST OPTIMIZATION PASS
  #
  # IMPORTANT:
  #
  # Do NOT simply retain the final pass.
  #
  # First identify the minimum objective. Then consider all
  # solutions with effectively equivalent likelihoods and
  # select the one having the smallest maximum gradient.
  # ==========================================================
  
  objectives <- vapply(
    
    optimization_history,
    
    function(x) {
      x$objective
    },
    
    numeric(1)
  )
  
  
  gradients <- vapply(
    
    optimization_history,
    
    function(x) {
      x$max_abs_gradient
    },
    
    numeric(1)
  )
  
  
  finite_objective <-
    is.finite(objectives)
  
  
  if (!any(finite_objective)) {
    
    stop(
      "All optimization passes returned non-finite objectives."
    )
  }
  
  
  min_objective <-
    min(
      objectives[
        finite_objective
      ]
    )
  
  
  equivalent_to_best <- vapply(
    
    seq_along(objectives),
    
    function(i) {
      
      is.finite(
        objectives[i]
      ) &&
        objectives_equivalent(
          objectives[i],
          min_objective
        )
    },
    
    logical(1)
  )
  
  
  candidate_passes <- which(
    
    equivalent_to_best &
      is.finite(gradients)
  )
  
  
  if (
    length(candidate_passes) == 0L
  ) {
    
    # Defensive fallback:
    # use minimum finite objective.
    selected_pass <- which.min(
      ifelse(
        finite_objective,
        objectives,
        Inf
      )
    )
    
  } else {
    
    selected_pass <-
      candidate_passes[
        which.min(
          gradients[
            candidate_passes
          ]
        )
      ]
  }
  
  
  selected <-
    optimization_history[[selected_pass]]
  
  
  opt <- selected$opt
  
  
  # Re-evaluate TMB at selected solution.
  invisible(
    obj$fn(opt$par)
  )
  
  
  final_grad_info <-
    evaluate_gradient(
      opt$par
    )
  
  
  final_gradient <-
    final_grad_info$gradient
  
  gradient_finite <-
    final_grad_info$finite
  
  max_abs_gradient <-
    final_grad_info$max_abs_gradient
  
  gradient_ok <-
    final_grad_info$gradient_ok
  
  
  optimizer_ok <-
    identical(
      as.integer(
        opt$convergence
      ),
      0L
    )
  
  
  if (verbose) {
    
    message("")
    message(
      "Selected optimization pass ",
      selected_pass,
      " of ",
      length(
        optimization_history
      ),
      "."
    )
    
    message(
      "  Selected NLL: ",
      format(
        opt$objective,
        digits = 13
      )
    )
    
    message(
      "  Selected maximum absolute gradient: ",
      signif(
        max_abs_gradient,
        9
      )
    )
  }
  
  
  # ==========================================================
  # 30. GRADIENT STATUS
  # ==========================================================
  
  gradient_status <-
    
    if (!gradient_finite) {
      
      "FAIL"
      
    } else if (
      max_abs_gradient <=
      gradient_tol
    ) {
      
      "PASS"
      
    } else if (
      max_abs_gradient <=
      gradient_marginal_factor *
      gradient_tol
    ) {
      
      "MARGINAL"
      
    } else {
      
      "FAIL"
    }
  
  
  # ==========================================================
  # 31. GRADIENT TABLE
  # ==========================================================
  
  gradient_table <- data.frame(
    
    parameter =
      names(final_gradient),
    
    gradient =
      as.numeric(
        final_gradient
      ),
    
    abs_gradient =
      abs(
        as.numeric(
          final_gradient
        )
      ),
    
    stringsAsFactors =
      FALSE
  )
  
  
  gradient_table <-
    gradient_table[
      order(
        gradient_table$abs_gradient,
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]
  
  
  n_show <- min(
    n_gradient_report,
    nrow(
      gradient_table
    )
  )
  
  
  largest_gradients <-
    
    if (n_show > 0L) {
      
      gradient_table[
        seq_len(n_show),
        ,
        drop = FALSE
      ]
      
    } else {
      
      gradient_table
    }
  
  
  # ==========================================================
  # 32. SDREPORT
  # ==========================================================
  
  if (verbose) {
    
    message("")
    message(
      "Optimization completed."
    )
    
    message(
      "Starting TMB::sdreport()..."
    )
    
    message(
      "  getReportCovariance = ",
      get_report_covariance
    )
    
    message(
      "  getJointPrecision = ",
      get_joint_precision
    )
    
    message(
      "This step may take substantially longer than ",
      "optimization for large random-effects models."
    )
  }
  
  
  sdreport_error <- NULL
  
  
  sdreport_start_time <-
    proc.time()[["elapsed"]]
  
  
  sdr <- tryCatch(
    
    TMB::sdreport(
      
      obj,
      
      par.fixed =
        opt$par,
      
      getJointPrecision =
        isTRUE(
          get_joint_precision
        ) &&
        length(
          random_effects
        ) > 0L,
      
      getReportCovariance =
        isTRUE(
          get_report_covariance
        )
    ),
    
    error = function(e) {
      
      sdreport_error <<-
        conditionMessage(e)
      
      NULL
    }
  )
  
  
  sdreport_elapsed <-
    proc.time()[["elapsed"]] -
    sdreport_start_time
  
  
  sdreport_ok <-
    !is.null(sdr)
  
  
  if (verbose) {
    
    if (sdreport_ok) {
      
      message(
        "TMB::sdreport() completed successfully."
      )
      
      message(
        "  Elapsed time: ",
        round(
          sdreport_elapsed,
          2
        ),
        " seconds"
      )
      
    } else {
      
      message(
        "TMB::sdreport() FAILED."
      )
      
      if (
        !is.null(
          sdreport_error
        )
      ) {
        
        message(
          "  Error: ",
          sdreport_error
        )
      }
    }
  }
  
  
  # ==========================================================
  # 33. HESSIAN
  # ==========================================================
  
  pd_hessian <-
    
    if (
      sdreport_ok &&
      !is.null(
        sdr$pdHess
      )
    ) {
      
      isTRUE(
        sdr$pdHess
      )
      
    } else {
      
      FALSE
    }
  
  
  # ==========================================================
  # 34. TMB REPORT
  # ==========================================================
  
  report <- tryCatch(
    
    obj$report(),
    
    error = function(e) {
      
      warning(
        "TMB report failed: ",
        conditionMessage(e),
        call. = FALSE
      )
      
      NULL
    }
  )
  
  
  # ==========================================================
  # 35. FIXED PARAMETER SUMMARY
  # ==========================================================
  
  if (sdreport_ok) {
    
    fixed_matrix <- tryCatch(
      
      summary(
        sdr,
        "fixed"
      ),
      
      error =
        function(e) NULL
    )
    
    
    if (!is.null(fixed_matrix)) {
      
      fixed_summary <-
        as.data.frame(
          fixed_matrix
        )
      
      fixed_summary$parameter <-
        rownames(
          fixed_summary
        )
      
      rownames(
        fixed_summary
      ) <- NULL
      
    } else {
      
      fixed_summary <-
        data.frame()
    }
    
  } else {
    
    fixed_summary <-
      data.frame()
  }
  
  
  # ==========================================================
  # 36. ADREPORT SUMMARY
  # ==========================================================
  
  if (sdreport_ok) {
    
    derived_matrix <- tryCatch(
      
      summary(
        sdr,
        "report"
      ),
      
      error =
        function(e) NULL
    )
    
    
    if (!is.null(derived_matrix)) {
      
      derived_summary <-
        as.data.frame(
          derived_matrix
        )
      
      derived_summary$parameter <-
        rownames(
          derived_summary
        )
      
      rownames(
        derived_summary
      ) <- NULL
      
    } else {
      
      derived_summary <-
        data.frame()
    }
    
  } else {
    
    derived_summary <-
      data.frame()
  }
  
  
  # ==========================================================
  # 37. FINITE STANDARD ERRORS
  # ==========================================================
  
  finite_se <-
    
    if (
      nrow(
        fixed_summary
      ) > 0L &&
      "Std. Error" %in%
      names(
        fixed_summary
      )
    ) {
      
      all(
        is.finite(
          fixed_summary[["Std. Error"]]
        )
      )
      
    } else {
      
      FALSE
    }
  
  
  # ==========================================================
  # 38. HEURISTIC: NEAR-ZERO RANDOM-EFFECT SD
  # ==========================================================
  
  near_zero_sd <-
    data.frame()
  
  
  if (
    nrow(fixed_summary) > 0L &&
    all(
      c(
        "Estimate",
        "parameter"
      ) %in%
      names(
        fixed_summary
      )
    )
  ) {
    
    log_sd_rows <- grepl(
      "^log_sd",
      fixed_summary$parameter
    )
    
    
    if (any(log_sd_rows)) {
      
      tmp_sd <-
        fixed_summary[
          log_sd_rows,
          ,
          drop = FALSE
        ]
      
      
      tmp_sd$sd_estimate <-
        exp(
          tmp_sd$Estimate
        )
      
      
      tmp_sd$near_zero <-
        tmp_sd$Estimate <
        log_sd_warning_threshold
      
      
      near_zero_sd <-
        tmp_sd[
          tmp_sd$near_zero,
          ,
          drop = FALSE
        ]
    }
  }
  
  
  variance_status <-
    
    if (
      nrow(
        near_zero_sd
      ) == 0L
    ) {
      
      "PASS"
      
    } else {
      
      "WARNING"
    }
  
  
  # ==========================================================
  # 39. HEURISTIC: LARGE SE / ESTIMATE
  # ==========================================================
  
  large_se_parameters <-
    data.frame()
  
  
  if (
    nrow(fixed_summary) > 0L &&
    all(
      c(
        "Estimate",
        "Std. Error",
        "parameter"
      ) %in%
      names(
        fixed_summary
      )
    )
  ) {
    
    tmp <- fixed_summary
    
    
    tmp$se_to_estimate_ratio <-
      NA_real_
    
    
    valid_ratio <-
      
      is.finite(
        tmp$Estimate
      ) &
      
      is.finite(
        tmp[["Std. Error"]]
      ) &
      
      abs(
        tmp$Estimate
      ) >
      estimate_zero_tol
    
    
    tmp$se_to_estimate_ratio[
      valid_ratio
    ] <-
      
      tmp[["Std. Error"]][valid_ratio] /
      
      abs(
        tmp$Estimate[
          valid_ratio
        ]
      )
    
    
    tmp$large_relative_se <-
      
      !is.na(
        tmp$se_to_estimate_ratio
      ) &
      
      tmp$se_to_estimate_ratio >
      se_estimate_ratio_threshold
    
    
    large_se_parameters <-
      
      tmp[
        tmp$large_relative_se,
        ,
        drop = FALSE
      ]
  }
  
  
  uncertainty_status <-
    
    if (
      nrow(
        large_se_parameters
      ) == 0L
    ) {
      
      "PASS"
      
    } else {
      
      "WARNING"
    }
  
  
  # ==========================================================
  # 40. INDIVIDUAL STATUS
  # ==========================================================
  
  optimizer_status <-
    if (optimizer_ok) {
      "PASS"
    } else {
      "FAIL"
    }
  
  
  sdreport_status <-
    if (sdreport_ok) {
      "PASS"
    } else {
      "FAIL"
    }
  
  
  hessian_status <-
    if (pd_hessian) {
      "PASS"
    } else {
      "FAIL"
    }
  
  
  se_status <-
    if (finite_se) {
      "PASS"
    } else {
      "FAIL"
    }
  
  
  # ==========================================================
  # 41. STRICT + ACCEPTABLE CONVERGENCE
  # ==========================================================
  
  strict_convergence <-
    
    optimizer_ok &&
    gradient_status == "PASS" &&
    sdreport_ok &&
    pd_hessian &&
    finite_se
  
  
  acceptable_convergence <-
    
    optimizer_ok &&
    gradient_status %in%
    c(
      "PASS",
      "MARGINAL"
    ) &&
    sdreport_ok &&
    pd_hessian &&
    finite_se
  
  
  overall_status <-
    
    if (strict_convergence) {
      
      "PASS"
      
    } else if (
      acceptable_convergence
    ) {
      
      "MARGINAL"
      
    } else {
      
      "FAIL"
    }
  
  
  # Keep old field for compatibility.
  converged <-
    acceptable_convergence
  
  
  # ==========================================================
  # 42. FAILURE REASONS
  # ==========================================================
  
  failure_reasons <-
    character(0)
  
  
  if (!optimizer_ok) {
    
    failure_reasons <- c(
      failure_reasons,
      paste0(
        "Optimizer returned convergence code ",
        opt$convergence,
        "."
      )
    )
  }
  
  
  if (
    gradient_status == "MARGINAL"
  ) {
    
    failure_reasons <- c(
      failure_reasons,
      paste0(
        "Maximum absolute gradient exceeds the strict ",
        "tolerance but remains within the marginal range."
      )
    )
  }
  
  
  if (
    gradient_status == "FAIL"
  ) {
    
    failure_reasons <- c(
      failure_reasons,
      paste0(
        "Maximum absolute gradient exceeds the ",
        "marginal gradient threshold."
      )
    )
  }
  
  
  if (!sdreport_ok) {
    
    failure_reasons <- c(
      failure_reasons,
      "TMB::sdreport() was not successful."
    )
  }
  
  
  if (
    sdreport_ok &&
    !pd_hessian
  ) {
    
    failure_reasons <- c(
      failure_reasons,
      "The Hessian is not positive definite."
    )
  }
  
  
  if (
    sdreport_ok &&
    !finite_se
  ) {
    
    failure_reasons <- c(
      failure_reasons,
      "One or more fixed-parameter standard errors are non-finite."
    )
  }
  
  
  if (
    length(
      failure_reasons
    ) == 0L
  ) {
    
    failure_reasons <-
      "All numerical convergence criteria were satisfied."
  }
  
  
  # ==========================================================
  # 43. DIAGNOSTIC TABLE
  # ==========================================================
  
  diagnostic_table <- data.frame(
    
    Diagnostic = c(
      
      "Optimizer convergence",
      
      "Maximum absolute gradient",
      
      "sdreport",
      
      "Positive-definite Hessian",
      
      "Finite fixed-parameter SEs",
      
      "Near-zero variance components",
      
      "Large SE relative to estimate",
      
      "Overall numerical convergence"
    ),
    
    Value = c(
      
      paste0(
        "code ",
        opt$convergence
      ),
      
      format(
        max_abs_gradient,
        digits = 8
      ),
      
      if (
        sdreport_ok
      ) {
        "successful"
      } else {
        "failed"
      },
      
      as.character(
        pd_hessian
      ),
      
      as.character(
        finite_se
      ),
      
      as.character(
        nrow(
          near_zero_sd
        )
      ),
      
      as.character(
        nrow(
          large_se_parameters
        )
      ),
      
      as.character(
        acceptable_convergence
      )
    ),
    
    Status = c(
      
      optimizer_status,
      
      gradient_status,
      
      sdreport_status,
      
      hessian_status,
      
      se_status,
      
      variance_status,
      
      uncertainty_status,
      
      overall_status
    ),
    
    stringsAsFactors =
      FALSE
  )
  
  
  # ==========================================================
  # 44. FINAL VERBOSE REPORT
  # ==========================================================
  
  if (verbose) {
    
    cat("\n")
    cat(
      "----- Numerical convergence diagnostics -----\n"
    )
    
    cat(
      "Selected optimization pass: ",
      selected_pass,
      " / ",
      length(
        optimization_history
      ),
      "\n",
      sep = ""
    )
    
    cat(
      "Final negative log-likelihood: ",
      format(
        opt$objective,
        digits = 13
      ),
      "\n\n",
      sep = ""
    )
    
    cat(
      "Optimizer convergence code: ",
      opt$convergence,
      "\n",
      sep = ""
    )
    
    cat(
      "Optimizer status:           ",
      optimizer_status,
      "\n\n",
      sep = ""
    )
    
    cat(
      "Maximum absolute gradient:  ",
      format(
        max_abs_gradient,
        digits = 8
      ),
      "\n",
      sep = ""
    )
    
    cat(
      "Strict gradient tolerance:  ",
      gradient_tol,
      "\n",
      sep = ""
    )
    
    cat(
      "Marginal gradient limit:    ",
      gradient_marginal_factor *
        gradient_tol,
      "\n",
      sep = ""
    )
    
    cat(
      "Gradient status:            ",
      gradient_status,
      "\n\n",
      sep = ""
    )
    
    cat(
      "sdreport status:            ",
      sdreport_status,
      "\n",
      sep = ""
    )
    
    cat(
      "Positive-definite Hessian:  ",
      pd_hessian,
      " [",
      hessian_status,
      "]\n",
      sep = ""
    )
    
    cat(
      "Finite fixed-parameter SEs: ",
      finite_se,
      " [",
      se_status,
      "]\n\n",
      sep = ""
    )
    
    cat(
      "Strict convergence:         ",
      strict_convergence,
      "\n",
      sep = ""
    )
    
    cat(
      "Acceptable convergence:     ",
      acceptable_convergence,
      "\n",
      sep = ""
    )
    
    cat(
      "Overall status:             ",
      overall_status,
      "\n",
      sep = ""
    )
    
    
    if (
      gradient_status != "PASS" &&
      nrow(
        largest_gradients
      ) > 0L
    ) {
      
      cat("\n")
      cat(
        "Largest gradient parameters:\n"
      )
      
      for (
        i in seq_len(
          nrow(
            largest_gradients
          )
        )
      ) {
        
        cat(
          sprintf(
            "  %-20s = %.8g\n",
            
            largest_gradients$parameter[i],
            
            largest_gradients$abs_gradient[i]
          )
        )
      }
    }
    
    
    cat(
      "---------------------------------------------\n"
    )
  }
  
  
  # ==========================================================
  # 45. WARNING
  # ==========================================================
  
  if (
    overall_status == "MARGINAL"
  ) {
    
    warning(
      
      sprintf(
        paste0(
          "Fit is numerically MARGINAL: optimizer, ",
          "Hessian, sdreport and SE diagnostics passed, ",
          "but maximum absolute gradient %.6g exceeds ",
          "the strict tolerance %.6g."
        ),
        
        max_abs_gradient,
        gradient_tol
      ),
      
      call. = FALSE
    )
    
  } else if (
    overall_status == "FAIL"
  ) {
    
    warning(
      
      paste0(
        "Fit did not satisfy numerical convergence criteria. ",
        paste(
          failure_reasons,
          collapse = " "
        )
      ),
      
      call. = FALSE
    )
  }
  
  
  # ==========================================================
  # 46. RETURN OBJECT
  # ==========================================================
  
  out <- list(
    
    fit = opt,
    
    tmb_object = obj,
    
    sdreport = sdr,
    
    report = report,
    
    fixed_effects =
      fixed_summary,
    
    derived =
      derived_summary,
    
    site_data =
      site_df,
    
    sample_data =
      sample_df,
    
    long_df =
      dat,
    
    otu_stats =
      otu_stats,
    
    retained_otus =
      retained_otus,
    
    data_tmb =
      data_tmb,
    
    parameters_start =
      parameters,
    
    parameter_map =
      map,
    
    random_effects =
      random_effects,
    
    formulas = list(
      
      occupancy =
        occupancy_formula,
      
      capture =
        capture_formula,
      
      abundance =
        abundance_formula
    ),
    
    abundance_family =
      abundance_family,
    
    
    convergence = list(
      
      # Backward-compatible overall result
      converged =
        converged,
      
      strict_convergence =
        strict_convergence,
      
      acceptable_convergence =
        acceptable_convergence,
      
      overall_status =
        overall_status,
      
      optimizer_ok =
        optimizer_ok,
      
      optimizer_status =
        optimizer_status,
      
      optimizer_code =
        opt$convergence,
      
      optimizer_message =
        opt$message,
      
      code =
        opt$convergence,
      
      message =
        opt$message,
      
      objective =
        opt$objective,
      
      gradient =
        final_gradient,
      
      max_abs_gradient =
        max_abs_gradient,
      
      gradient_tolerance =
        gradient_tol,
      
      gradient_marginal_factor =
        gradient_marginal_factor,
      
      gradient_ok =
        gradient_ok,
      
      gradient_status =
        gradient_status,
      
      gradient_finite =
        gradient_finite,
      
      gradient_table =
        gradient_table,
      
      largest_gradients =
        largest_gradients,
      
      sdreport_ok =
        sdreport_ok,
      
      sdreport_status =
        sdreport_status,
      
      sdreport_error =
        sdreport_error,
      
      sdreport_elapsed_seconds =
        sdreport_elapsed,
      
      pd_hessian =
        pd_hessian,
      
      hessian_status =
        hessian_status,
      
      finite_standard_errors =
        finite_se,
      
      finite_fixed_parameter_standard_errors =
        finite_se,
      
      se_status =
        se_status,
      
      failure_reasons =
        failure_reasons,
      
      optimization_history =
        optimization_history,
      
      n_optimization_passes =
        length(
          optimization_history
        ),
      
      selected_pass =
        selected_pass,
      
      max_restarts =
        max_restarts,
      
      objective_rel_tol =
        objective_rel_tol,
      
      gradient_improvement_tol =
        gradient_improvement_tol,
      
      stalled =
        stalled
    ),
    
    
    diagnostics = list(
      
      table =
        diagnostic_table,
      
      near_zero_sd =
        near_zero_sd,
      
      near_zero_sd_status =
        variance_status,
      
      log_sd_warning_threshold =
        log_sd_warning_threshold,
      
      large_se_parameters =
        large_se_parameters,
      
      large_se_status =
        uncertainty_status,
      
      se_estimate_ratio_threshold =
        se_estimate_ratio_threshold
    ),
    
    
    sdreport_settings = list(
      
      getReportCovariance =
        get_report_covariance,
      
      getJointPrecision =
        get_joint_precision
    ),
    
    
    zero_inflation_structure = list(
      
      enabled =
        abundance_family %in%
        c(
          "zip",
          "zinb"
        ),
      
      random_zi_otu =
        use_zi_otu
    ),
    
    
    note = paste(
      
      "Joint observed-data likelihood.",
      
      "Latent occupancy Z and biological-sample capture A",
      
      "are analytically marginalized.",
      
      "Enabled Gaussian random effects are integrated",
      
      "using the TMB Laplace approximation.",
      
      "Strict numerical convergence requires optimizer code 0,",
      
      "a finite outer gradient no greater than gradient_tol,",
      
      "successful sdreport, a positive-definite Hessian,",
      
      "and finite fixed-parameter standard errors.",
      
      "A MARGINAL gradient status is reported separately",
      
      "when the gradient exceeds the strict tolerance but",
      
      "remains below gradient_marginal_factor * gradient_tol.",
      
      "Optimization restarts retain the smallest-gradient",
      
      "solution among effectively equivalent likelihoods.",
      
      "getReportCovariance and getJointPrecision are FALSE",
      
      "by default to reduce memory use."
    )
  )
  
  
  class(out) <- c(
    "eDNAModel_joint",
    "list"
  )
  
  
  return(out)
}
        
