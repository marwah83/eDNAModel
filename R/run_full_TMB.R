#' @title Run Full TMB Model
#' @description Fits a hierarchical multi-species occupancy-abundance model using TMB...
#'
#' @param data_array_filtered A 3D array (Species x Sites x Replicates) after filtering.
#' @param X A data frame containing covariates for model matrices.
#' @param a.formula Formula for the abundance model.
#' @param o.formula Formula for the occupancy model.
#' @param linko Link function for occupancy: 1 = logit, 2 = probit.
#' @param linka Link function for abundance: 0 = log, 1 = logit, 2 = probit, 3 = cloglog.
#' @param family Distribution: 0 = ZIP, 1 = ZINB, 2 = Binomial.
#' @param Ntrials Matrix of trials for Binomial models.
#' @param offset Offset matrix (same dimension as y)
#' @param control Optimization control list.
#' @param data_array_filtered Filtered 3D data array ready for TMB model fitting.
#'
#' @return A list with optimization results, fitted TMB object, and estimated probabilities.
#'
#' @export
#'
#' @useDynLib eDNAModel, .registration = TRUE
#' @importFrom TMB MakeADFun
#' @importFrom TMB openmp
#' @importFrom methods as
#' @importFrom stats aggregate model.matrix optim plogis
#' @import Matrix
#' @seealso \code{\link{summary.eDNAModel}}

run_full_TMB <- function(y,
                         X,
                         a.formula = ~1,
                         o.formula = ~1,
                         linko = 1,
                         linka = 0,
                         family = 1,
                         Ntrials = matrix(0), offset = NULL,
                         control = list(trace = TRUE,
                         startOptcontrol = list(maxit = 200),
                         optControl = list(maxit = 10e3))) {
method = "LBFGS"
 stopifnot(is.matrix(y))
  stopifnot(all(rownames(y) %in% rownames(X)))
  X <- X[rownames(y), , drop = FALSE]

  # Ensure required columns
  required_cols <- c("Site", "Sample", "Replicate")
  stopifnot(all(required_cols %in% colnames(X)))

  # Ensure factor types
  X$Site <- as.factor(X$Site)
  X$Sample <- as.factor(X$Sample)
  X$Replicate <- as.factor(X$Replicate)

  # Site index for TMB (0-based)
  sites <- as.numeric(X$Site) - 1

  # Aggregate y over sites for occupancy
  ysites <- aggregate(y, by = list(site = sites), FUN = sum)
  ysites <- as.matrix(ysites[, -1])

  if (!linko %in% c(1, 2)) stop("linko must be 1 (logit) or 2 (probit)")
  if (!linka %in% c(0, 1, 2, 3)) stop("linka must be 0 (log), 1 (logit), 2 (probit), or 3 (cloglog)")
  if (!family %in% c(0, 1, 2)) stop("family must be 0 (ZIP), 1 (ZINB), or 2 (Binomial)")
  if (!is.matrix(Ntrials)) stop("Ntrials must be a matrix")

  # Abundance model matrices
  Xa <- model.matrix(gllvm:::nobars1_(a.formula), X)
  if (gllvm:::anyBars(a.formula)) {
    Zalist <- gllvm:::mkReTrms1(gllvm:::findbars1(a.formula), X)
    Za <- Matrix::t(Zalist$Zt)
    csa <- Zalist$cs
  } else {
    Za <- as(matrix(0), "TsparseMatrix")
    csa <- matrix(0)
  }

  # Occupancy model matrices
  site_levels <- sort(unique(sites))
  row_matches <- match(site_levels, sites)
  Xocc <- X[row_matches, , drop = FALSE]
  if (nrow(Xocc) != length(site_levels)) stop("Mismatch between site-level covariates and unique site IDs")

  Xo <- model.matrix(gllvm:::nobars1_(o.formula), Xocc)
  if (gllvm:::anyBars(o.formula)) {
    Zolist <- gllvm:::mkReTrms1(gllvm:::findbars1(o.formula), Xocc)
    Zo <- Matrix::t(Zolist$Zt)
    cso <- Zolist$cs
  } else {
    Zo <- as(matrix(0), "TsparseMatrix")
    cso <- matrix(0)
  }

  # Starting parameters
  Ba <- matrix(0, nrow = ncol(Xa), ncol = ncol(y))
  Bo <- matrix(0, nrow = ncol(Xo), ncol = ncol(y))
  Ua <- matrix(0)
  Uo <- matrix(0)
  logphi <- rep(0, ncol(y))
  logsda <- rep(0, nrow(Ua))
  corsa <- rep(1e-5, nrow(csa))
  logsdo <- rep(0, nrow(Uo))
  corso <- rep(1e-5, nrow(cso))

  if(is.null(offset)){
    offset <- matrix(0, ncol = ncol(y), nrow = nrow(y))
  }

  maplist <- list()

  # Initial TMB object
  fit <- TMB::MakeADFun(
    data = list(
      Y = y, Ysites = ysites, Xa = Xa, Xo = Xo,
      Za = as(matrix(0), "TsparseMatrix"), Zo = as(matrix(0),"TsparseMatrix"), family = family, sites = sites,
      csa = csa, cso = cso, NTrials = Ntrials,
      linka = linka, linko = linko, offset = offset
    ),
    parameters = list(
      Ba = Ba, Bo = Bo, Ua = Ua, Uo = Uo,
      logphi = logphi, logsda = logsda, corsa = corsa,
      logsdo = logsdo, corso = corso
    ),
    DLL = "eDNAModel",
    map = maplist
  )

  # First optimization
  opt <- minic::rnewton(fit$par, fit$fn, fit$gr, method = method,
               control = control$startOptControl, verbose = control$trace)

  # Reconstruct estimates for second run
  Ba2 <- matrix(opt$par[names(opt$par) == "Ba"], nrow = ncol(Xa))
  Bo2 <- matrix(opt$par[names(opt$par) == "Bo"], nrow = ncol(Xo))

  Ua <- if (nrow(Za) > 1) matrix(0, nrow = ncol(Za), ncol = ncol(y)) else matrix(0)
  Uo <- if (nrow(Zo) > 1) matrix(0, nrow = ncol(Zo), ncol = ncol(y)) else matrix(0)

  logsda <- rep(0, nrow(Ua))
  corsa <- rep(1e-5, nrow(csa))
  logsdo <- rep(0, nrow(Uo))
  corso <- rep(1e-5, nrow(cso))
  logphi <- opt$par[names(opt$par) == "logphi"]

  random <- NULL
  if (nrow(Za) > 1) random <- c(random, "Ua")
  if (nrow(Zo) > 1) random <- c(random, "Uo")

  # Final TMB fit
  fit <- TMB::MakeADFun(
    data = list(Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = Za, Zo = Zo,
                family = family, sites = sites, csa = csa, cso = cso,
                NTrials = Ntrials, linka = linka, linko = linko, offset = offset),
    parameters = list(Ba = Ba2, Bo = Bo2, Ua = Ua, Uo = Uo, logphi = logphi,
                      logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso),
    random = random,
    DLL = "eDNAModel", silent = FALSE,
    inner.control = list(mgcmax = 1e+200, tol10 = 0.01),
    map = maplist
  )

  opt2 <- minic::rnewton(fit$par, fit$fn, fit$gr, method = method,
                                control = control$optControl, verbose = control$trace)


  # Occupancy & Detection probability calculation
  etao <- fit$report(fit$env$last.par.best)$etao
  occup.prob <- 1 - plogis(etao)

  occup.prob_site <- occup.prob  # Already site-level?

  lambda_matrix <- matrix(exp(fit$report(fit$env$last.par.best)$etaa), ncol = ncol(occup.prob))
  site_ids <- as.numeric(sites)
  lambda_prod <- aggregate(exp(-lambda_matrix), by = list(site = site_ids), FUN = prod)[, -1]

  prob.detect <- 1 - occup.prob_site * lambda_prod

  out <- list(
    start = opt,
    optimization = opt2,
    occupancy_probability = as.matrix(occup.prob),
    detection_probability = as.matrix(prob.detect),
    TMBobj = fit
  )
  class(out) <- "eDNAModel"
  return(out)
}

