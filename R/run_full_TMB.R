#' Run Full TMB Model
#'
#' Fits a hierarchical multi-species occupancy-abundance model using TMB,
#' with support for covariate formulas and random effects.
#'
#' @param data_array_filtered A 3D array (Species x Sites x Replicates) after filtering.
#' @param covariate_data A data frame containing covariates for model matrices.
#' @param a.formula A formula for the abundance (observation) model (default `~ site + (1 | replicate)`).
#' @param o.formula A formula for the occurrence (occupancy) model (default `~ site`).
#'
#' @return A list containing the fitted TMB object and optimization result.
#' @export
#'@name run_full_TMB
run_full_TMB <- function(data_array_filtered,
                         covariate_data,
                         a.formula = ~ 1,
                         o.formula = ~ 1) {

  ## -------------------------------
  ## Helper: Convert 3D array to long
  ## -------------------------------
  to2D <- function(array_3d) {
    species <- dim(array_3d)[1]
    sites <- dim(array_3d)[2]
    replicates <- dim(array_3d)[3]
    data_values <- matrix(array_3d, nrow = sites * replicates, ncol = species, byrow = FALSE)
    site_column <- rep(dimnames(array_3d)$Sites, each = replicates)
    replicate_column <- rep(dimnames(array_3d)$Replicates, times = sites)
    final_data <- data.frame(site = factor(site_column), replicate = factor(replicate_column), data_values)
    colnames(final_data)[3:(2 + species)] <- dimnames(array_3d)$Species
    return(final_data)
  }

  ## -------------------------------
  ## Prepare data
  ## -------------------------------
  Y_long <- to2D(data_array_filtered)
  y <- as.matrix(Y_long[, -(1:2)])  # Remove site & replicate columns
  x <- covariate_data

  sites <- as.numeric(x$site) - 1
  ysites <- as.matrix(stats::aggregate(y, FUN = sum, by = list(sites)))[, -1]

  ## -------------------------------
  ## Model matrices
  ## -------------------------------
  Xa <- model.matrix(a.formula, x)
  Xo <- model.matrix(o.formula, x)

  ## Random effects for abundance
  Zalist <- gllvm:::mkReTrms1(gllvm:::findbars1(a.formula), x)
  Za <- Matrix::t(Zalist$Zt)
  csa <- Zalist$cs

  ## No random effects for occupancy
  Zo <- as(matrix(0), "TsparseMatrix")
  cso <- matrix(0)

  ## -------------------------------
  ## TMB model setup
  ## -------------------------------
  family <- 1
  linka <- 0
  linko <- 1
  NTrials <- matrix(0, nrow(y), ncol(y))

  Ba <- matrix(0, nrow = ncol(Xa), ncol = ncol(y))
  Bo <- matrix(0, nrow = ncol(Xo), ncol = ncol(y))
  Ua <- matrix(0)
  Uo <- matrix(0)
  logphi <- rep(0, ncol(y))
  logsda <- rep(0, nrow(Ua))
  corsa <- rep(1e-5, nrow(csa))
  logsdo <- rep(0, nrow(Uo))
  corso <- 0

  maplist <- list()

  fit <- TMB::MakeADFun(
    data = list(
      Y = y, Ysites = ysites, Xa = Xa, Xo = Xo,
      Za = Za, Zo = Zo,
      family = family, sites = sites,
      csa = csa, cso = cso,
      NTrials = NTrials,
      linka = linka, linko = linko
    ),
    parameters = list(
      Ba = Ba, Bo = Bo,
      Ua = Ua, Uo = Uo,
      logphi = logphi,
      logsda = logsda,
      corsa = corsa,
      logsdo = logsdo,
      corso = corso
    ),
    DLL = "eDNAModel",
    map = maplist
  )

  ## Optimization
  opt <- optim(fit$par, fit$fn, fit$gr, method = "L-BFGS-B", control = list(trace = 1, maxit = 5000))

  ## ADD THIS PART EXACTLY AS YOU PROVIDED ##
  Bas <- opt$par[names(opt$par)=="Ba"]
  if("Ba"%in%names(maplist)){
    Bas2 <- Bas[fit$env$map$Ba]
  }else{
    Bas2 <- Bas
  }
  Ba2 = matrix(Bas2,nrow=ncol(Xa))

  Bos <- opt$par[names(opt$par)=="Bo"]
  if("Bo"%in%names(maplist)){
    Bos2 <- Bos[fit$env$map$Bo]
  }else{
    Bos2 <- Bos
  }
  Bo2 = matrix(Bos2,nrow=ncol(Xo))

  if(nrow(Za)>1){
    Ua = matrix(0,nrow=ncol(Za),ncol=ncol(y))
  }else{
    Ua = matrix(0)
  }

  if(nrow(Zo)>1){
    Uo = matrix(0,nrow=ncol(Zo),ncol=ncol(y))
  }else{
    Uo = matrix(0)
  }

  logsda <- rep(0,nrow(Ua))
  corsa <- rep(1e-5,nrow(csa))
  logsdo <- rep(0, nrow(Uo))
  corso <- rep(1e-5,nrow(cso))
  if(family==1)logphi=opt$par[names(opt$par)=="logphi"]

  random = NULL
  if(nrow(Za)>1)random <- c(random, "Ua")
  if(nrow(Zo)>1)random <- c(random, "Uo")

  fit <- MakeADFun(
    data = list(Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = Za, Zo = Zo, family = family, sites = sites,
                csa = csa, cso = cso, NTrials = NTrials, linka = linka, linko = linko),
    parameters = list(Ba = Ba2, Bo = Bo2, Ua = Ua, Uo = Uo, logphi = logphi, logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso),
    random = random,
    DLL = "eDNAModel", silent = FALSE,
    inner.control = list(mgcmax = 1e+200, tol10 = 0.01),
    map = maplist
  )

  opt2 <- optim(fit$par, fit$fn, fit$gr, method = "L-BFGS-B", control = list(trace = 1, maxit = 1e6)) # optimize

  ## End modeling ##

  ## Some results ##
  occup.prob <- (1 - plogis(fit$report(fit$env$last.par.best)$etao))
  lambda <- exp(fit$report(fit$env$last.par.best)$etaa)
  prob.detect <- 1 - occup.prob * as.matrix(aggregate(exp(-lambda), FUN = prod, list(sites + 1)))[,-1]

  # Checks
  any_violations <- any(ysites[occup.prob[,] > 0.8] == 0)
  all_correct <- all(ysites[occup.prob[,] < 0.8] == 0)

  message(" Final TMB model fitting completed.")

  return(list(
    optimization = opt2,
    occupancy_probability = occup.prob,
    lambda = lambda,
    detection_probability = prob.detect,
    check_large_occupancy = any_violations,
    check_small_occupancy = all_correct
  ))
}
