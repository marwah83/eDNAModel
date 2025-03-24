#' Full TMB Pipeline with Results Extraction (Complete & Preserved Code)
#'
#' This function runs a full TMB model fitting pipeline and returns occupancy, lambda, and detection probabilities.
#'
#' @useDynLib eDNAModel, .registration = TRUE
#' @importFrom TMB MakeADFun
#' @importFrom TMB openmp
#' @import Matrix
#'
#' @param data_array_filtered Filtered 3D data array ready for TMB model fitting.
#' @return List containing optimized object opt2 and final extracted results (occupancy, lambda, prob.detect)
#' @export
#'
library(TMB)

run_full_TMB <- function(data_array_filtered) {

  ## TMB stuff ##

  #compile("~/Desktop/Diversea/occ.cpp", framework = "TMBad", flags = "-O0 -g") # Only use flags for debugging
  #dyn.load("~/Desktop/Diversea/occ.so")
  ## end TMB stuff ##

  ## Data processing ##
  y <- data_array_filtered

  to2D <- function(array_3d) {
    species <- dim(array_3d)[1]
    sites <- dim(array_3d)[2]
    replicates <- dim(array_3d)[3]
    data_values <- matrix(array_3d, nrow = sites * replicates, ncol = species, byrow = FALSE)
    site_column <- rep(dimnames(array_3d)$Sites, each = replicates)
    replicate_column <- rep(dimnames(array_3d)$Replicates, times = sites)
    final_data <- data.frame(Site = site_column, Replicate = replicate_column, data_values)
    colnames(final_data)[3:(2 + species)] <- dimnames(array_3d)$Species
    return(final_data)
  }

  Y <- to2D(y)
  ## End Data ##

  ## Modeling ##
  y <- as.matrix(Y[,-c(1:2)])
  sites <- as.numeric(Y[,1]) - 1
  ysites <- as.matrix(aggregate(y, FUN = sum, list(sites)))[,-1]
  Xa <- model.matrix(~Site, Y)
  Zalist <- gllvm:::mkReTrms1(gllvm:::findbars1(~diag(1|Replicate)), Y)
  csa <- Zalist$cs
  Za <- Matrix::t(Zalist$Zt)
  Xo <- model.matrix(~0 + Site, data.frame(Site = as.factor(1:max(as.numeric(Y$Site)))))
  Zo <- as(matrix(0), "TsparseMatrix")
  cso <- matrix(0)

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
  #TMB::openmp(parallel::detectCores() - 1, autopar = TRUE, DLL = "eDNAModel")

  fit <- MakeADFun(
    data = list(Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = as(matrix(0), "TsparseMatrix"), Zo = as(matrix(0), "TsparseMatrix"),
                family = family, sites = sites, csa = matrix(0), cso = matrix(0), NTrials = NTrials, linka = linka, linko = linko),
    parameters = list(Ba = Ba, Bo = Bo, Ua = Ua, Uo = Uo, logphi = logphi, logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso),
    DLL = "eDNAModel", map = maplist
  )
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
