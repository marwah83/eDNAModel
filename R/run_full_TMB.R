#' Full TMB Pipeline with Results Extraction (Complete & Preserved Code)
#' Fits a hierarchical multi-species occupancy-abundance model using TMB,
#' with support for covariate formulas and random effects.
#' @useDynLib eDNAModel, .registration = TRUE
#' @importFrom TMB MakeADFun
#' @importFrom TMB openmp
#' @import Matrix
#' #' @param data_array_filtered Filtered 3D data array ready for TMB model fitting.
#' @return List containing optimized object opt2 and final extracted results (occupancy, lambda, prob.detect)
#' 
#' @export
#'@name run_full_TMB
run_full_TMB <- function(data_array_filtered) {
  
  ## TMB stuff ##
  
  #compile("~/Desktop/Diversea/occ.cpp", framework = "TMBad", flags = "-O0 -g") # Only use flags for debugging
  #dyn.load("~/Desktop/Diversea/occ.so")
  ## end TMB stuff ##
  
  ## Data processing ##
  y <- data_array_filtered
  ## -------------------------------
  ## Helper: Convert 3D array to long
  ## -------------------------------
  # function to go from 3d array to matrix of replicates and saites by species
  to2D <- function(array_3d) {
    # Get dimensions
    species <- dim(array_3d)[1]
    sites <- dim(array_3d)[2]
    replicates <- dim(array_3d)[3]
    
    # Reshape the data into a 2D matrix
    data_values <- matrix(array_3d, nrow = sites * replicates, ncol = species, byrow = FALSE)
    
    # Create Sites and Replicates columns
    site_column <- rep(dimnames(array_3d)$Sites, each = replicates)
    replicate_column <- rep(dimnames(array_3d)$Replicates, times = sites)
    
    # Combine into a data frame
    final_data <- data.frame(
      Site = site_column,
      Replicate = replicate_column,
      data_values
    )
    
    # Set column names for species
    colnames(final_data)[3:(2 + species)] <- dimnames(array_3d)$Species
    
    return(final_data)
  }
  
  Y <- to2D(y)
  ## end Data stuff ##
  
  ## Start modeling! ##
  # readying list format for TMB
  # First "data" components
  y <- as.matrix(Y[,-c(1:2)])
  sites = as.numeric(Y[,1])-1
  ysites <- as.matrix(aggregate(y,FUN=sum,list(sites)))[,-1]#sum over replicates
  Xa <- model.matrix(~Site,Y)
  Zalist <- gllvm:::mkReTrms1(gllvm:::findbars1(~diag(1|Replicate)), Y)
  csa <- Zalist$cs # set to single column matrix if no REs
  Za <- Matrix::t(Zalist$Zt)
  Xo <- model.matrix(~0+Site,data.frame(Site=as.factor(1:max(as.numeric(Y$Site)))))
  Zo <- as(matrix(0), "TsparseMatrix")
  cso <- matrix(0)# set to single column matrix if no REs
  
  family = 1#0 is ZIP, 1 is ZINB, 2 is BINOM with Ntrials
  linka = 0 # abundance link function; 0 = log, 1 = logit,  2 = probit, 3 = cloglog
  linko = 1 # occupancy link function; 0 = log, 1 = logit,  2 = probit, 3 = cloglog
  NTrials = matrix(0)  # nrow(y) by ncol(y) matrix of size for binomial. should separately be specified with family = 2
  
  # Readying "parameters" component TMB list function
  Ba = matrix(0,nrow=ncol(Xa),ncol=ncol(y))
  Bo = matrix(0,nrow=ncol(Xo),ncol=ncol(y))
  Ua = matrix(0)
  Uo = matrix(0)
  logphi=rep(0,ncol(y))
  logsda <- rep(0,nrow(Ua))
  corsa <- rep(1e-5,nrow(csa))
  logsdo <- rep(0, nrow(Uo))
  corso <- 0
  
  # get some starting values
  # map formatted like this makes the model into u_i + b_j for eta^a and u_i for etaa^o
  # this should probably be interpreted via a formula interface in terms of an interaction with a species covariate
  # so that absence of an interaction gives a community-level effect
  # We could also do it without this, but this might speed things up a little later
  maplist = list(Ba = factor(c(rbind(1:ncol(y),replicate(ncol(y),1:(nrow(Ba)-1)+ncol(y))))), Bo=factor(rep(1:nrow(Ba),times=ncol(y))))
  maplist <- list() # just keeping map empty for now
  
  #TMB::openmp(parallel::detectCores()-1,autopar=TRUE, DLL = "occ")
  
  fit <- MakeADFun(data = list(Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = as(matrix(0), "TsparseMatrix"), Zo = as(matrix(0),"TsparseMatrix"), family = family,sites = sites, csa = matrix(0), cso = matrix(0), NTrials = NTrials, linka = linka, linko = linko), 
                   parameters = list(Ba=Ba, Bo=Bo, Ua=Ua, Uo=Uo, logphi=logphi, logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso),
                   DLL = "eDNAModel",
                   map = maplist)
  opt <- optim(fit$par, fit$fn, fit$gr, method="L-BFGS-B", control = list(trace = 1, maxit =5000)) # optimize
  
  #opt <- nlminb(start = fit$par, objective = fit$fn, gradient = fit$gr,
  # control = list(trace = 1, eval.max = 5000, iter.max = 5000))
  
  # Get the estimates in their correct format again
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
  
  # Final TMB call with integration over REs
  random = NULL
  if(nrow(Za)>1)random <- c(random, "Ua")
  if(nrow(Zo)>1)random <- c(random, "Uo")
  fit <- MakeADFun(data = list(Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = Za, Zo = Zo,family = family,sites = sites, csa = csa, cso = cso, NTrials = NTrials, linka = linka, linko = linko), 
                   parameters = list(Ba=Ba2, Bo=Bo2, Ua=Ua, Uo=Uo, logphi=logphi, logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso),
                   random = random,
                   DLL = "eDNAModel", silent = FALSE,
                   inner.control=list(mgcmax = 1e+200,tol10=0.01),
                   map = maplist)
  opt2 <- optim(fit$par, fit$fn, fit$gr, method="L-BFGS-B", control = list(trace = 1, maxit = 1e6)) # optimize
  
  ## End modeling ##
  
  ## Some results ##
  # make sure you get the link functions right here, I haven't bothered to code this up in detail
  occup.prob <- (1-plogis(fit$report(fit$env$last.par.best)$etao))
  #abundance mean
  lambda <- exp(fit$report(fit$env$last.par.best)$etaa)
  #prob. detect. at site
  prob.detect <- 1-occup.prob*as.matrix(aggregate(exp(-lambda), FUN=prod, list(sites+1)))[,-1]
  
  #let's check this is correct
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
