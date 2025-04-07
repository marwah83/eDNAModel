#' Run Full TMB Model
#'
#' Fits a hierarchical multi-species occupancy-abundance model using TMB,
#' with support for covariate formulas and random effects.
#'
#' @param data_array_filtered A 3D array (Species x Sites x Replicates) after filtering.
#' @param X A data frame containing covariates for model matrices.
#' @param a.formula A formula for the abundance model. Default is \code{~ site + (1 | replicate)}.
#' @param o.formula A formula for the occupancy model. Default is \code{~ site}.
#'
#' @return A list containing the fitted TMB object and optimization result.
#' @export
run_full_TMB <- function(data_array_filtered,
                         X,
                         a.formula = ~ 1,
                         o.formula = ~ 1, linko = 1, linka = 0, family = 1, Ntrials = matrix(0), control = list(maxit = 10e3, trace = 1)) {
  
   
    Y<- to2D(data_array_filtered)
    
    # Inspect dimensions
    dim(data_array_filtered)
    # [1] n_species n_sites n_replicates
    
    n_sites <- dim(data_array_filtered)[2]
    n_replicates <- dim(data_array_filtered)[3]    
    
    Y<- to2D(data_array_filtered)
    y <- as.matrix(Y[, -(1:2)])
    
    sites = as.numeric(Y[,1])-1
    ysites <- as.matrix(aggregate(y,FUN=sum,list(sites)))[,-1]#sum over replicates
    Xa <- model.matrix(gllvm:::nobars1_(a.formula), X)
    if(gllvm:::anyBars(a.formula)){
    Zalist <- gllvm:::mkReTrms1(gllvm:::findbars1(a.formula), X)
    csa <- Zalist$cs # set to single column matrix if no REs
    Za <- Matrix::t(Zalist$Zt)
    }else{
    Za = as(matrix(0), "TsparseMatrix")
    csa = matrix(0)
    }
    
    # Here create a new Xocc at the site-level
    
    #
    Xo <- model.matrix(gllvm:::nobars1_(o.formula), Xocc) # This is a problem: here the  X should be at the site-level
    if(gllvm:::anyBars(o.formula)){
    Zolist <- gllvm:::mkReTrms1(gllvm:::findbars1(o.formula), Xocc)
    cso <- Zolist$cs # set to single column matrix if no REs
    Zo <- Matrix::t(Zolist$Zt)
    }else{
    Zo <- as(matrix(0), "TsparseMatrix")
    cso <- matrix(0)# set to single column matrix if no REs
    }
    
    # Add defensive programming for family still
    #family = 1#0 is ZIP, 1 is ZINB, 2 is BINOM with Ntrials
    # Add defensive coding for linka and linko here to check that linka is one of 0,1,2,3 and linko  is one of 1,2
    # linka = 0 # abundance link function; 0 = log, 1 = logit,  2 = probit, 3 = cloglog
    # linko = 1 # occupancy link function; 1 = logit,  2 = probit
    # add defensive coding on NTrials
    
    # Readying "parameters" component TMB list function
    Ba = matrix(0,nrow=ncol(Xa),ncol=ncol(y))
    Bo = matrix(0,nrow=ncol(Xo),ncol=ncol(y))
    Ua = matrix(0)
    Uo = matrix(0)
    logphi=rep(0,ncol(y))
    logsda <- rep(0,nrow(Ua))
    corsa <- rep(1e-5,nrow(csa))
    logsdo <- rep(0, nrow(Uo))
    corso <- rep(1e-5, nrow(cso))
    
    # get some starting values
    # map formatted like this makes the model into u_i + b_j for eta^a and u_i for etaa^o
    # this should probably be interpreted via a formula interface in terms of an interaction with a species covariate
    # so that absence of an interaction gives a community-level effect
    # We could also do it without this, but this might speed things up a little later
    maplist = list(Ba = factor(c(rbind(1:ncol(y),replicate(ncol(y),1:(nrow(Ba)-1)+ncol(y))))), Bo=factor(rep(1:nrow(Ba),times=ncol(y))))
    maplist <- list() # just keeping map empty for now
    
    #TMB::openmp(parallel::detectCores()-1,autopar=TRUE, DLL = "occ")
    #setwd("~/Documents/eDNAModel/eDNAModel")
    #TMB::compile("src/eDNAModel.cpp")
    #dyn.load(TMB::dynlib("src/eDNAModel"))
     
    fit <- TMB::MakeADFun(data = list(Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = as(matrix(0), "TsparseMatrix"), Zo = as(matrix(0),"TsparseMatrix"), family = family,sites = sites, csa = matrix(0), cso = matrix(0), NTrials = NTrials, linka = linka, linko = linko),
                          parameters = list(Ba=Ba, Bo=Bo, Ua=Ua, Uo=Uo, logphi=logphi, logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso),
                          DLL = "eDNAModel",
                          map = maplist)
    opt <- optim(fit$par, fit$fn, fit$gr, method="L-BFGS-B", control = list(trace = control$trace, maxit = 5e3)) # optimize
    
    Bas <- opt$par[names(opt$par) == "Ba"]
    Ba2 <- matrix(Bas, nrow = ncol(Xa))
    
    Bos <- opt$par[names(opt$par) == "Bo"]
    Bo2 <- matrix(Bos, nrow = ncol(Xo))
    
    if (nrow(Za) > 1) {
      Ua <- matrix(0, nrow = ncol(Za), ncol = ncol(y))
    } else {
      Ua <- matrix(0)
    }
    
    if (nrow(Zo) > 1) {
      Uo <- matrix(0, nrow = ncol(Zo), ncol = ncol(y))
    } else {
      Uo <- matrix(0)
    }
    
    logsda <- rep(0, nrow(Ua))
    corsa <- rep(1e-5, nrow(csa))
    logsdo <- rep(0, nrow(Uo))
    corso <- rep(1e-5, nrow(cso))
    logphi <- opt$par[names(opt$par) == "logphi"]
    
    random <- NULL
    if (nrow(Za) > 1) random <- c(random, "Ua")
    if (nrow(Zo) > 1) random <- c(random, "Uo")
    
    fit <- TMB::MakeADFun(
      data = list(Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = Za, Zo = Zo, family = family, sites = sites,
                  csa = csa, cso = cso, NTrials = NTrials, linka = linka, linko = linko),
      parameters = list(Ba = Ba2, Bo = Bo2, Ua = Ua, Uo = Uo, logphi = logphi,
                        logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso),
      random = random,
      DLL = "eDNAModel", silent = FALSE,
      inner.control = list(mgcmax = 1e+200, tol10 = 0.01),
      map = maplist
    )
    
    opt2 <- optim(fit$par, fit$fn, fit$gr, method = "L-BFGS-B", control = list(trace = control$trace, maxit = control$maxit))
    
    site_ids <- as.numeric(sites)
    etao <- fit$report(fit$env$last.par.best)$etao
    occup.prob <- 1 - plogis(etao)
    occup.prob_site <- aggregate(occup.prob, by = list(site = site_ids), FUN = mean)[, -1]
    
    lambda_matrix <- matrix(exp(fit$report(fit$env$last.par.best)$etaa), ncol = ncol(occup.prob))
    lambda_prod <- aggregate(exp(-lambda_matrix), by = list(site = site_ids), FUN = prod)[, -1]
    
    prob.detect <- 1 - occup.prob_site * lambda_prod
    
    #Separate plotting into its own function
    #hist(occup.prob, main = "occupancy probability", xlab="")
    #hist(unlist(prob.detect), main = "Probability of Detection", xlab = "")
    
    #let's check this is correct
    #any(ysites[occup.prob[,]>0.8]==0) # NOPE, everything with large occup. prob. has counts larger than 0
    #all(ysites[occup.prob[,]<0.8]==0) # YUP, everything with small occup. prob. has count 0
    
    #any_violations <- any(ysites[occup.prob[,] > 0.8] == 0)
    #all_correct <- all(ysites[occup.prob[,] < 0.8] == 0)

    return(list(start = opt, 
      optimization = opt2,
      occupancy_probability = occup.prob,
      detection_probability = prob.detect,
      TMBobj = fit
    ))
    
  }
  
  
