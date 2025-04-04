#' Run Full TMB Model
#'
#' Fits a hierarchical multi-species occupancy-abundance model using TMB,
#' with support for covariate formulas and random effects.
#'
#' @param data_array_filtered A 3D array (Species x Sites x Replicates) after filtering.
#' @param covariate_data A data frame containing covariates for model matrices.
#' @param a.formula A formula for the abundance model. Default is \code{~ site + (1 | replicate)}.
#' @param o.formula A formula for the occupancy model. Default is \code{~ site}.
#'
#' @return A list containing the fitted TMB object and optimization result.
#' @export
run_full_TMB <- function(data_array_filtered,
                         covariate_data,
                         a.formula = ~ site + (1 | replicate),
                         o.formula = ~ site) {

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

  Y<- to2D(data_array_filtered)

  # Inspect dimensions
  dim(data_array_filtered)
  # [1] n_species n_sites n_replicates

  n_sites <- dim(data_array_filtered)[2]
  n_replicates <- dim(data_array_filtered)[3]

  # Construct covariate_data to match row count of Y = n_sites * n_replicates
  covariate_data <- data.frame(
    site = factor(rep(1:n_sites, each = n_replicates)),
    replicate = factor(rep(1:n_replicates, times = n_sites))
  )



  Y<- to2D(data_array_filtered)
  y <- as.matrix(Y[, -(1:2)])
  x <- covariate_data

  a.formula = ~ site + (1 | replicate)
  o.formula = ~ site


  y <- as.matrix(Y[,-c(1:2)])
  sites = as.numeric(Y[,1])-1
  ysites <- as.matrix(aggregate(y,FUN=sum,list(sites)))[,-1]#sum over replicates
  Xa <- model.matrix(~ site, covariate_data)
  Zalist <- gllvm:::mkReTrms1(gllvm:::findbars1(a.formula), covariate_data)
  csa <- Zalist$cs # set to single column matrix if no REs
  Za <- Matrix::t(Zalist$Zt)
  Xo <- model.matrix(o.formula, covariate_data)
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
  #setwd("~/Documents/eDNAModel/eDNAModel")
  #TMB::compile("src/eDNAModel.cpp")
  #dyn.load(TMB::dynlib("src/eDNAModel"))



  fit <- TMB::MakeADFun(data = list(Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = as(matrix(0), "TsparseMatrix"), Zo = as(matrix(0),"TsparseMatrix"), family = family,sites = sites, csa = matrix(0), cso = matrix(0), NTrials = NTrials, linka = linka, linko = linko),
                        parameters = list(Ba=Ba, Bo=Bo, Ua=Ua, Uo=Uo, logphi=logphi, logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso),
                        DLL = "eDNAModel",
                        map = maplist)
  opt <- optim(fit$par, fit$fn, fit$gr, method="L-BFGS-B", control = list(trace = 1, maxit =5000)) # optimize

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

  opt2 <- optim(fit$par, fit$fn, fit$gr, method = "L-BFGS-B", control = list(trace = 1, maxit = 1e6))

  site_ids <- as.numeric(sites)
  etao <- fit$report(fit$env$last.par.best)$etao
  occup.prob <- 1 - plogis(etao)
  occup.prob_site <- aggregate(occup.prob, by = list(site = site_ids), FUN = mean)[, -1]

  lambda_matrix <- matrix(exp(fit$report(fit$env$last.par.best)$etaa), ncol = ncol(occup.prob))
  lambda_prod <- aggregate(exp(-lambda_matrix), by = list(site = site_ids), FUN = prod)[, -1]

  prob.detect <- 1 - occup.prob_site * lambda_prod

  hist(occup.prob, main = "occupancy probability", xlab="")
  hist(unlist(prob.detect), main = "Probability of Detection", xlab = "")

  #let's check this is correct
  #any(ysites[occup.prob[,]>0.8]==0) # NOPE, everything with large occup. prob. has counts larger than 0
  #all(ysites[occup.prob[,]<0.8]==0) # YUP, everything with small occup. prob. has count 0

  #any_violations <- any(ysites[occup.prob[,] > 0.8] == 0)
  #all_correct <- all(ysites[occup.prob[,] < 0.8] == 0)


  message(" Final TMB model fitting completed.")

  return(list(
    optimization = opt2,
    occupancy_probability = occup.prob,
    detection_probability = prob.detect
  ))

}


