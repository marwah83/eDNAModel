#' Run Full TMB Model
#'
#' Fits a hierarchical multi-species occupancy-abundance model using TMB,
#' with support for covariate formulas and random effects.
#'
#' @param data_array_filtered A 3D array (Species x Sites x Replicates) after filtering.
#' @param X A data frame containing the covariates for model matrices.
#' @param a.formula A formula for the abundance model. Default is \code{~ 1}.
#' @param o.formula A formula for the occupancy model. Default is \code{~ 1}.
#' @param linko Link function for occupancy. 1 = logit, 2 = probit.
#' @param linka Link function for abundance. 0 = log, 1 = logit, 2 = probit, 3 = cloglog.
#' @param family Family of distribution. 0 = ZIP, 1 = ZINB, 2 = Binomial.
#' @param Ntrials Trials matrix for Binomial. Default is matrix(0).
#' @param control List of optimization controls. Default: list(maxit = 10e3, trace = 1).
#'
#' @return A list containing the fitted TMB object and optimization result.
#' @export
run_full_TMB <- function(data_array_filtered,
                         X=NULL,
                         a.formula = ~ 1,
                         o.formula = ~ 1,
                         linko = 1,
                         linka = 0,
                         family = 1,
                         Ntrials = matrix(0),
                         control = list(maxit = 10e3, trace = 1)) {
  
  Y <- to2D(data_array_filtered)
  y <- as.matrix(Y[, -(1:2)])
  sites <- as.numeric(Y[, 1]) - 1
  ysites <- as.matrix(aggregate(y, FUN = sum, list(sites)))[, -1]
  
  X<- data.frame(
    site = factor(rep(1:3, each = 4)),
    replicate = factor(rep(1:4, times = 3))
  )
  
  
  if (!linko %in% c(1, 2)) stop("linko must be 1 (logit) or 2 (probit)")
  if (!linka %in% c(0, 1, 2, 3)) stop("linka must be 0 (log), 1 (logit), 2 (probit), or 3 (cloglog)")
  if (!family %in% c(0, 1, 2)) stop("family must be 0 (ZIP), 1 (ZINB), or 2 (Binomial)")
  if (!is.matrix(Ntrials)) stop("Ntrials must be a matrix")
  
  Xa <- model.matrix(gllvm:::nobars1_(a.formula), X)
  if (gllvm:::anyBars(a.formula)) {
    Zalist <- gllvm:::mkReTrms1(gllvm:::findbars1(a.formula), X)
    Za <- Matrix::t(Zalist$Zt)
    csa <- Zalist$cs
  } else {
    Za <- as(matrix(0), "TsparseMatrix")
    csa <- matrix(0)
  }
  
  
  Xocc <- X[!duplicated(sites + 1), , drop = FALSE]
  Xo <- model.matrix(gllvm:::nobars1_(o.formula), Xocc)
  if(nrow(Xo)!=length(unique(sites)))stop("Something went horribly wrong!")
  
  if (gllvm:::anyBars(o.formula)) {
    Zolist <- gllvm:::mkReTrms1(gllvm:::findbars1(o.formula), Xocc)
    Zo <- Matrix::t(Zolist$Zt)
    cso <- Zolist$cs
  } else {
    Zo <- as(matrix(0), "TsparseMatrix")
    cso <- matrix(0)
  }
  
  Ba <- matrix(0, nrow = ncol(Xa), ncol = ncol(y))
  Bo <- matrix(0, nrow = ncol(Xo), ncol = ncol(y))
  Ua <- matrix(0)
  Uo <- matrix(0)
  logphi <- rep(0, ncol(y))
  logsda <- rep(0, nrow(Ua))
  corsa <- rep(1e-5, nrow(csa))
  logsdo <- rep(0, nrow(Uo))
  corso <- rep(1e-5, nrow(cso))
  
  maplist <- list()
  
  fit <- TMB::MakeADFun(
    data = list(
      Y = y, Ysites = ysites, Xa = Xa, Xo = Xo,
      Za = Za, Zo = Zo, family = family, sites = sites,
      csa = csa, cso = cso, NTrials = Ntrials,
      linka = linka, linko = linko
    ),
    parameters = list(
      Ba = Ba, Bo = Bo, Ua = Ua, Uo = Uo,
      logphi = logphi, logsda = logsda, corsa = corsa,
      logsdo = logsdo, corso = corso
    ),
    DLL = "eDNAModel",
    map = maplist
  )
  
  opt <- optim(fit$par, fit$fn, fit$gr, method = "L-BFGS-B", control = list(trace = control$trace, maxit = control$maxit))
  
  return(list(fit = fit, opt = opt))


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

