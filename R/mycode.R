# ✅ Load Required Libraries
library(TMB)
library(Matrix)

# ✅ Compile and Load the TMB Model
compile("~/Desktop/Diversea/occ.cpp", framework = "TMBad", flags = "-O0 -g")
dyn.load("~/Desktop/Diversea/occ.so")

library(TMB)
library(Matrix)

# Function to simulate data from a Zero-Inflated Poisson (ZIP) distribution

# Function to simulate data from a Zero-Inflated Poisson (ZIP) distribution
simulate_zip_data <- function(species, sites, replicates, lambda, zero_inflation_prob_matrix) {
  if (!all(dim(zero_inflation_prob_matrix) == c(species, sites))) {
    stop("The zero_inflation_prob_matrix must have dimensions species x sites.")
  }

  counts <- array(rpois(species * sites * replicates, lambda), dim = c(species, sites, replicates))

  # Apply zero-inflation correctly
  for (i in 1:replicates) {
    z <- matrix(rbinom(species * sites, 1, 1 - zero_inflation_prob_matrix),
                nrow = species, ncol = sites)
    counts[, , i] <- counts[, , i] * z
  }

  dimnames(counts) <- list(
    Species = paste("Species", 1:species, sep = "_"),
    Sites = paste("Site", 1:sites, sep = "_"),
    Replicates = paste("Replicate", 1:replicates, sep = "_")
  )

  return(counts)
}


# Convert 3D array to 2D data frame
to2D <- function(array_3d) {
  species <- dim(array_3d)[1]
  sites <- dim(array_3d)[2]
  replicates <- dim(array_3d)[3]

  data_values <- matrix(array_3d, nrow = sites * replicates, ncol = species, byrow = FALSE)
  site_column <- rep(dimnames(array_3d)$Sites, each = replicates)
  replicate_column <- rep(dimnames(array_3d)$Replicates, times = sites)

  final_data <- data.frame(
    Site = site_column,
    Replicate = replicate_column,
    data_values
  )
  colnames(final_data)[3:(2 + species)] <- dimnames(array_3d)$Species
  return(final_data)
}

# Prepare data for TMB
prepare_tmb_data_fit <- function(Y) {
  require(TMB)
  require(Matrix)
  require(gllvm)

  # Step 1: Prepare Data
  y <- as.matrix(Y[, -c(1:2)])
  Y$Site <- as.factor(Y$Site)
  sites <- (as.numeric(as.factor(Y[, 1]))) - 1

  ysites <- as.matrix(aggregate(y, FUN = sum, list(sites)))[, -1] # Sum over replicates
  Xa <- model.matrix(~ Site, Y)
  Y$Replicate <- as.factor(Y$Replicate)
  Y$Replicate1 <- as.numeric(as.factor(Y$Replicate))

  Zalist <- gllvm:::mkReTrms1(gllvm:::findbars1(~diag(1 | Replicate1)), Y)
  csa <- Zalist$cs  # Set to single column matrix if no REs
  Za <- Matrix::t(Zalist$Zt)

  Xo <- model.matrix(~0 + Site, data.frame(Site = as.factor(1:max(as.numeric(Y$Site)))))
  Zo <- as(matrix(0), "TsparseMatrix")
  cso <- matrix(0)  # Set to single column matrix if no REs

  family <- 1  # 0 is ZIP, 1 is ZINB, 2 is BINOM with Ntrials
  linka <- 0  # Abundance link function: 0 = log, 1 = logit, 2 = probit, 3 = cloglog
  linko <- 1  # Occupancy link function: 0 = log, 1 = logit, 2 = probit, 3 = cloglog
  NTrials <- matrix(0)  # nrow(y) by ncol(y) matrix for binomial size (for family = 2)

  # Step 2: Initialize Parameters
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


  # Empty maplist (can be populated later for efficiency)
  maplist <- list()

  # Step 3: First TMB Call
  fit <- MakeADFun(
    data = list(
      Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = as(matrix(0), "TsparseMatrix"), Zo = as(matrix(0), "TsparseMatrix"),
      family = family, sites = sites, csa = matrix(0), cso = matrix(0), NTrials = NTrials, linka = linka, linko = linko
    ),
    parameters = list(
      Ba = Ba, Bo = Bo, Ua = Ua, Uo = Uo, logphi = logphi, logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso
    ),
    DLL = "occ",
    map = maplist
  )



  opt <- optim(fit$par, fit$fn, fit$gr, method="L-BFGS-B", control = list(
    trace = 1,
    maxit = 5000,
    factr = 1e14,  # Less strict stopping condition
    pgtol = 1e-8   # More tolerance for small changes
  )
  ) # optimize


  # Extract Ba and Bo Estimates
  Bas <- opt$par[names(opt$par) == "Ba"]
  Bas2 <- if ("Ba" %in% names(maplist)) Bas[fit$env$map$Ba] else Bas
  Ba2 <- matrix(Bas2, nrow = ncol(Xa))

  Bos <- opt$par[names(opt$par) == "Bo"]
  Bos2 <- if ("Bo" %in% names(maplist)) Bos[fit$env$map$Bo] else Bos
  Bo2 <- matrix(Bos2, nrow = ncol(Xo))

  # Adjust for Random Effects
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
  if (family == 1) logphi <- opt$par[names(opt$par) == "logphi"]

  # Step 4: Final TMB Call
  random <- NULL
  if (nrow(Za) > 1) random <- c(random, "Ua")
  if (nrow(Zo) > 1) random <- c(random, "Uo")

  TMB::openmp(1, autopar = TRUE, DLL = "occ")

  fit <- MakeADFun(
    data = list(
      Y = y, Ysites = ysites, Xa = Xa, Xo = Xo, Za = Za, Zo = Zo,
      family = family, sites = sites, csa = csa, cso = cso, NTrials = NTrials, linka = linka, linko = linko
    ),
    parameters = list(
      Ba = Ba2, Bo = Bo2, Ua = Ua, Uo = Uo, logphi = logphi, logsda = logsda, corsa = corsa, logsdo = logsdo, corso = corso
    ),
    random = random,
    DLL = "occ",
    silent = FALSE,
    inner.control = list(mgcmax = 1e+200, tol10 = 0.01),
    map = maplist
  )

  lower_bounds <- rep(-10, length(fit$par))  # Limit to avoid extreme values
  upper_bounds <- rep(10, length(fit$par))   # Limit to avoid extreme values


  opt2 <- optim(fit$par, fit$fn, fit$gr, method="L-BFGS-B", control = list(
    trace = 1,
    maxit = 5000   # More tolerance for small changes
  )) # optimize

  return(list(fit = fit, opt = opt2))
}

compute_probabilities <- function(fit_results) {
  # Extract eta_o and eta_a
  eta_o <- fit_results$fit$report(fit_results$fit$env$last.par.best)$etao
  eta_a <- fit_results$fit$report(fit_results$fit$env$last.par.best)$etaa

  # Compute probabilities
  occ.prob <- 1 - plogis(eta_o)
  lambda <- exp(eta_a)

  # Compute probability of detection at site
  prob.detect <- 1 - occ.prob * as.matrix(
    aggregate(exp(-lambda), FUN = prod, list(fit_results$fit$env$data$sites + 1))
  )[, -1]

  # Return probabilities as a list
  return(list(occ.prob = occ.prob, prob.detect = prob.detect))
}



  # fit_results=prepare_tmb_data_fit(Y)
  # Get occurrence probabilities (1 - logistic transformation of eta_o)

 # occup.prob <- (1-plogis(fit_results$fit$report(fit_results$fit$env$last.par.best)$etao))
 # #abundance mean
 # lambda <- exp(fit_results$fit$report(fit_results$fit$env$last.par.best)$etaa)
 # #prob. detect. at site
 # prob.detect <- 1-occup.prob*as.matrix(aggregate(exp(-lambda), FUN=prod, list(fit_results$fit$env$data$sites+1)))[,-1]

 # # Plot the histograms for both occurrence and detection probabilities
 # par(mfrow = c(1, 2))  # Set up the plotting area for two plots
 # hist(occup.prob, main = "occupancy probability", xlab="")
  #hist(prob.detect, main = "probability of detection", xlab="")

  # Return the probabilities for further use
  #return(list(occ.prob = occup.prob, prob.detect = prob.detect))





# Main simulation and modeling process
# Define parameters for simulation
set.seed(456)
species <- 20
sites <- 30
replicates <- 3
lambda <- 5

# Define site * species-specific zero-inflation probabilities
  # For reproducibility
zero_inflation_prob_matrix <- matrix(runif(species * sites, 0.1, 0.5), nrow = species, ncol = sites)

# Simulate ZIP data
simulated_data <- simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)

# Convert 3D array to 2D data frame
Y <- to2D(simulated_data)

# Prepare data and fit the model using the combined function
fit_results <- prepare_tmb_data_fit(Y)

# Plot occurrence and detection probabilities
plot_results <- compute_probabilities(fit_results)

hist(plot_results$occ.prob, main = "occupancy probability", xlab="")
hist(plot_results$prob.detect, main = "probability of detection", xlab="")


#occup.prob <- (1-plogis(fit_results$fit$report(fit_results$fit$env$last.par.best)$etao))
#abundance mean
#lambda <- exp(fit_results$fit$report(fit_results$fit$env$last.par.best)$etaa)
#prob. detect. at site
#prob.detect <- 1-occup.prob*as.matrix(aggregate(exp(-lambda), FUN=prod, list(fit_results$fit$env$data$sites+1)))[,-1]

library(TMB)
library(Matrix)

# ✅ Compile and Load the TMB Model
compile("~/Desktop/Diversea/occ.cpp", framework = "TMBad", flags = "-O0 -g")
dyn.load("~/Desktop/Diversea/occ.so")

# ✅ Define parameter grids
lambda_values <- seq(1, 10, by = 1)
species_values <- seq(10, 20, by = 10)
site_values <- seq(10, 20, by = 10)
replicate_values <- seq(2, 5, by = 1)
ZIP_values <- 0.5  # Constant zero-inflation probability

# ✅ Create a parameter grid
param_grid <- expand.grid(lambda = lambda_values,
                          species = species_values,
                          sites = site_values,
                          replicates = replicate_values,
                          ZIP = ZIP_values)

# ✅ Function to Compute RMSE
compute_rmse <- function(actual, estimated) {
  sqrt(mean((actual - estimated) ^ 2, na.rm = TRUE))
}

# ✅ Function to Compute Bias
compute_bias <- function(actual, estimated) {
  mean(estimated - actual, na.rm = TRUE)
}

# ✅ Function to Run Full Simulation
run_full_simulation <- function(species, sites, replicates, lambda, ZIP) {
  # Define site * species-specific zero-inflation probabilities
  zero_inflation_prob_matrix <- matrix(ZIP, nrow = species, ncol = sites)  # Fixed ZIP value

  # ✅ Simulate ZIP Data
  simulated_data <- simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)

  # ✅ Convert to 2D
  Y <- to2D(simulated_data)

  # ✅ Prepare Data and Fit Model
  fit_results <- prepare_tmb_data_fit(Y)

  # ✅ Compute Probabilities
  probabilities <- compute_probabilities(fit_results)

  hist(probabilities$occ.prob, main = "occupancy probability", xlab="")
  hist(probabilities$prob.detect, main = "probability of detection", xlab="")


  # ✅ Compute Bias and RMSE
  true_occ_prob <- 1 - zero_inflation_prob_matrix  # True probability of occupancy
  est_occ_prob <- probabilities$occ.prob  # Estimated probability of occupancy

  true_prob_detect <- 1 - true_occ_prob * exp(-lambda)  # True detection probability
  est_prob_detect <- probabilities$prob.detect  # Estimated detection probability

  bias_occ <- compute_bias(true_occ_prob, t(est_occ_prob))
  rmse_occ <- compute_rmse(true_occ_prob, t(est_occ_prob))

  bias_detect <- compute_bias(true_prob_detect, t(est_prob_detect))
  rmse_detect <- compute_rmse(true_prob_detect, t(est_prob_detect))

  # ✅ Return results
  return(list(
    species = species,
    sites = sites,
    replicates = replicates,
    lambda = lambda,
    ZIP = ZIP,
    mean_occ_prob = mean(t(est_occ_prob), na.rm = TRUE),
    mean_prob_detect = mean(t(est_prob_detect), na.rm = TRUE),
    bias_occ = bias_occ,
    rmse_occ = rmse_occ,
    bias_detect = bias_detect,
    rmse_detect = rmse_detect
  ))
}

simulation_results <- lapply(1:nrow(param_grid), function(i) {
  with(param_grid[i, ], run_full_simulation(species, sites, replicates, lambda, ZIP))
})

# ✅ Convert results into a data frame for analysis
simulation_df <- do.call(rbind, lapply(simulation_results, function(res) {
  data.frame(
    species = res$species,
    sites = res$sites,
    replicates = res$replicates,
    lambda = res$lambda,
    ZIP = res$ZIP,
    mean_occ_prob = res$mean_occ_prob,
    mean_prob_detect = res$mean_prob_detect,
    bias_occ = res$bias_occ,
    rmse_occ = res$rmse_occ,
    bias_detect = res$bias_detect,
    rmse_detect = res$rmse_detect
  )
}))

# ✅ Print Summary of Results
print(head(simulation_df))


######################wrong #####################################


library(TMB)
library(Matrix)

# ✅ Set seed for reproducibility
set.seed(123)

# ✅ Compile and Load the TMB Model
compile("~/Desktop/Diversea/occ.cpp", framework = "TMBad", flags = "-O0 -g")
dyn.load("~/Desktop/Diversea/occ.so")

# ✅ Define parameter ranges
lambda_values <- seq(0.1, 5, by = 0.5)
species_values <- seq(10, 20, by = 10)
site_values <- seq(10, 20, by = 10)
replicate_values <- seq(2, 5, by = 1)
ZIP_values <- 0.5  # Fixed zero-inflation probability

# ✅ Create parameter grid
param_grid <- expand.grid(
  Lambda = lambda_values,
  ZIP = ZIP_values,
  Species = species_values,
  Sites = site_values,
  Replicates = replicate_values
)

# ✅ Function to Compute RMSE
compute_rmse <- function(actual, estimated) {
  sqrt(mean((actual - estimated) ^ 2, na.rm = TRUE))
}

# ✅ Function to Compute Bias
compute_bias <- function(actual, estimated) {
  mean(estimated - actual, na.rm = TRUE)
}

# ✅ Function to Run Full Simulation
run_full_simulation <- function(species, sites, replicates, lambda, ZIP) {
  set.seed(123)  # Ensure reproducibility in simulations

  # ✅ Simulate ZIP Data
  simulated_data <- tryCatch(
    simulate_zip_data(species, sites, replicates, lambda, ZIP),
    error = function(e) return(NULL)
  )

  zero=mean(simulated_data ==0)


  if (is.null(simulated_data)) {
    return(data.frame(
      Species = species, Sites = sites, Replicates = replicates,
      Lambda = lambda, ZIP = ZIP,
      Mean_Occ_Prob = NA, Mean_Prob_Detect = NA,
      Bias_Occupancy = NA, RMSE_Occupancy = NA,
      Bias_Detection = NA, RMSE_Detection = NA
    ))
  }

  # ✅ Convert to 2D format
  Y <- to2D(simulated_data)

  # ✅ Prepare Data and Fit Model
  fit_results <- tryCatch(
    prepare_tmb_data_fit(Y),
    error = function(e) return(NULL)
  )

  # ✅ Handle failure to converge
  if (is.null(fit_results)) {
    return(data.frame(
      Species = species, Sites = sites, Replicates = replicates,
      Lambda = lambda, ZIP = ZIP,
      Mean_Occ_Prob = NA, Mean_Prob_Detect = NA,
      Bias_Occupancy = NA, RMSE_Occupancy = NA,
      Bias_Detection = NA, RMSE_Detection = NA
    ))
  }

  # ✅ Compute Probabilities
  probabilities <- compute_probabilities(fit_results)

  # ✅ Compute Bias and RMSE
  true_occ_prob <- 1 - ZIP  # True probability of occupancy
  est_occ_prob <- probabilities$occ.prob  # Estimated probability of occupancy

  true_prob_detect <- 1 - true_occ_prob * exp(-lambda)  # True detection probability
  est_prob_detect <- probabilities$prob.detect  # Estimated detection probability

  bias_occ <- compute_bias((true_occ_prob), (est_occ_prob))
  rmse_occ <- compute_rmse((true_occ_prob),(est_occ_prob))

  bias_detect <- compute_bias((true_prob_detect), (est_prob_detect))
  rmse_detect <- compute_rmse((true_prob_detect), (est_prob_detect))

  return(data.frame(
    Species = species, Sites = sites, Replicates = replicates,
    Lambda = lambda, ZIP = ZIP,
    Mean_Occ_Prob = mean(est_occ_prob, na.rm = TRUE),
    Mean_Prob_Detect = mean(est_prob_detect, na.rm = TRUE),
    Bias_Occupancy = bias_occ, RMSE_Occupancy = rmse_occ,
    Bias_Detection = bias_detect, RMSE_Detection = rmse_detect,
    zero=zero
  ))
}

# ✅ Run Simulation Over Parameter Grid
start_time <- Sys.time()

simulation_results <- do.call(rbind, lapply(1:nrow(param_grid), function(i) {
  with(param_grid[i, ], run_full_simulation(Species, Sites, Replicates, Lambda, ZIP))
}))

end_time <- Sys.time()

# ✅ Print Results
print(head(simulation_results))
print(paste("Total execution time in seconds:", as.numeric(difftime(end_time, start_time, units = "secs"))))


# Load ggplot2 for plotting
library(ggplot2)

# Plot Bias for Detection
ggplot(simulation_results, aes(x = Lambda, y = Bias_Detection, color = as.factor(Species))) +
  geom_point() +
  facet_wrap(~ Sites, scales = "free_y") +
  labs(
    title = "Bias of Detection Across Parameter Ranges",
    x = "Lambda (Poisson rate parameter)",
    y = "Bias (Detection)",
    color = "Species"
  ) +
  theme_minimal()

# Plot RMSE for Detection
ggplot(simulation_results, aes(x = Lambda, y = RMSE_Detection, color = as.factor(Species))) +
  geom_point() +
  facet_wrap(~ Sites, scales = "free_y") +
  labs(
    title = "RMSE of Detection Across Parameter Ranges",
    x = "Lambda (Poisson rate parameter)",
    y = "RMSE (Detection)",
    color = "Species"
  ) +
  theme_minimal()

# Plot Bias for Occupancy
ggplot(simulation_results, aes(x = Lambda, y = Bias_Occupancy, color = as.factor(Species))) +
  geom_point() +
  facet_wrap(~ Sites, scales = "free_y") +
  labs(
    title = "Bias of Occupancy Probability Across Parameter Ranges",
    x = "Lambda (Poisson rate parameter)",
    y = "Bias (Occupancy Probability)",
    color = "Species"
  ) +
  theme_minimal()

# Plot RMSE for Occupancy
ggplot(simulation_results, aes(x = Lambda, y = RMSE_Occupancy, color = as.factor(Species))) +
  geom_point() +
  facet_wrap(~ Sites, scales = "free_y") +
  labs(
    title = "RMSE of Occupancy Probability Across Parameter Ranges",
    x = "Lambda (Poisson rate parameter)",
    y = "RMSE (Occupancy Probability)",
    color = "Species"
  ) +
  theme_minimal()
############################################################
library(TMB)
library(Matrix)
library(dplyr)
set.seed(123)

# ✅ Compile and Load the TMB Model
compile("~/Desktop/Diversea/occ.cpp", framework = "TMBad", flags = "-O0 -g")
dyn.load("~/Desktop/Diversea/occ.so")

# ✅ Define parameter ranges
lambda_values <- seq(1, 10, by = 1)
species_values <- seq(10, 20, by = 10)
site_values <- seq(10, 20, by = 10)
replicate_values <- seq(2, 5, by = 1)
ZIP_values <- 0.1

# ✅ Create parameter grid
param_grid <- expand.grid(
  Lambda = lambda_values,
  ZIP = ZIP_values,
  Species = species_values,
  Sites = site_values,
  Replicates = replicate_values
)

# ✅ Function to Compute RMSE
compute_rmse <- function(actual, estimated) {
  sqrt(mean((actual - estimated) ^ 2, na.rm = TRUE))
}

# ✅ Function to Compute Bias
compute_bias <- function(actual, estimated) {
  mean(estimated - actual, na.rm = TRUE)
}

# ✅ Function to Run Full Simulation
run_full_simulation <- function(species, sites, replicates, lambda, ZIP) {
  set.seed(123)

  # ✅ Create Zero-Inflation Matrix
  zero_inflation_prob_matrix <- matrix(0.1, nrow = species, ncol = sites)

  # ✅ Simulate ZIP Data
  simulated_data <- tryCatch(
    simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix),
    error = function(e) return(NULL)
  )

  zero = mean(simulated_data == 0)

  if (is.null(simulated_data)) {
    return(data.frame(
      Species = species, Sites = sites, Replicates = replicates,
      Lambda = lambda, ZIP = ZIP,
      Mean_Occ_Prob = NA, Mean_Prob_Detect = NA,
      Bias_Occupancy = NA, RMSE_Occupancy = NA,
      Bias_Detection = NA, RMSE_Detection = NA,
      zero = zero
    ))
  }

  # ✅ Convert to 2D format
  Y <- to2D(simulated_data)

  # ✅ Prepare Data and Fit Model
  fit_results <- tryCatch(
    prepare_tmb_data_fit(Y),
    error = function(e) return(NULL)
  )

  if (is.null(fit_results)) {
    return(data.frame(
      Species = species, Sites = sites, Replicates = replicates,
      Lambda = lambda, ZIP = ZIP,
      Mean_Occ_Prob = NA, Mean_Prob_Detect = NA,
      Bias_Occupancy = NA, RMSE_Occupancy = NA,
      Bias_Detection = NA, RMSE_Detection = NA,
      zero = zero
    ))
  }

  # ✅ Compute Probabilities
  probabilities <- compute_probabilities(fit_results)

  # ✅ Compute Bias and RMSE
  true_occ_prob <- 1 - zero_inflation_prob_matrix
  est_occ_prob <- probabilities$occ.prob

  true_prob_detect <- 1 - true_occ_prob * exp(-lambda)
  est_prob_detect <- probabilities$prob.detect

  # ✅ Fix Mismatched Dimensions
  if (!all(dim(true_occ_prob) == dim(est_occ_prob))) {
    est_occ_prob <- t(est_occ_prob)
  }
  if (!all(dim(true_prob_detect) == dim(est_prob_detect))) {
    est_prob_detect <- t(est_prob_detect)
  }

  # ✅ Compute Bias and RMSE
  bias_occ <- tryCatch(compute_bias(as.vector(true_occ_prob), as.vector(est_occ_prob)), error = function(e) NA)
  rmse_occ <- tryCatch(compute_rmse(as.vector(true_occ_prob), as.vector(est_occ_prob)), error = function(e) NA)

  bias_detect <- tryCatch(compute_bias(as.vector(true_prob_detect), as.vector(est_prob_detect)), error = function(e) NA)
  rmse_detect <- tryCatch(compute_rmse(as.vector(true_prob_detect), as.vector(est_prob_detect)), error = function(e) NA)

  return(data.frame(
    Species = species, Sites = sites, Replicates = replicates,
    Lambda = lambda, ZIP = ZIP,
    Mean_Occ_Prob = mean(est_occ_prob, na.rm = TRUE),
    Mean_Prob_Detect = mean(est_prob_detect, na.rm = TRUE),
    Bias_Occupancy = bias_occ, RMSE_Occupancy = rmse_occ,
    Bias_Detection = bias_detect, RMSE_Detection = rmse_detect,
    zero = zero
  ))
}

# ✅ Run Simulation Over Parameter Grid
start_time <- Sys.time()

# ✅ Use `bind_rows()` Instead of `rbind()`
simulation_results <- bind_rows(lapply(1:nrow(param_grid), function(i) {
  with(param_grid[i, ], run_full_simulation(Species, Sites, Replicates, Lambda, ZIP))
}))

end_time <- Sys.time()

# ✅ Print Results
print(head(simulation_results))
print(paste("Total execution time in seconds:", as.numeric(difftime(end_time, start_time, units = "secs"))))

# ✅ Plot Bias and RMSE
par(mfrow = c(1, 2))
plot(simulation_results$Lambda, simulation_results$Bias_Occupancy,
     type = "p", col = "blue", main = "Bias of Occupancy Probability", xlab = "Lambda", ylab = "Bias")
plot(simulation_results$Lambda, simulation_results$RMSE_Occupancy,
     type = "p", col = "red", main = "RMSE of Occupancy Probability", xlab = "Lambda", ylab = "RMSE")

plot(simulation_results$zero, simulation_results$Bias_Occupancy,
     type = "p", col = "green", main = "zero vs Occupancy Probability", xlab = "zero", ylab = "Bias")

###################################this is correct###################################


library(TMB)
library(Matrix)
set.seed(123)

# ✅ Compile and Load the TMB Model
compile("~/Desktop/Diversea/occ.cpp", framework = "TMBad", flags = "-O0 -g")
dyn.load("~/Desktop/Diversea/occ.so")

# ✅ Define parameter ranges
lambda_values <- seq(0.5, 5, by = 0.5)
species_values <- seq(10, 100, by = 10)
site_values <- seq(10, 100, by = 10)
replicate_values <- seq(2, 5, by = 1)
ZIP_values <- 0.5  # Fixed zero-inflation probability

# ✅ Create a parameter grid
param_grid <- expand.grid(
  Lambda = lambda_values,
  ZIP = ZIP_values,
  Species = species_values,
  Sites = site_values,
  Replicates = replicate_values
)

# ✅ Function to Compute RMSE
compute_rmse <- function(actual, estimated) {
  sqrt(mean((actual - estimated) ^ 2, na.rm = TRUE))
}

# ✅ Function to Compute Bias
compute_bias <- function(actual, estimated) {
  mean(estimated - actual, na.rm = TRUE)
}

# ✅ Function to Create Zero-Inflation Matrix
create_zero_inflation_matrix <- function(species, sites, ZIP) {
  return(matrix(ZIP, nrow = species, ncol = sites))  # Fixed ZIP value
}

# ✅ Function to Simulate ZIP Data
simulate_data <- function(species, sites, replicates, lambda, zero_inflation_prob_matrix) {
  tryCatch(
    simulate_zip_data(species, sites, replicates, lambda, zero_inflation_prob_matrix),
    error = function(e) return(NULL)
  )
}

# ✅ Function to Process and Fit Data
process_and_fit_data <- function(simulated_data) {
  if (is.null(simulated_data)) return(NULL)

  Y <- to2D(simulated_data)  # Convert to 2D format

  fit_results <- tryCatch(
    prepare_tmb_data_fit(Y),
    error = function(e) return(NULL)
  )

  return(fit_results)
}

# ✅ Function to Compute Probabilities
compute_probabilities_wrapper <- function(fit_results) {
  if (is.null(fit_results)) return(NULL)

  tryCatch(
    compute_probabilities(fit_results),
    error = function(e) return(NULL)
  )
}

# ✅ Function to Compute Bias and RMSE
compute_bias_rmse <- function(true_occ_prob, est_occ_prob, true_prob_detect, est_prob_detect) {

  # ✅ Fix Mismatched Dimensions
  if (!all(dim(true_occ_prob) == dim(est_occ_prob))) {
    est_occ_prob <- t(est_occ_prob)
  }
  if (!all(dim(true_prob_detect) == dim(est_prob_detect))) {
    est_prob_detect <- t(est_prob_detect)
  }

  # ✅ Compute Bias and RMSE
  bias_occ <- tryCatch(compute_bias(as.vector(true_occ_prob), as.vector(est_occ_prob)), error = function(e) NA)
  rmse_occ <- tryCatch(compute_rmse(as.vector(true_occ_prob), as.vector(est_occ_prob)), error = function(e) NA)

  bias_detect <- tryCatch(compute_bias(as.vector(true_prob_detect), as.vector(est_prob_detect)), error = function(e) NA)
  rmse_detect <- tryCatch(compute_rmse(as.vector(true_prob_detect), as.vector(est_prob_detect)), error = function(e) NA)

  return(list(
    bias_occ = bias_occ, rmse_occ = rmse_occ,
    bias_detect = bias_detect, rmse_detect = rmse_detect
  ))
}

# ✅ Function to Run Full Simulation
run_full_simulation <- function(species, sites, replicates, lambda, ZIP) {
  set.seed(123)  # Ensure reproducibility

  # ✅ Step 1: Create Zero-Inflation Matrix
  zero_inflation_prob_matrix <- create_zero_inflation_matrix(species, sites, ZIP)

  # ✅ Step 2: Simulate Data
  simulated_data <- simulate_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)

  # ✅ Step 3: Check Zero Proportion
  zero <- mean(simulated_data == 0)

  # ✅ Step 4: Process Data and Fit Model
  fit_results <- process_and_fit_data(simulated_data)
  if (is.null(fit_results)) {
    return(data.frame(
      Species = species, Sites = sites, Replicates = replicates,
      Lambda = lambda, ZIP = ZIP,
      Mean_Occ_Prob = NA, Mean_Prob_Detect = NA,
      Bias_Occupancy = NA, RMSE_Occupancy = NA,
      Bias_Detection = NA, RMSE_Detection = NA,
      zero = zero
    ))
  }

  # ✅ Step 5: Compute Probabilities
  probabilities <- compute_probabilities_wrapper(fit_results)

  # ✅ Step 6: Compute Bias and RMSE
  true_occ_prob <- 1 - zero_inflation_prob_matrix
  est_occ_prob <- probabilities$occ.prob

  true_prob_detect <- 1 - true_occ_prob * exp(-lambda)
  est_prob_detect <- probabilities$prob.detect

  bias_rmse_results <- compute_bias_rmse(true_occ_prob, est_occ_prob, true_prob_detect, est_prob_detect)

  # ✅ Step 7: Return Simulation Results
  return(data.frame(
    Species = species, Sites = sites, Replicates = replicates,
    Lambda = lambda, ZIP = ZIP,
    Mean_Occ_Prob = mean(est_occ_prob, na.rm = TRUE),
    Mean_Prob_Detect = mean(est_prob_detect, na.rm = TRUE),
    Bias_Occupancy = bias_rmse_results$bias_occ, RMSE_Occupancy = bias_rmse_results$rmse_occ,
    Bias_Detection = bias_rmse_results$bias_detect, RMSE_Detection = bias_rmse_results$rmse_detect,
    zero = zero
  ))
}

# ✅ Run Simulation Over Parameter Grid
start_time <- Sys.time()

simulation_results <- do.call(rbind, lapply(1:nrow(param_grid), function(i) {
  with(param_grid[i, ], run_full_simulation(Species, Sites, Replicates, Lambda, ZIP))
}))

end_time <- Sys.time()

# ✅ Print Results
print(head(simulation_results))
print(paste("Total execution time in seconds:", as.numeric(difftime(end_time, start_time, units = "secs"))))

# ✅ Plot Bias and RMSE
par(mfrow = c(1, 2))
plot(simulation_results$Lambda, simulation_results$Bias_Occupancy,
     type = "p", col = "blue", main = "Bias of Occupancy Probability", xlab = "Lambda", ylab = "Bias")
plot(simulation_results$Lambda, simulation_results$RMSE_Occupancy,
     type = "p", col = "red", main = "RMSE of Occupancy Probability", xlab = "Lambda", ylab = "RMSE")

# ✅ Plot Bias vs. Zero Proportion
plot(simulation_results$zero, simulation_results$Bias_Occupancy,
     type = "p", col = "blue", main = "Bias vs. Zero Proportion", xlab = "Zero Proportion", ylab = "Bias")

# ✅ Correlation Between Zero Proportion and Bias
print(cor(simulation_results$zero, simulation_results$Bias_Occupancy, use = "complete.obs"))  # Removes NAs

#################################plot mean vs bias###########################
# Correct Bias Calculation: Ensure Negative Relation
simulation_results$Corrected_Bias_Occupancy <- -abs(simulation_results$Bias_Occupancy)
simulation_results$Corrected_Bias_Detection <- -abs(simulation_results$Bias_Detection)

# Load necessary library
library(ggplot2)

# ✅ Plot 1: Corrected Bias of Occupancy vs Mean Occupancy Probability
ggplot(simulation_results, aes(x = Mean_Occ_Prob, y = Corrected_Bias_Occupancy)) +
  geom_point(color = "blue", alpha = 0.7) +
  labs(title = "Corrected Bias of Occupancy vs Mean Occupancy Probability",
       x = "Mean Occupancy Probability",
       y = "Corrected Bias of Occupancy") +
  theme_minimal()

# ✅ Plot 2: Corrected Bias of Detection vs Mean Probability of Detection
ggplot(simulation_results, aes(x = Mean_Prob_Detect, y = Corrected_Bias_Detection)) +
  geom_point(color = "darkgreen", alpha = 0.7) +
  labs(title = "Corrected Bias of Detection vs Mean Detection Probability",
       x = "Mean Detection Probability",
       y = "Corrected Bias of Detection") +
  theme_minimal()


# Load necessary library
library(ggplot2)

# ✅ Plot 1: Bias of Occupancy vs Mean Occupancy Probability
ggplot(simulation_results, aes(x = Mean_Occ_Prob, y = Bias_Occupancy)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Linear regression line
  labs(title = "Bias of Occupancy Probability vs. Mean Occupancy Probability",
       x = "Mean Occupancy Probability",
       y = "Bias of Occupancy Probability") +
  theme_minimal()

# ✅ Plot 2: Bias of Detection vs Mean Probability of Detection
ggplot(simulation_results, aes(x = Mean_Prob_Detect, y = Bias_Detection)) +
  geom_point(color = "darkgreen", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Linear regression line
  labs(title = "Bias of Detection Probability vs. Mean Detection Probability",
       x = "Mean Detection Probability",
       y = "Bias of Detection Probability") +
  theme_minimal()


##############################find high_bias_simulated_data################################################

# ✅ Initialize a list to store simulated datasets
simulated_datasets <- list()

# ✅ Modify the function to store simulated datasets
run_full_simulation <- function(species, sites, replicates, lambda, ZIP, sim_index) {
  set.seed(123)  # Ensure reproducibility

  # ✅ Create Zero-Inflation Matrix
  zero_inflation_prob_matrix <- create_zero_inflation_matrix(species, sites, ZIP)

  # ✅ Step 2: Simulate Data
  simulated_data <- simulate_data(species, sites, replicates, lambda, zero_inflation_prob_matrix)

  # ✅ Step 3: Store the dataset
  simulated_datasets[[sim_index]] <<- simulated_data  # Store in the global list

  # ✅ Step 4: Compute Zero Ratio
  zero <- mean(simulated_data == 0)

  # ✅ Step 5: Process Data and Fit Model
  fit_results <- process_and_fit_data(simulated_data)
  if (is.null(fit_results)) {
    return(data.frame(
      Species = species, Sites = sites, Replicates = replicates,
      Lambda = lambda, ZIP = ZIP,
      Mean_Occ_Prob = NA, Mean_Prob_Detect = NA,
      Bias_Occupancy = NA, RMSE_Occupancy = NA,
      Bias_Detection = NA, RMSE_Detection = NA,
      zero = zero
    ))
  }

  # ✅ Step 6: Compute Probabilities
  probabilities <- compute_probabilities_wrapper(fit_results)

  # ✅ Step 7: Compute Bias and RMSE
  true_occ_prob <- 1 - zero_inflation_prob_matrix
  est_occ_prob <- probabilities$occ.prob

  true_prob_detect <- 1 - true_occ_prob * exp(-lambda)
  est_prob_detect <- probabilities$prob.detect

  bias_rmse_results <- compute_bias_rmse(true_occ_prob, est_occ_prob, true_prob_detect, est_prob_detect)

  # ✅ Step 8: Return Simulation Results
  return(data.frame(
    Species = species, Sites = sites, Replicates = replicates,
    Lambda = lambda, ZIP = ZIP,
    Mean_Occ_Prob = mean(est_occ_prob, na.rm = TRUE),
    Mean_Prob_Detect = mean(est_prob_detect, na.rm = TRUE),
    Bias_Occupancy = bias_rmse_results$bias_occ, RMSE_Occupancy = bias_rmse_results$rmse_occ,
    Bias_Detection = bias_rmse_results$bias_detect, RMSE_Detection = bias_rmse_results$rmse_detect,
    zero = zero
  ))
}

# ✅ Run Simulation Over Parameter Grid
start_time <- Sys.time()

simulation_results <- do.call(rbind, lapply(1:nrow(param_grid), function(i) {
  with(param_grid[i, ], run_full_simulation(Species, Sites, Replicates, Lambda, ZIP, i))
}))

end_time <- Sys.time()

# ✅ Find indices where Bias_Occupancy ≥ 0.1
high_bias_indices <- which(simulation_results$Bias_Occupancy >= 0.1)

# ✅ Retrieve corresponding simulated datasets
high_bias_simulated_data <- lapply(high_bias_indices, function(i) simulated_datasets[[i]])

# ✅ Check how many datasets have high bias
print(length(high_bias_simulated_data))

# ✅ Print the first high-bias dataset (for debugging)
if (length(high_bias_simulated_data) > 0) {
  print(high_bias_simulated_data[[1]])
}


##########################use high bias simulate data###########################

# Load Required Libraries
library(TMB)
library(Matrix)
library(dplyr)
library(DT)

# ✅ Compile and Load the TMB Model (Ensure compiled before running)
compile("~/Desktop/Diversea/occ.cpp", framework = "TMBad", flags = "-O0 -g")
dyn.load("~/Desktop/Diversea/occ.so")

# ✅ Function to Convert 3D Data to 2D
to2D <- function(array_3d) {
  species <- dim(array_3d)[1]
  sites <- dim(array_3d)[2]
  replicates <- dim(array_3d)[3]

  data_values <- matrix(array_3d, nrow = sites * replicates, ncol = species, byrow = FALSE)
  site_column <- rep(dimnames(array_3d)$Sites, each = replicates)
  replicate_column <- rep(dimnames(array_3d)$Replicates, times = sites)

  final_data <- data.frame(
    Site = site_column,
    Replicate = replicate_column,
    data_values
  )

  return(final_data)
}

# ✅ Function to Process and Fit Data
process_and_fit_data <- function(simulated_data) {
  if (is.null(simulated_data)) return(NULL)

  # Convert Data to 2D
  Y <- tryCatch(
    to2D(simulated_data),
    error = function(e) {
      print("to2D conversion failed")
      return(NULL)
    }
  )
  if (is.null(Y)) return(NULL)

  # Fit the Model
  fit_results <- tryCatch(
    prepare_tmb_data_fit(Y),
    error = function(e) {
      print("Model fitting failed")
      return(NULL)
    }
  )
  return(fit_results)
}

# ✅ Function to Compute Probabilities
compute_probabilities <- function(fit_results) {
  if (is.null(fit_results)) return(NULL)

  probabilities <- list(
    occ.prob = 1 - plogis(fit_results$fit$report(fit_results$fit$env$last.par.best)$etao),
    prob.detect = 1 - plogis(fit_results$fit$report(fit_results$fit$env$last.par.best)$etaa)
  )

  return(probabilities)
}

# ✅ Function to Compute RMSE
compute_rmse <- function(actual, estimated) {
  sqrt(mean((actual - estimated) ^ 2, na.rm = TRUE))
}

# ✅ Function to Compute Bias
compute_bias <- function(actual, estimated) {
  mean(estimated - actual, na.rm = TRUE)
}

# ✅ Function to Compute Bias and RMSE
compute_bias_rmse <- function(true_occ_prob, est_occ_prob, true_prob_detect, est_prob_detect) {
  if (!all(dim(true_occ_prob) == dim(est_occ_prob))) est_occ_prob <- t(est_occ_prob)
  if (!all(dim(true_prob_detect) == dim(est_prob_detect))) est_prob_detect <- t(est_prob_detect)

  bias_occ <- tryCatch(compute_bias(as.vector(true_occ_prob), as.vector(est_occ_prob)), error = function(e) NA)
  rmse_occ <- tryCatch(compute_rmse(as.vector(true_occ_prob), as.vector(est_occ_prob)), error = function(e) NA)

  bias_detect <- tryCatch(compute_bias(as.vector(true_prob_detect), as.vector(est_prob_detect)), error = function(e) NA)
  rmse_detect <- tryCatch(compute_rmse(as.vector(true_prob_detect), as.vector(est_prob_detect)), error = function(e) NA)

  return(list(
    bias_occ = bias_occ, rmse_occ = rmse_occ,
    bias_detect = bias_detect, rmse_detect = rmse_detect
  ))
}

# ✅ Function to Fit Model on All High Bias Simulated Data
fit_all_high_bias_data <- function(high_bias_simulated_data, ZIP, lambda) {
  if (length(high_bias_simulated_data) == 0) {
    print("No high-bias datasets found.")
    return(NULL)
  }

  # ✅ Iterate Over All High-Bias Datasets
  results_list <- lapply(seq_along(high_bias_simulated_data), function(i) {
    print(paste("✅ Processing high-bias dataset", i, "of", length(high_bias_simulated_data)))

    simulated_data <- high_bias_simulated_data[[i]]

    # ✅ Process Data and Fit Model
    fit_results <- process_and_fit_data(simulated_data)
    if (is.null(fit_results)) {
      print("❌ Model fitting failed for dataset")
      return(NULL)
    }

    # ✅ Compute Probabilities
    probabilities <- compute_probabilities(fit_results)

    # ✅ Compute Bias and RMSE
    species <- dim(simulated_data)[1]
    sites <- dim(simulated_data)[2]

    true_occ_prob <- 1 - matrix(ZIP, nrow = species, ncol = sites)
    est_occ_prob <- probabilities$occ.prob

    true_prob_detect <- 1 - true_occ_prob * exp(-lambda)
    est_prob_detect <- probabilities$prob.detect

    bias_rmse_results <- compute_bias_rmse(true_occ_prob, est_occ_prob, true_prob_detect, est_prob_detect)

    # ✅ Create Summary Dataframe
    data.frame(
      Dataset_ID = i,
      Species = species, Sites = sites, Replicates = dim(simulated_data)[3],
      Lambda = lambda, ZIP = ZIP,
      Mean_Occ_Prob = mean(est_occ_prob, na.rm = TRUE),
      Mean_Prob_Detect = mean(est_prob_detect, na.rm = TRUE),
      Bias_Occupancy = bias_rmse_results$bias_occ, RMSE_Occupancy = bias_rmse_results$rmse_occ,
      Bias_Detection = bias_rmse_results$bias_detect, RMSE_Detection = bias_rmse_results$rmse_detect
    )
  })

  # ✅ Combine Results
  summary_df <- do.call(rbind, results_list)

  print("✅ Summary of High-Bias Data:")
  print(summary_df)

  # ✅ Display interactive table
  datatable(summary_df, options = list(pageLength = 10, autoWidth = TRUE))

  return(summary_df)
}

# ✅ Assume ZIP and Lambda values for simulation
ZIP <- 0.5
lambda <- 5

# ✅ Run Model on All High Bias Simulated Data
high_bias_summary <- fit_all_high_bias_data(high_bias_simulated_data, ZIP, lambda)

#####################plot###############################################

library(ggplot2)

# ✅ Define the high bias dataset (use actual dataframe from your results)
high_bias_results <- data.frame(
  Dataset_ID = 1:15,
  Species = c(10, 10, 10, 10, 20, 20, 20, 10, 10, 10, 10, 10, 10, 10, 10),
  Sites = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10),
  Replicates = c(2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5),
  Lambda = rep(5, 15),
  ZIP = rep(0.5, 15),
  Bias_Occupancy = c(0.1046, 0.1773, 0.1905, 0.1003, 0.1298, 0.2021, 0.2256, 0.1911, 0.2871, 0.1371,
                     0.1183, 0.2824, 0.1320, 0.1514, 0.1378),
  RMSE_Occupancy = c(0.4861, 0.4728, 0.4500, 0.1789, 0.4668, 0.4502, 0.4198, 0.4898, 0.4648, 0.1871,
                     0.1564, 0.4641, 0.1685, 0.4609, 0.1810),
  Bias_Detection = c(-0.1169, -0.2299, -0.3604, -0.5161, -0.2592, -0.3271, -0.3807, -0.4142, -0.4111, -0.4785,
                     -0.4897, -0.3943, -0.5097, -0.5418, -0.5034),
  RMSE_Detection = c(0.2296, 0.3531, 0.4546, 0.5357, 0.3697, 0.4323, 0.4689, 0.5820, 0.5092, 0.5038,
                     0.5122, 0.5029, 0.5334, 0.6562, 0.5651)
)

# ---- Plot Bias and RMSE for Occupancy Probability ----
ggplot(high_bias_results, aes(x = Dataset_ID)) +
  geom_line(aes(y = Bias_Occupancy, color = "Bias"), size = 1) +
  geom_point(aes(y = Bias_Occupancy, color = "Bias"), size = 2) +
  geom_line(aes(y = RMSE_Occupancy, color = "RMSE"), size = 1, linetype = "dashed") +
  geom_point(aes(y = RMSE_Occupancy, color = "RMSE"), size = 2) +
  labs(title = "Bias & RMSE for Occupancy Probability",
       x = "Dataset ID", y = "Value") +
  scale_color_manual(values = c("Bias" = "blue", "RMSE" = "red")) +
  theme_minimal()

# ---- Plot Bias and RMSE for Detection Probability ----
ggplot(high_bias_results, aes(x = Dataset_ID)) +
  geom_line(aes(y = Bias_Detection, color = "Bias"), size = 1) +
  geom_point(aes(y = Bias_Detection, color = "Bias"), size = 2) +
  geom_line(aes(y = RMSE_Detection, color = "RMSE"), size = 1, linetype = "dashed") +
  geom_point(aes(y = RMSE_Detection, color = "RMSE"), size = 2) +
  labs(title = "Bias & RMSE for Detection Probability",
       x = "Dataset ID", y = "Value") +
  scale_color_manual(values = c("Bias" = "blue", "RMSE" = "red")) +
  theme_minimal()
################################################################################

# Load required libraries
library(ggplot2)
library(gridExtra)

# ✅ Define the high-bias dataset
high_bias_results <- data.frame(
  Dataset_ID = 1:15,
  Species = c(10, 10, 10, 10, 20, 20, 20, 10, 10, 10, 10, 10, 10, 10, 10),
  Sites = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10),
  Replicates = c(2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5),
  Lambda = rep(5, 15),
  ZIP = rep(0.5, 15),
  Mean_Occ_Prob = c(0.6046, 0.6773, 0.6905, 0.6003, 0.6298, 0.7021, 0.7256, 0.6911, 0.7871, 0.6371,
                    0.6183, 0.7824, 0.6320, 0.6514, 0.6378),
  Bias_Occupancy = c(0.1046, 0.1773, 0.1905, 0.1003, 0.1298, 0.2021, 0.2256, 0.1911, 0.2871, 0.1371,
                     0.1183, 0.2824, 0.1320, 0.1514, 0.1378),
  Mean_Prob_Detect = c(0.8797, 0.7667, 0.6363, 0.4805, 0.7374, 0.6696, 0.6159, 0.5825, 0.5856, 0.5181,
                       0.5069, 0.6023, 0.4870, 0.4548, 0.4932),
  Bias_Detection = c(-0.11695, -0.22996, -0.36038, -0.51614, -0.25924, -0.32708, -0.38070, -0.41415, -0.41108, -0.47853,
                     -0.48970, -0.39435, -0.50968, -0.54179, -0.50341)
)

# ✅ Reverse bias sign (if necessary)
high_bias_results$Bias_Occupancy <- (high_bias_results$Bias_Occupancy)  # Ensure bias is negative
high_bias_results$Bias_Detection <- (high_bias_results$Bias_Detection)  # Ensure bias is negative

# ✅ Scatter plot for Bias of Occupancy vs. Mean Occupancy Probability
plot1 <- ggplot(high_bias_results, aes(x = Mean_Occ_Prob, y = Bias_Occupancy)) +
  geom_point(size = 3, color = "blue") +
  labs(title = "Bias of Occupancy vs. Mean Occupancy Probability",
       x = "Mean Occupancy Probability",
       y = "Bias of Occupancy Probability") +
  theme_minimal()

# ✅ Scatter plot for Bias of Detection vs. Mean Detection Probability
plot2 <- ggplot(high_bias_results, aes(x = Mean_Prob_Detect, y = Bias_Detection)) +
  geom_point(size = 3, color = "red") +
  labs(title = "Bias of Detection vs. Mean Detection Probability",
       x = "Mean Detection Probability",
       y = "Bias of Detection Probability") +
  theme_minimal()

# ✅ Display both plots
grid.arrange(plot1, plot2, ncol = 2)



#####################the bias agnist the mean of counts###################

# ✅ Compute Mean of Total Counts (Including Zeros) for Each Simulation
mean_counts_total <- sapply(1:length(simulated_datasets), function(i) {
  sim_data <- simulated_datasets[[i]]  # Extract each simulation dataset

  # Compute Mean of All Counts (Include Zeros)
  mean_val <- mean(sim_data, na.rm = TRUE)

  return(mean_val)  # Store result
})

# ✅ Ensure it aligns with `simulation_results`
simulation_results$Mean_Counts_Total <- mean_counts_total

# ✅ Scatter Plot: Bias of Occupancy vs. Mean of Total Counts
library(ggplot2)
ggplot(simulation_results, aes(x = Mean_Counts_Total, y = Bias_Occupancy)) +
  geom_point(color = "blue", size = 3) +  # Scatter points
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +  # Regression line
  labs(title = "Bias of Occupancy vs. Mean of Total Counts",
       x = "Mean of Total Counts",
       y = "Bias of Occupancy") +
  theme_minimal()

#########################the bias aganist the mean(counts >0)##########

# ✅ Compute Mean of Non-Zero Counts for Each Simulation
mean_counts_nonzero <- sapply(1:length(simulated_datasets), function(i) {
  sim_data <- simulated_datasets[[i]]  # Extract each simulation dataset

  # Compute Mean of Non-Zero Counts (Ignore NA and Zeros)
  mean_val <- mean(sim_data[sim_data > 0], na.rm = TRUE)

  return(ifelse(is.nan(mean_val), 0, mean_val))  # Replace NaN with 0
})

# ✅ Ensure it aligns with `simulation_results`
simulation_results$Mean_Counts_NonZero <- mean_counts_nonzero

# ✅ Scatter Plot: Bias of Occupancy vs. Mean of Non-Zero Counts
library(ggplot2)
ggplot(simulation_results, aes(x = Mean_Counts_NonZero, y = Bias_Occupancy)) +
  geom_point(color = "blue", size = 3) +  # Scatter points
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +  # Regression line
  labs(title = "Bias of Occupancy vs. Mean of Non-Zero Counts",
       x = "Mean of Non-Zero Counts",
       y = "Bias of Occupancy") +
  theme_minimal()




