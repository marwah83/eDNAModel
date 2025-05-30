\name{fit.phyloseq}
\alias{fit.phyloseq}
\title{Fit eDNA Occupancy-Abundance Model from a Phyloseq Object}
\description{
Fits a hierarchical multispecies occupancy-abundance model using the TMB engine, starting directly from a \code{phyloseq} object.
}
\usage{
fit.phyloseq(phyloseq_obj,
             a.formula = ~1,
             o.formula = ~1,
             linko = 1,
             linka = 0,
             family = 1,
             Ntrials = matrix(0),
             offset = NULL,
             control = list(maxit = 10000, trace = 1))
}
\arguments{
  \item{phyloseq_obj}{A \code{phyloseq} object containing an OTU table and sample metadata.}
  \item{a.formula}{Formula for the abundance component. Default is \code{~1}.}
  \item{o.formula}{Formula for the occupancy component. Default is \code{~1}.}
  \item{linko}{Link function for occupancy. 1 = logit (default), 2 = probit.}
  \item{linka}{Link function for abundance. 0 = log (default), 1 = logit, 2 = probit, 3 = cloglog.}
  \item{family}{Model family. 0 = ZIP, 1 = ZINB (default), 2 = Binomial.}
  \item{Ntrials}{Matrix of trials (only for binomial family). Default is a matrix of 0s.}
  \item{offset}{Optional offset matrix, same dimension as \code{y}. Default is \code{NULL}.}
  \item{control}{A list of control arguments passed to the optimizer (e.g. \code{maxit}, \code{trace}).}
}
\details{
This function handles the preprocessing of \code{phyloseq} objects, including:
\itemize{
  \item Transposing the OTU table (if taxa are rows)
  \item Filtering for biological samples (based on \code{sampletype})
  \item Inferring \code{Site}, \code{Sample}, and \code{Replicate} IDs from metadata and rownames
}
It then prepares the \code{y} and \code{X} matrices and calls \code{\link{run_full_TMB}} to fit the model.
}
\value{
An object of class \code{eDNAModel} containing:
\itemize{
  \item Model fit object from TMB
  \item Estimated occupancy and abundance probabilities
  \item Optimization results
  \item Input \code{y} and \code{X} used for fitting
}
}
\seealso{
\code{\link{run_full_TMB}}, \code{\link{phyloseq}}, \code{\link{summary.eDNAModel}}, \code{\link{compute_residuals_TMB}}
}
\examples{
# Load example phyloseq object
physeq <- readRDS(system.file("extdata", "longdataexample.RDS", package = "eDNAModel"))

# Fit model with default settings
model <- fit.phyloseq(physeq)

# Summarize results
summary(model)
}
