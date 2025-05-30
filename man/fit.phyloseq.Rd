\name{fit.phyloseq}
\alias{fit.phyloseq}
\title{Fit eDNA Occupancy-Abundance Model from a Phyloseq Object}
\description{
Fits a hierarchical multispecies occupancy-abundance model using the TMB engine, starting directly from a \code{phyloseq} object.
}
\usage{
fit.phyloseq(
  phyloseq_obj,
  a.formula = ~1,
  o.formula = ~1,
  linko = 1,
  linka = 0,
  family = 1,
  Ntrials = matrix(0),
  offset = NULL,
  control = list(maxit = 10000, trace = 1),
  verbose = TRUE,
  ...
)
}
\arguments{
  \item{phyloseq_obj}{A \code{phyloseq} object containing an OTU table and sample metadata.}
  \item{a.formula}{Formula for the abundance component. Default is \code{~1}.}
  \item{o.formula}{Formula for the occupancy component. Default is \code{~1}.}
  \item{linko}{Link function for occupancy. \code{1 = logit} (default), \code{2 = probit}.}
  \item{linka}{Link function for abundance. \code{0 = log} (default), \code{1 = logit}, \code{2 = probit}, \code{3 = cloglog}.}
  \item{family}{Model family. \code{0 = ZIP}, \code{1 = ZINB} (default), \code{2 = Binomial}.}
  \item{Ntrials}{Matrix of trials (only for binomial family). Default is a matrix of 0s.}
  \item{offset}{Optional offset matrix, same dimension as \code{y}. Default is \code{NULL}.}
  \item{control}{A list of control arguments passed to the optimizer (e.g. \code{maxit}, \code{trace}).}
  \item{verbose}{Logical. If \code{TRUE}, prints progress messages during preprocessing and fitting.}
  \item{...}{Additional arguments passed to \code{\link{run_full_TMB}}.}
}
\details{
This function handles preprocessing of the \code{phyloseq} object and extracts necessary inputs for hierarchical modeling. It performs:
\itemize{
  \item Transposing the OTU table if taxa are rows,
  \item Filtering for biological samples based on \code{sampletype},
  \item Inferring \code{Site} from \code{location}, and parsing \code{Sample} and \code{Replicate} from sample names.
}
The function returns a fitted model using \code{\link{run_full_TMB}}.
}
\value{
An object of class \code{eDNAModel} containing:
\itemize{
  \item Model fit object from TMB (\code{$TMBobj}),
  \item Estimated occupancy and abundance probabilities,
  \item Optimization result,
  \item Input \code{y} and \code{X} matrices used for modeling.
}
}
\seealso{
\code{\link{run_full_TMB}}, \code{\link{compute_residuals_TMB}}, \code{\link{summary.eDNAModel}}, \code{\link{phyloseq}}
}
\examples{
# Load example phyloseq object
physeq <- readRDS(system.file("extdata", "longdataexample.RDS", package = "eDNAModel"))

# Fit model with default settings
model <- fit.phyloseq(physeq)

# Summarize results
summary(model)
}
