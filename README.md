[README_eDNAModel.md](https://github.com/user-attachments/files/34268662/README_eDNAModel.md)
# eDNAModel

**Hierarchical occupancy, capture, and abundance modelling for environmental DNA data**

`eDNAModel` is an R package for analysing environmental DNA (eDNA) metabarcoding data with hierarchical statistical models. It provides two complementary approaches:

- **`FitModel_joint()`**: joint estimation of occupancy, capture, and abundance processes using Template Model Builder (TMB).
- **`FitModel_gllvm()`**: hierarchical modelling with generalized linear latent variable models (GLLVMs) to explore residual multivariate structure in species occupancy.

eDNA sequencing counts are affected by ecological occurrence, imperfect capture or detection, and variation in the observation process. In particular, a zero count does not necessarily imply that a species is absent. The package helps distinguish these processes while accommodating differences among operational taxonomic units (OTUs).

## Installation

Install **`phyloseq` first** from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
```

Then install `eDNAModel` from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("marwah83/eDNAModel")
```

Load the packages:

```r
library(phyloseq)
library(eDNAModel)
```

**Note:** TMB-based models require an appropriate C++ compilation toolchain for your operating system.

## Statistical framework

The package distinguishes three processes:

1. **Occupancy:** whether OTU *k* is present at site *i*, represented by a latent state \(Z_{ik} \sim \mathrm{Bernoulli}(\psi_{ik})\).
2. **Capture:** whether an OTU present at a site is captured or detected in a replicate sample.
3. **Abundance:** variation in observed sequence counts, modelled using a supported count distribution.

These components address the difference between ecological absence and a failure to observe DNA from a species that is present.

## Example dataset: Marine invasive species in Trondheim

The included `physeq_new` object contains a dataset named `Marine Invasive species Trondheim`. In this dataset, the following structure has been verified:

| Characteristic | Value |
| --- | ---: |
| Unique sites (`Name`) | 75 |
| Samples | 150 |
| Replicates per site | 2 |
| Unique `Name`–`Replicate` combinations | 150 |
| Distinct sampling months within each site | 1 |

Load the dataset and create a unique sample identifier:

```r
library(phyloseq)
library(eDNAModel)

data("physeq_new", package = "eDNAModel")
physeq_one <- physeq_new$`Marine Invasive species Trondheim`

# Preserve the original unique phyloseq sample names as metadata.
phyloseq::sample_data(physeq_one)$SampleID <-
  phyloseq::sample_names(physeq_one)

# Inspect the sample metadata.
sample_df <- as(phyloseq::sample_data(physeq_one), "data.frame")

stopifnot(
  nrow(sample_df) == 150L,
  length(unique(sample_df$Name)) == 75L,
  length(unique(sample_df$SampleID)) == 150L,
  all(table(sample_df$Name) == 2L),
  !anyDuplicated(sample_df[c("Name", "Replicate")])
)
```

The examples below use the same sampling hierarchy:

| Function argument | Column | Interpretation |
| --- | --- | --- |
| `site_col` | `"Name"` | Site-level occupancy unit |
| `sample_col` | `"SampleID"` | Unique sample identifier |
| `replicate_col` | `"Replicate"` | Replicate within a site |

The two replicates must represent repeated observations of the same ecological site for the occupancy–detection interpretation to apply.

## 1. `FitModel_joint()`: Joint hierarchical model using TMB

`FitModel_joint()` jointly estimates occupancy, capture, and abundance components with TMB. It supports OTU-specific random effects and several count distributions, including Poisson, negative binomial, zero-inflated Poisson, and zero-inflated negative binomial.

### Fit a model

The following example includes `Samplingmonth` as a site-level occupancy predictor. Because each site has one sampling month, the covariate is consistent across its two replicates.

```r
set.seed(123)

out_joint <- FitModel_joint(
  phyloseq = physeq_one,
  site_col = "Name",
  sample_col = "SampleID",
  replicate_col = "Replicate",
  occupancy_formula = ~ Samplingmonth,
  capture_formula = ~ 1,
  abundance_formula = ~ 1,
  abundance_family = "poisson",
  min_species_sum = 10,
  min_detection_replicates = 1,
  random_occ_otu = TRUE,
  random_capture_otu = TRUE,
  random_abund_otu = TRUE,
  random_sample = TRUE,
  random_sample_otu = FALSE,
  gradient_tol = 1e-3,
  max_restarts = 2,
  verbose = TRUE
)
```

To fit an intercept-only baseline occupancy model, change `occupancy_formula = ~ Samplingmonth` to `occupancy_formula = ~ 1`.

### Inspect results

```r
names(out_joint)
str(out_joint, max.level = 1)
```

### Other abundance families

Change `abundance_family` to one of:

```r
"poisson"  # Poisson
"nbinom"  # Negative binomial
"zip"     # Zero-inflated Poisson
"zinb"    # Zero-inflated negative binomial
```

The appropriate family depends on the observed count distribution and model diagnostics. More complex families may require additional computation or stronger data support.

### Selected optimization controls

| Argument | Purpose |
| --- | --- |
| `gradient_tol` | Gradient convergence tolerance |
| `max_restarts` | Maximum optimization restarts |
| `iter_max` | Maximum optimizer iterations |
| `eval_max` | Maximum objective-function evaluations |
| `otu_intervals` | Request OTU-level uncertainty intervals |
| `otu_conf_level` | Confidence level for OTU intervals |
| `get_report_covariance` | Request covariance for reported quantities |
| `get_joint_precision` | Request joint precision information |

See `?FitModel_joint` for all supported arguments and returned components.

## 2. `FitModel_gllvm()`: Hierarchical model with latent variables

`FitModel_gllvm()` incorporates a GLLVM into the occupancy component to represent residual multivariate variation among OTUs. In a general GLLVM occupancy formulation,

### Generalized Linear Latent Variable Model (GLLVM)

The `eDNAModel` package incorporates a Generalized Linear Latent Variable Model (GLLVM) into the occupancy component to represent residual multivariate variation among OTUs.

In a general GLLVM occupancy formulation,

$$
\mathrm{logit}(\psi_{ik}) =
\mathbf{x}_i^{\top}\boldsymbol{\beta}_k +
\mathbf{u}_i^{\top}\boldsymbol{\lambda}_k
$$

where:

- $\psi_{ik}$ is the occupancy probability of OTU $k$ at site $i$.
- $\mathbf{x}_i$ is the vector of measured site-level covariates.
- $\boldsymbol{\beta}_k$ is the vector of OTU-specific regression coefficients.
- $\mathbf{u}_i$ is the vector of latent site scores.
- $\boldsymbol{\lambda}_k$ is the vector of OTU-specific latent factor loadings.

The latent variables capture residual associations among OTUs that are not explained by the measured environmental predictors. This allows the GLLVM to identify underlying community patterns and dependencies among OTUs beyond those explained by the observed environmental covariates.


### Fit a model

```r
set.seed(123)

out_gllvm <- FitModel_gllvm(
  phyloseq = physeq_one,
  site_col = "Name",
  sample_col = "SampleID",
  replicate_col = "Replicate",
  abundance_rhs = (1 | Samplingmonth / OTU) + (1 | OTU),
  occupancy_covars = "Samplingmonth",
  abundance_family = "poisson",
  min_species_sum = 10,
  min_detection_replicates = 1,
  n_iter = 50,
  burn_in = 10,
  num_lv_c = 0,
  num_lv = 1,
  gllvm_n_init = 10,
  gllvm_gradient_tol = 1e-3,
  gllvm_strict = FALSE,
  gllvm_condition_tol = 1e8,
  gllvm_eigen_tol = 1e-8,
  gllvm_warn_separation = TRUE,
  gllvm_auto_exclude = TRUE,
  gllvm_min_occupied_sites = 2,
  gllvm_excluded_prior = c(1, 1),
  verbose = TRUE
)
```

**Version note:** This example uses the updated `FitModel_gllvm()` interface discussed in the development version. Check `args(eDNAModel::FitModel_gllvm)` if your installed version reports an unused argument. The `abundance_rhs` syntax and required metadata should also be checked against the current function documentation.

### Inspect results and diagnostics

```r
names(out_gllvm)
head(out_gllvm$summary)

# Numerical diagnostics and screened taxa
out_gllvm$gllvm_diagnostics
out_gllvm$gllvm_excluded_otus
out_gllvm$gllvm_separation

# Latent-variable summaries
head(out_gllvm$mean_lv_sites)
head(out_gllvm$mean_lv_species)
```

### Important GLLVM controls

| Argument | Purpose |
| --- | --- |
| `num_lv` | Number of unconstrained latent variables |
| `num_lv_c` | Number of constrained latent variables |
| `gllvm_n_init` | Number of initialization attempts per GLLVM fit |
| `gllvm_gradient_tol` | Gradient convergence tolerance |
| `gllvm_strict` | Controls strict treatment of numerical diagnostics |
| `gllvm_condition_tol` | Hessian conditioning threshold |
| `gllvm_eigen_tol` | Hessian eigenvalue tolerance |
| `gllvm_warn_separation` | Enable separation warnings |
| `gllvm_auto_exclude` | Screen problematic OTUs in the GLLVM occupancy step |
| `gllvm_min_occupied_sites` | Occupancy-variation screening threshold |
| `gllvm_excluded_prior` | Beta-prior smoothing parameters for screened OTUs |

For example, `num_lv = 1` fits one unconstrained latent dimension. Increasing to `num_lv = 2` may capture more complex community structure but can make estimation more difficult. Always inspect convergence diagnostics before interpreting the latent scores or loadings.

### Probability transformations

The interpretation of a stored linear predictor depends on the link used by that component:

```r
# Inverse logit for logit-linked occupancy/capture predictors
psi <- plogis(eta_occupancy)
capture <- plogis(eta_capture)

# Inverse complementary log-log for cloglog-linked detection predictors
p_detect <- -expm1(-exp(eta_detection))
```

These are illustrations of the inverse-link functions, **not** direct commands on `out_gllvm`: replace the `eta_*` placeholders with the appropriate numeric linear predictors from the returned model objects. Do not apply a link inversion twice if a returned value is already a probability.

## Which function should I use?

| Feature | `FitModel_joint()` | `FitModel_gllvm()` |
| --- | --- | --- |
| Hierarchical occupancy component | Yes | Yes |
| Capture/detection component | Yes | Yes |
| Abundance component | Yes | Yes |
| TMB joint likelihood | Yes | Different fitting framework |
| GLLVM latent occupancy structure | No | Yes |
| OTU-specific effects | Yes | Supported through model specification |
| GLLVM separation screening | Not applicable | Yes |
| Primary use | Joint estimation of ecological and observation processes | Multivariate occupancy and residual community structure |

Choose **`FitModel_joint()`** when your primary goal is joint estimation of occupancy, capture, and abundance. Choose **`FitModel_gllvm()`** when you also want to explore latent multivariate community structure.

## Practical modelling guidance

- **Check replication:** confirm that the site identifier groups genuine replicate observations.
- **Check covariates:** occupancy covariates should be defined at the intended occupancy-unit level.
- **Inspect OTU prevalence:** very rare or nearly ubiquitous OTUs can cause separation or weak identification.
- **Check numerical diagnostics:** optimizer success alone does not guarantee a well-identified model.
- **Compare specifications carefully:** when comparing functions, align covariates, filtering rules, and count families where supported.
- **Start simple:** begin with an intercept-only or one-latent-variable model before adding complexity.

## Documentation

```r
?FitModel_joint
?FitModel_gllvm

vignette(package = "eDNAModel")
help(package = "eDNAModel")
```

## Development and issues

Source code, updates, and issue reporting:

**https://github.com/marwah83/eDNAModel**

## Citation

If you use `eDNAModel` in your research, please cite the package and relevant methodological publications. If citation metadata are included in the installed package, you can retrieve them with:

```r
citation("eDNAModel")
```

## License

Refer to the repository's `DESCRIPTION` and license files for the current licensing terms.
