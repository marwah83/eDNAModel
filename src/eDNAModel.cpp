#include <TMB.hpp>

// ============================================================
// Numerically stable logistic helper functions
// ============================================================
//
// TMB already provides logspace_add(a,b):
//
//     logspace_add(a,b)
//     = log(exp(a) + exp(b))
//
// This is preferable to explicitly exponentiating probabilities
// because likelihood terms can become extremely small.
//
// We use these helper functions repeatedly for occupancy,
// capture, and zero-inflation probabilities.
// ============================================================

template<class Type>
Type log_invlogit(Type x)
{
// ----------------------------------------------------------
// Computes:
//
//     log(invlogit(x))
//
// where:
//
//     invlogit(x)
//     =
//     1 / (1 + exp(-x))
//
// Therefore:
//
//     log(invlogit(x))
//     =
//     -log(1 + exp(-x))
//
// Using logspace_add makes this numerically stable.
// ----------------------------------------------------------

return -logspace_add(
Type(0),
-x
);
}



template<class Type>
Type log1m_invlogit(Type x)
{
// ----------------------------------------------------------
// Computes:
//
//     log(1 - invlogit(x))
//
// Because:
//
//     1 - invlogit(x)
//     =
//     1 / (1 + exp(x))
//
// therefore:
//
//     log(1-p)
//     =
//     -log(1 + exp(x)).
// ----------------------------------------------------------

return -logspace_add(
Type(0),
x
);
}



// ============================================================
// Joint observed-data eDNA likelihood
// ============================================================
//
// Hierarchy:
//
// ------------------------------------------------------------
// OCCUPANCY LEVEL
//
//     Z_im ~ Bernoulli(psi_im)
//
// i = site
// m = OTU
//
// There is ONE occupancy state per site x OTU.
//
// ------------------------------------------------------------
// CAPTURE LEVEL
//
//     A_ijm | Z_im = 1
//       ~ Bernoulli(p_ijm)
//
// j = biological sample
//
// ------------------------------------------------------------
// PCR / READ LEVEL
//
//     Y_ijkm | A_ijm = 1
//       ~ Count(lambda_ijkm)
//
// k = PCR replicate
//
// The count distribution may be:
//
//     Poisson
//     NB2
//     ZIP
//     ZINB
//
// ------------------------------------------------------------
//
// IMPORTANT:
//
// Z and A are NOT sampled or optimized directly.
//
// Instead:
//
//     A is analytically marginalized
//
// and then:
//
//     Z is analytically marginalized.
//
// Gaussian random effects are integrated by TMB using the
// Laplace approximation when passed through the random
// argument of MakeADFun().
// ============================================================

template<class Type>
Type objective_function<Type>::operator() ()
{

// ==========================================================
// 1. DATA
// ==========================================================

DATA_VECTOR(y);

// Observed sequencing counts.
//
// Each y(r) is one PCR/read observation.

// ----------------------------------------------------------
// Design matrices
// ----------------------------------------------------------

DATA_MATRIX(X_occ);

// Occupancy design matrix.
//
// One row per:
//
//     SITE x OTU
//
// NOT one row per PCR observation.

DATA_MATRIX(X_cap);

// Capture design matrix.
//
// One row per:
//
//     BIOLOGICAL SAMPLE x OTU

DATA_MATRIX(X_abund);

// Abundance design matrix.
//
// One row per PCR/read observation.

// ----------------------------------------------------------
// Abundance offset
// ----------------------------------------------------------

DATA_VECTOR(offset_abund);

// For example:
//
//     log(total_reads)
//
// representing sequencing depth.

// ----------------------------------------------------------
// Hierarchical mapping indices
// ----------------------------------------------------------

DATA_IVECTOR(occ_otu);

// For each SITE x OTU occupancy group g,
// gives the OTU index m.
//
// Used by:
//
//     b_occ_otu(m)

DATA_IVECTOR(sample_otu);

// For every SAMPLE x OTU group s,
// gives the OTU index m.
//
// Used by:
//
//     b_cap_otu(m)

DATA_IVECTOR(row_sample_group);

// Maps:
//
//     PCR/read row r
//
// to:
//
//     SAMPLE x OTU group s.

DATA_IVECTOR(sample_site_group);

// Maps:
//
//     SAMPLE x OTU group s
//
// to:
//
//     SITE x OTU group g.
//
// This ensures multiple samples from a site share
// the same occupancy state.

DATA_IVECTOR(row_otu);

// Maps each PCR/read observation r to its OTU m.
//
// Used for:
//
//     abundance OTU random effects
//
// and now also:
//
//     zero-inflation OTU random effects.

DATA_IVECTOR(row_sample_id);

// Maps every PCR/read observation to its biological sample.

DATA_IVECTOR(row_sample_otu_re);

// Maps every PCR/read observation to a unique
// sample x OTU random-effect level.

// ----------------------------------------------------------
// Positive-history indicators
// ----------------------------------------------------------

DATA_IVECTOR(sample_positive);

// sample_positive(s) = 1 if at least one PCR replicate
// for sample x OTU group s is positive.
//
// If positive:
//
//     A_ijm MUST equal 1.

DATA_IVECTOR(site_positive);

// site_positive(g) = 1 if at least one observation anywhere
// inside site x OTU group g is positive.
//
// If positive:
//
//     Z_im MUST equal 1.

// ----------------------------------------------------------
// Number of hierarchical groups
// ----------------------------------------------------------

DATA_INTEGER(n_site_groups);

// Number of SITE x OTU combinations.

DATA_INTEGER(n_sample_groups);

// Number of SAMPLE x OTU combinations.

// ----------------------------------------------------------
// Count family
// ----------------------------------------------------------
//
// 0 = Poisson
// 1 = NB2
// 2 = ZIP
// 3 = ZINB
// ----------------------------------------------------------

DATA_INTEGER(family_code);

// ==========================================================
// 2. RANDOM-EFFECT SWITCHES
// ==========================================================

DATA_INTEGER(use_occ_otu);

// OTU random effect for occupancy.

DATA_INTEGER(use_cap_otu);

// OTU random effect for capture.

DATA_INTEGER(use_abund_otu);

// OTU random effect for abundance.

DATA_INTEGER(use_sample_re);

// Biological-sample abundance random effect.

DATA_INTEGER(use_sample_otu_re);

// Sample x OTU abundance random effect.

DATA_INTEGER(use_zi_otu);

// NEW:
//
// OTU-specific random effect for ZIP/ZINB structural zeros.
//
// 1 = estimate taxon-specific PCR dropout heterogeneity.
// 0 = all OTUs share one zero-inflation probability.



// ==========================================================
// 3. FIXED-EFFECT PARAMETERS
// ==========================================================

PARAMETER_VECTOR(beta_occ);

// Occupancy fixed effects.

PARAMETER_VECTOR(beta_cap);

// Capture fixed effects.

PARAMETER_VECTOR(beta_abund);

// Abundance fixed effects.



// ==========================================================
// 4. RANDOM EFFECTS
// ==========================================================

PARAMETER_VECTOR(b_occ_otu);

// OTU-specific occupancy deviations.

PARAMETER_VECTOR(b_cap_otu);

// OTU-specific capture deviations.

PARAMETER_VECTOR(b_abund_otu);

// OTU-specific abundance deviations.

PARAMETER_VECTOR(b_sample);

// Biological-sample abundance random effects.

PARAMETER_VECTOR(b_sample_otu);

// Sample x OTU abundance random effects.

PARAMETER_VECTOR(b_zi_otu);

// NEW:
//
// OTU-specific deviations in structural-zero probability.
//
// Mathematically:
//
//     b_zi_otu(m)
//       ~ Normal(0, sigma_zi^2)
//
// and:
//
//     logit(pi_m)
//       =
//       zi_intercept
//       +
//       b_zi_otu(m).



// ==========================================================
// 5. RANDOM-EFFECT STANDARD DEVIATIONS
// ==========================================================

PARAMETER(log_sd_occ_otu);

PARAMETER(log_sd_cap_otu);

PARAMETER(log_sd_abund_otu);

PARAMETER(log_sd_sample);

PARAMETER(log_sd_sample_otu);

PARAMETER(log_sd_zi_otu);

// NEW:
//
// log standard deviation of OTU-specific
// zero-inflation random effects.

Type sd_occ_otu =
exp(
log_sd_occ_otu
);

Type sd_cap_otu =
exp(
log_sd_cap_otu
);

Type sd_abund_otu =
exp(
log_sd_abund_otu
);

Type sd_sample =
exp(
log_sd_sample
);

Type sd_sample_otu =
exp(
log_sd_sample_otu
);

Type sd_zi_otu =
exp(
log_sd_zi_otu
);

// exp() guarantees positive standard deviations.



// ==========================================================
// 6. ABUNDANCE-DISTRIBUTION PARAMETERS
// ==========================================================

PARAMETER(log_theta);

// Shared NB2 dispersion parameter.

Type theta =
exp(
log_theta
);

// One theta is currently shared across all OTUs:
//
//     Var(Y)
//       =
//       mu + mu^2 / theta.

PARAMETER(zi_intercept);

// Global zero-inflation intercept.
//
// This now represents the population-average PCR dropout
// tendency across OTUs.
//
// OTUs can deviate from this through b_zi_otu.



// ==========================================================
// 7. NEGATIVE LOG-LIKELIHOOD
// ==========================================================

Type nll =
Type(0);

// TMB minimizes nll:
//
//     nll = -log L.



// ==========================================================
// 8. GAUSSIAN RANDOM-EFFECT DENSITIES
// ==========================================================

if (use_occ_otu == 1)
{

for (
  int m = 0;
  m < b_occ_otu.size();
  m++
)
{

  nll -= dnorm(
    b_occ_otu(m),
    Type(0),
    sd_occ_otu,
    true
  );

  // b_occ_otu(m)
  // ~ N(0, sigma_occ^2).
}

}



if (use_cap_otu == 1)
{

for (
  int m = 0;
  m < b_cap_otu.size();
  m++
)
{

  nll -= dnorm(
    b_cap_otu(m),
    Type(0),
    sd_cap_otu,
    true
  );
}

}



if (use_abund_otu == 1)
{

for (
  int m = 0;
  m < b_abund_otu.size();
  m++
)
{

  nll -= dnorm(
    b_abund_otu(m),
    Type(0),
    sd_abund_otu,
    true
  );
}

}



if (use_sample_re == 1)
{

for (
  int j = 0;
  j < b_sample.size();
  j++
)
{

  nll -= dnorm(
    b_sample(j),
    Type(0),
    sd_sample,
    true
  );
}

}



if (use_sample_otu_re == 1)
{

for (
  int q = 0;
  q < b_sample_otu.size();
  q++
)
{

  nll -= dnorm(
    b_sample_otu(q),
    Type(0),
    sd_sample_otu,
    true
  );
}

}



// ----------------------------------------------------------
// NEW: zero-inflation OTU random effects
// ----------------------------------------------------------

if (
use_zi_otu == 1 &&
(
family_code == 2 ||
family_code == 3
)
)
{

for (
  int m = 0;
  m < b_zi_otu.size();
  m++
)
{

  nll -= dnorm(
    b_zi_otu(m),
    Type(0),
    sd_zi_otu,
    true
  );

  // Assumption:
  //
  //     b_zi_m
  //     ~ Normal(
  //         0,
  //         sigma_zi^2
  //       )
  //
  // This allows each OTU to have a different
  // PCR-level dropout probability.
}

}



// ==========================================================
// 9. OCCUPANCY LINEAR PREDICTOR
// ==========================================================

vector<Type> eta_occ =
X_occ *
beta_occ;

// Fixed-effect component:
//
//     eta_occ_im
//       =
//       X_occ_im beta_occ.

if (use_occ_otu == 1)
{

for (
  int g = 0;
  g < n_site_groups;
  g++
)
{

  eta_occ(g) +=
    b_occ_otu(
      occ_otu(g)
    );

  // Full occupancy model:
  //
  // logit(psi_im)
  // =
  // X_occ_im beta_occ
  // +
  // b_occ_m.
}

}



// ==========================================================
// 10. CAPTURE LINEAR PREDICTOR
// ==========================================================

vector<Type> eta_cap =
X_cap *
beta_cap;

if (use_cap_otu == 1)
{

for (
  int s = 0;
  s < n_sample_groups;
  s++
)
{

  eta_cap(s) +=
    b_cap_otu(
      sample_otu(s)
    );

  // logit(p_ijm)
  // =
  // X_cap beta_cap
  // +
  // b_cap_m.
}

}



// ==========================================================
// 11. ABUNDANCE LINEAR PREDICTOR
// ==========================================================

vector<Type> eta_abund =
X_abund *
beta_abund;

for (
int r = 0;
r < y.size();
r++
)
{

// --------------------------------------------------------
// Offset
// --------------------------------------------------------

eta_abund(r) +=
  offset_abund(r);


// --------------------------------------------------------
// OTU abundance random effect
// --------------------------------------------------------

if (use_abund_otu == 1)
{

  eta_abund(r) +=
    b_abund_otu(
      row_otu(r)
    );
}


// --------------------------------------------------------
// Sample abundance random effect
// --------------------------------------------------------

if (use_sample_re == 1)
{

  eta_abund(r) +=
    b_sample(
      row_sample_id(r)
    );
}


// --------------------------------------------------------
// Sample x OTU abundance random effect
// --------------------------------------------------------

if (use_sample_otu_re == 1)
{

  eta_abund(r) +=
    b_sample_otu(
      row_sample_otu_re(r)
    );
}

}



// ==========================================================
// 12. ZERO-INFLATION LINEAR PREDICTOR
// ==========================================================
//
// NEW MODEL:
//
//     logit(pi_m)
//       =
//       zi_intercept
//       +
//       b_zi_otu(m)
//
// where:
//
//     b_zi_otu(m)
//       ~ N(0, sigma_zi^2).
//
// We construct one pi per OTU.
// ==========================================================

vector<Type> pi_zi_otu(
b_zi_otu.size()
);

for (
int m = 0;
m < b_zi_otu.size();
m++
)
{

Type eta_zi =
  zi_intercept;


if (
  use_zi_otu == 1 &&
  (
    family_code == 2 ||
    family_code == 3
  )
)
{

  eta_zi +=
    b_zi_otu(m);
}


pi_zi_otu(m) =
  Type(1) /
  (
    Type(1) +
    exp(
      -eta_zi
    )
  );

// pi_zi_otu(m)
//
// is the OTU-specific structural-zero probability.

}



// ==========================================================
// 13. READ-LEVEL COUNT LOG LIKELIHOOD
// ==========================================================

vector<Type> log_count(
y.size()
);

for (
int r = 0;
r < y.size();
r++
)
{

// --------------------------------------------------------
// Mean abundance
// --------------------------------------------------------

Type mu =
  exp(
    eta_abund(r)
  );

// Because abundance uses a log link:
//
//     lambda_ijkm
//       =
//       exp(eta_abund_ijkm).


Type lp =
  Type(0);


// --------------------------------------------------------
// Get OTU-specific zero-inflation probability
// --------------------------------------------------------

Type pi_zi =
  pi_zi_otu(
    row_otu(r)
  );

// Each observation uses the dropout probability
// corresponding to its OTU.



// ========================================================
// POISSON
// ========================================================

if (
  family_code == 0
)
{

  lp =
    dpois(
      y(r),
      mu,
      true
    );

  // log P(Y=y | mu).
}



// ========================================================
// NEGATIVE BINOMIAL 2
// ========================================================

else if (
  family_code == 1
)
{

  Type variance =
    mu +
    mu *
    mu /
    theta;


  lp =
    dnbinom2(
      y(r),
      mu,
      variance,
      true
    );

  // NB2:
  //
  // E(Y) = mu
  //
  // Var(Y)
  // =
  // mu + mu^2/theta.
}



// ========================================================
// ZERO-INFLATED POISSON
// ========================================================

else if (
  family_code == 2
)
{

  Type log_count_component =
    dpois(
      y(r),
      mu,
      true
    );


  // ------------------------------------------------------
  // Observed zero
  // ------------------------------------------------------

  if (
    y(r) ==
    Type(0)
  )
  {

    lp =
      logspace_add(

        log(
          pi_zi
        ),

        log(
          Type(1) -
          pi_zi
        ) +
        log_count_component

      );

    // For y = 0:
    //
    // P(Y=0)
    //
    // =
    //
    // pi_m
    //
    // +
    //
    // (1-pi_m)
    // P_Poisson(Y=0 | mu).
    //
    // pi_m is now OTU-specific.
  }


  // ------------------------------------------------------
  // Positive observation
  // ------------------------------------------------------

  else
  {

    lp =
      log(
        Type(1) -
        pi_zi
      ) +
      log_count_component;

    // A positive count cannot originate from the
    // structural-zero component.
  }
}



// ========================================================
// ZERO-INFLATED NEGATIVE BINOMIAL 2
// ========================================================

else if (
  family_code == 3
)
{

  Type variance =
    mu +
    mu *
    mu /
    theta;


  Type log_count_component =
    dnbinom2(
      y(r),
      mu,
      variance,
      true
    );


  if (
    y(r) ==
    Type(0)
  )
  {

    lp =
      logspace_add(

        log(
          pi_zi
        ),

        log(
          Type(1) -
          pi_zi
        ) +
        log_count_component

      );

    // ZINB zero:
    //
    // P(Y=0)
    //
    // =
    //
    // pi_m
    //
    // +
    //
    // (1-pi_m)
    // P_NB(Y=0 | mu, theta).
  }


  else
  {

    lp =
      log(
        Type(1) -
        pi_zi
      ) +
      log_count_component;
  }
}


// --------------------------------------------------------
// Store observation-level count log likelihood
// --------------------------------------------------------

log_count(r) =
  lp;

}



// ==========================================================
// 14. SAMPLE x OTU CONDITIONAL COUNT LIKELIHOOD
// ==========================================================
//
// First combine PCR replicates conditional on A = 1:
//
//     L_count_ijm
//
//     =
//
//     product_k
//     f(y_ijkm | A_ijm = 1)
//
// On log scale:
//
//     log L_count_ijm
//
//     =
//
//     sum_k
//     log f(y_ijkm).
// ==========================================================

vector<Type> log_sample_count(
n_sample_groups
);

log_sample_count.setZero();

for (
int r = 0;
r < y.size();
r++
)
{

int s =
  row_sample_group(r);


log_sample_count(s) +=
  log_count(r);

}



// ==========================================================
// 15. CAPTURE MARGINALIZATION
// ==========================================================

vector<Type> log_sample(
n_sample_groups
);

for (
int s = 0;
s < n_sample_groups;
s++
)
{

Type log_p =
  log_invlogit(
    eta_cap(s)
  );


Type log_1mp =
  log1m_invlogit(
    eta_cap(s)
  );


// --------------------------------------------------------
// Positive biological sample
// --------------------------------------------------------
//
// If at least one PCR observation is positive:
//
//     A_ijm = 1
//
// must hold.
//
// Therefore:
//
// L_ijm
//
// =
//
// p_ijm
// *
// product_k f(y_ijkm).
// --------------------------------------------------------


if (
  sample_positive(s) == 1
)
{

  log_sample(s) =
    log_p +
    log_sample_count(s);
}


// --------------------------------------------------------
// Entire biological sample is zero
// --------------------------------------------------------
//
// There are TWO possibilities:
//
// 1. Capture failed:
//
//        A = 0
//
//        probability = 1-p
//
// 2. Capture succeeded:
//
//        A = 1
//
//    but all PCR observations are zero.
//
// Therefore:
//
// L_ijm
//
// =
//
// (1-p_ijm)
//
// +
//
// p_ijm
// product_k f(0).
//
// f(0) already contains any ZIP/ZINB structural-zero
// probability pi_m.
// --------------------------------------------------------


else
{

  log_sample(s) =
    logspace_add(

      log_1mp,

      log_p +
      log_sample_count(s)

    );
}

}



// ==========================================================
// 16. SITE x OTU CONDITIONAL LIKELIHOOD GIVEN Z = 1
// ==========================================================
//
// Combine all biological sample likelihoods belonging
// to the same site x OTU:
//
//     L(data_im | Z_im=1)
//
//     =
//
//     product_j L_ijm.
// ==========================================================

vector<Type> log_site_conditional(
n_site_groups
);

log_site_conditional.setZero();

for (
int s = 0;
s < n_sample_groups;
s++
)
{

int g =
  sample_site_group(s);


log_site_conditional(g) +=
  log_sample(s);

}



// ==========================================================
// 17. OCCUPANCY MARGINALIZATION
// ==========================================================

vector<Type> psi(
n_site_groups
);

vector<Type> log_site(
n_site_groups
);

for (
int g = 0;
g < n_site_groups;
g++
)
{

// --------------------------------------------------------
// Occupancy probability
// --------------------------------------------------------

psi(g) =
  Type(1) /
  (
    Type(1) +
    exp(
      -eta_occ(g)
    )
  );


Type log_psi =
  log_invlogit(
    eta_occ(g)
  );


Type log_1mpsi =
  log1m_invlogit(
    eta_occ(g)
  );


// --------------------------------------------------------
// Positive site x OTU history
// --------------------------------------------------------
//
// If any observation is positive:
//
//     Z_im = 1
//
// is required.
//
// Therefore:
//
// L_im
//
// =
//
// psi_im
//
// *
//
// L(data_im | Z_im=1).
// --------------------------------------------------------


if (
  site_positive(g) == 1
)
{

  log_site(g) =
    log_psi +
    log_site_conditional(g);
}


// --------------------------------------------------------
// Entire site x OTU history is zero
// --------------------------------------------------------
//
// Two explanations are possible:
//
// 1. Species absent:
//
//        Z_im = 0
//
//        probability = 1-psi
//
// 2. Species present:
//
//        Z_im = 1
//
//    but all biological samples/PCR reads are zero.
//
// Therefore:
//
// L_im
//
// =
//
// (1-psi_im)
//
// +
//
// psi_im
// L(all-zero data | Z_im=1).
// --------------------------------------------------------


else
{

  log_site(g) =
    logspace_add(

      log_1mpsi,

      log_psi +
      log_site_conditional(g)

    );
}


// --------------------------------------------------------
// Add site x OTU contribution to total NLL
// --------------------------------------------------------

nll -=
  log_site(g);

}



// ==========================================================
// 18. DERIVED QUANTITIES
// ==========================================================

vector<Type> capture_prob(
n_sample_groups
);

for (
int s = 0;
s < n_sample_groups;
s++
)
{

capture_prob(s) =
  Type(1) /
  (
    Type(1) +
    exp(
      -eta_cap(s)
    )
  );

// p_ijm =
// invlogit(eta_cap_ijm).

}



vector<Type> lambda(
y.size()
);

for (
int r = 0;
r < y.size();
r++
)
{

lambda(r) =
  exp(
    eta_abund(r)
  );

// lambda_ijkm =
// exp(eta_abund_ijkm).

}



// ==========================================================
// 19. REPORT RESULTS TO R
// ==========================================================

REPORT(
eta_occ
);

REPORT(
eta_cap
);

REPORT(
eta_abund
);

REPORT(
psi
);

REPORT(
capture_prob
);

REPORT(
lambda
);

REPORT(
pi_zi_otu
);

// NEW:
//
// Returns estimated structural-zero probability for every OTU.
//
// Example in R:
//
//     fit$report$pi_zi_otu

REPORT(
log_sample
);

REPORT(
log_site
);



// ==========================================================
// 20. DELTA-METHOD DERIVED-PARAMETER STANDARD ERRORS
// ==========================================================

ADREPORT(
psi
);

// Approximate SEs for occupancy probabilities.

ADREPORT(
capture_prob
);

// Approximate SEs for capture probabilities.

ADREPORT(
pi_zi_otu
);

// NEW:
//
// Approximate SEs for OTU-specific zero-inflation
// probabilities using TMB automatic differentiation.

// ==========================================================
// Return negative log likelihood
// ==========================================================

return nll;
}
