#include <TMB.hpp>

// ============================================================
// Numerically stable logistic helper functions
// ============================================================
//
// TMB provides logspace_add(a,b), which evaluates
//
//     log(exp(a) + exp(b))
//
// without first exponentiating a and b.
// This avoids numerical underflow/overflow when probabilities
// are extremely small or large.
//
// These helper functions are useful because occupancy and
// capture probabilities are modeled on the logit scale.
// ============================================================


template<class Type>
Type log_invlogit(Type x)
{
  // invlogit(x) = 1 / (1 + exp(-x))
  //
  // Therefore:
  //
  // log(invlogit(x))
  // = -log(1 + exp(-x))
  //
  // logspace_add(0,-x)
  // = log(exp(0) + exp(-x))
  // = log(1 + exp(-x))
  //
  // so this returns log(p) in a numerically stable way.

  return -logspace_add(
    Type(0),
    -x
  );
}


template<class Type>
Type log1m_invlogit(Type x)
{
  // Calculates
  //
  // log(1 - invlogit(x)).
  //
  // Because
  //
  // 1 - invlogit(x)
  // = 1 / (1 + exp(x)),
  //
  // we have
  //
  // log(1-p)
  // = -log(1 + exp(x)).

  return -logspace_add(
    Type(0),
    x
  );
}


// ============================================================
// Joint observed-data eDNA likelihood
// ============================================================
//
// The model has THREE hierarchical levels:
//
// ------------------------------------------------------------
// 1. SITE x OTU LEVEL
//
//     Z_im ~ Bernoulli(psi_im)
//
// Z_im is the latent occupancy state.
// There is ONE Z for every site i and OTU m.
//
// ------------------------------------------------------------
// 2. BIOLOGICAL SAMPLE x OTU LEVEL
//
//     A_ijm | Z_im = 1 ~ Bernoulli(p_ijm)
//
// A_ijm represents successful capture of DNA in biological
// sample j, conditional on the OTU being present at site i.
//
// ------------------------------------------------------------
// 3. PCR / READ LEVEL
//
//     Y_ijkm | A_ijm = 1 ~ Count(lambda_ijkm)
//
// Y is the observed sequencing read count.
//
// ------------------------------------------------------------
//
// IMPORTANT:
//
// Z and A are NOT optimized as binary parameters.
//
// Instead:
//
//     A is analytically marginalized
//
// and then
//
//     Z is analytically marginalized.
//
// The resulting likelihood is therefore an observed-data
// likelihood conditional only on model parameters and any
// Gaussian random effects.
// ============================================================


template<class Type>
Type objective_function<Type>::operator() ()
{

  // ==========================================================
  // DATA
  // ==========================================================


  DATA_VECTOR(y);

  // y contains the observed sequencing counts.
  //
  // Each element:
  //
  //     y(r)
  //
  // corresponds to one PCR/read-level observation.
  //
  // Thus r indexes rows of the full long-format dataset.


  // ----------------------------------------------------------
  // Fixed-effect design matrices
  // ----------------------------------------------------------


  DATA_MATRIX(X_occ);

  // Design matrix for occupancy.
  //
  // IMPORTANT:
  //
  // rows correspond to SITE x OTU combinations,
  // NOT PCR/read observations.
  //
  // If there are G unique site x OTU combinations,
  //
  //     nrow(X_occ) = G.
  //
  // The occupancy linear predictor will be
  //
  //     eta_occ = X_occ * beta_occ.


  DATA_MATRIX(X_cap);

  // Design matrix for capture.
  //
  // Rows correspond to:
  //
  //     biological sample x OTU
  //
  // combinations.
  //
  // Thus capture lives one level below occupancy.


  DATA_MATRIX(X_abund);

  // Design matrix for abundance.
  //
  // Rows correspond to individual PCR/read count observations.
  //
  // This is the lowest observational level of the hierarchy.


  DATA_VECTOR(offset_abund);

  // Optional abundance offset.
  //
  // Typically this may be
  //
  //     log(total_reads)
  //
  // to account for sequencing depth.
  //
  // It enters the linear predictor with coefficient fixed at 1.


  // ----------------------------------------------------------
  // Mapping/index vectors
  // ----------------------------------------------------------


  DATA_IVECTOR(occ_otu);

  // For every SITE x OTU occupancy row g,
  // occ_otu(g) identifies which OTU that row belongs to.
  //
  // This is used to attach an OTU-specific occupancy
  // random effect.


  DATA_IVECTOR(sample_otu);

  // For every biological SAMPLE x OTU group s,
  // sample_otu(s) identifies the corresponding OTU.
  //
  // Used for the capture OTU random effect.


  DATA_IVECTOR(row_sample_group);

  // Maps every read/PCR row r to its SAMPLE x OTU group s.
  //
  // Symbolically:
  //
  //     r -> s
  //
  // This mapping allows all PCR replicates from one
  // biological sample and OTU to be combined.


  DATA_IVECTOR(sample_site_group);

  // Maps each SAMPLE x OTU group s to its SITE x OTU group g.
  //
  // Symbolically:
  //
  //     s -> g
  //
  // This is what connects capture to occupancy.
  //
  // Consequently:
  //
  // multiple biological samples from the same site share
  // the SAME site x OTU occupancy state.


  DATA_IVECTOR(row_otu);

  // Maps every read/PCR row r to the OTU represented by
  // that observation.
  //
  // Used for abundance OTU random effects.


  DATA_IVECTOR(row_sample_id);

  // Maps every read/PCR row to its biological sample.
  //
  // Used for the sample-level abundance random effect.


  DATA_IVECTOR(row_sample_otu_re);

  // Maps each PCR/read observation to its unique
  // sample x OTU random-effect level.
  //
  // This allows extra heterogeneity specific to a particular
  // OTU within a particular biological sample.


  DATA_IVECTOR(sample_positive);

  // Binary indicator:
  //
  //     sample_positive(s) = 1
  //
  // if at least one PCR/read belonging to sample x OTU
  // group s has a positive count.
  //
  // If this equals 1, then capture A MUST equal 1.


  DATA_IVECTOR(site_positive);

  // Binary indicator:
  //
  //     site_positive(g) = 1
  //
  // if at least one observation anywhere within site x OTU
  // group g is positive.
  //
  // If this equals 1, occupancy Z MUST equal 1.


  DATA_INTEGER(n_site_groups);

  // Number of unique SITE x OTU combinations.
  //
  // This determines the length of psi and eta_occ.


  DATA_INTEGER(n_sample_groups);

  // Number of unique biological SAMPLE x OTU combinations.
  //
  // This determines the length of p and eta_cap.


  // ----------------------------------------------------------
  // Count family
  // ----------------------------------------------------------
  //
  // Controls the distribution of sequencing counts:
  //
  //     family_code = 0  -> Poisson
  //     family_code = 1  -> Negative Binomial 2
  //     family_code = 2  -> Zero-Inflated Poisson
  //     family_code = 3  -> Zero-Inflated NB2
  // ----------------------------------------------------------


  DATA_INTEGER(family_code);


  // ----------------------------------------------------------
  // Random-effect switches
  // ----------------------------------------------------------
  //
  // Each variable is 0 or 1.
  //
  // These correspond to the R arguments such as
  //
  // random_occ_otu = TRUE/FALSE
  //
  // etc.
  // ----------------------------------------------------------


  DATA_INTEGER(use_occ_otu);

  // 1 -> include OTU random effect in occupancy.


  DATA_INTEGER(use_cap_otu);

  // 1 -> include OTU random effect in capture.


  DATA_INTEGER(use_abund_otu);

  // 1 -> include OTU random effect in abundance.


  DATA_INTEGER(use_sample_re);

  // 1 -> include biological-sample random effect in abundance.


  DATA_INTEGER(use_sample_otu_re);

  // 1 -> include sample x OTU random effect in abundance.



  // ==========================================================
  // FIXED EFFECTS
  // ==========================================================


  PARAMETER_VECTOR(beta_occ);

  // Fixed-effect regression coefficients for occupancy.
  //
  // Mathematically:
  //
  // eta_occ = X_occ beta_occ + random effects.


  PARAMETER_VECTOR(beta_cap);

  // Fixed-effect coefficients for capture.


  PARAMETER_VECTOR(beta_abund);

  // Fixed-effect coefficients for abundance.



  // ==========================================================
  // RANDOM EFFECTS
  // ==========================================================


  PARAMETER_VECTOR(b_occ_otu);

  // OTU-specific occupancy effects.
  //
  // Conceptually:
  //
  // b_occ_otu[m] ~ Normal(0, sigma_occ_otu^2)
  //
  // This allows some OTUs to have systematically higher or
  // lower occupancy than the population-level fixed effects.


  PARAMETER_VECTOR(b_cap_otu);

  // OTU-specific capture random effects.


  PARAMETER_VECTOR(b_abund_otu);

  // OTU-specific abundance random effects.


  PARAMETER_VECTOR(b_sample);

  // Biological-sample random effects in abundance.
  //
  // Captures sample-to-sample variation shared by all OTUs
  // within a sample.


  PARAMETER_VECTOR(b_sample_otu);

  // Sample x OTU random effect.
  //
  // Represents additional variability specific to one OTU
  // in one biological sample.



  // ==========================================================
  // RANDOM-EFFECT STANDARD DEVIATIONS
  // ==========================================================


  PARAMETER(log_sd_occ_otu);

  // log standard deviation of occupancy OTU random effects.


  PARAMETER(log_sd_cap_otu);

  // log standard deviation of capture OTU random effects.


  PARAMETER(log_sd_abund_otu);

  // log standard deviation of abundance OTU random effects.


  PARAMETER(log_sd_sample);

  // log standard deviation of biological-sample effects.


  PARAMETER(log_sd_sample_otu);

  // log standard deviation of sample x OTU effects.


  Type sd_occ_otu =
    exp(log_sd_occ_otu);

  // Standard deviation must be positive.
  //
  // Therefore the optimizer works with an unconstrained
  // log standard deviation and transforms back using exp().


  Type sd_cap_otu =
    exp(log_sd_cap_otu);


  Type sd_abund_otu =
    exp(log_sd_abund_otu);


  Type sd_sample =
    exp(log_sd_sample);


  Type sd_sample_otu =
    exp(log_sd_sample_otu);



  // ==========================================================
  // ABUNDANCE DISTRIBUTION PARAMETERS
  // ==========================================================


  PARAMETER(log_theta);

  // Log NB2 dispersion parameter.
  //
  // Only relevant for NB and ZINB.


  Type theta =
    exp(log_theta);

  // Forces theta > 0.


  PARAMETER(zi_intercept);

  // Zero-inflation intercept used by ZIP/ZINB.


  Type pi_zi =
    Type(1) /
    (
      Type(1) +
      exp(-zi_intercept)
    );

  // Converts zero-inflation intercept to probability:
  //
  //     pi_zi = invlogit(zi_intercept)
  //
  // pi_zi is the probability of belonging to the
  // structural-zero component.



  // ==========================================================
  // NEGATIVE LOG-LIKELIHOOD
  // ==========================================================


  Type nll =
    Type(0);

  // TMB minimizes an objective function.
  //
  // We therefore accumulate:
  //
  //     nll = - log L.
  //
  // Initial value is zero.



  // ==========================================================
  // RANDOM-EFFECT GAUSSIAN DENSITIES
  // ==========================================================
  //
  // Each enabled random effect receives a Gaussian density.
  //
  // These terms are part of the joint likelihood used by TMB.
  //
  // If the effects are declared as random in MakeADFun(),
  // TMB later integrates them using the Laplace approximation.
  // ==========================================================


  if (use_occ_otu == 1)
  {

    // Only evaluate this density if the occupancy OTU
    // random effect is enabled.


    for (int j = 0;
         j < b_occ_otu.size();
         j++)
    {

      // Loop through all OTU occupancy random effects.


      nll -= dnorm(
        b_occ_otu(j),
        Type(0),
        sd_occ_otu,
        true
      );

      // dnorm(..., true) returns:
      //
      //     log Normal(
      //         b_occ_otu[j];
      //         mean = 0,
      //         sd   = sd_occ_otu
      //     )
      //
      // We subtract it because nll is NEGATIVE log likelihood.
    }
  }



  if (use_cap_otu == 1)
  {

    for (int j = 0;
         j < b_cap_otu.size();
         j++)
    {

      nll -= dnorm(
        b_cap_otu(j),
        Type(0),
        sd_cap_otu,
        true
      );

      // Assumption:
      //
      // b_cap_otu[j]
      // ~ Normal(0, sd_cap_otu^2).
    }
  }



  if (use_abund_otu == 1)
  {

    for (int j = 0;
         j < b_abund_otu.size();
         j++)
    {

      nll -= dnorm(
        b_abund_otu(j),
        Type(0),
        sd_abund_otu,
        true
      );

      // OTU-specific abundance effect:
      //
      // b_abund_otu[j]
      // ~ Normal(0, sd_abund_otu^2).
    }
  }



  if (use_sample_re == 1)
  {

    for (int j = 0;
         j < b_sample.size();
         j++)
    {

      nll -= dnorm(
        b_sample(j),
        Type(0),
        sd_sample,
        true
      );

      // Biological sample effects:
      //
      // b_sample[j]
      // ~ Normal(0, sd_sample^2).
    }
  }



  if (use_sample_otu_re == 1)
  {

    for (int j = 0;
         j < b_sample_otu.size();
         j++)
    {

      nll -= dnorm(
        b_sample_otu(j),
        Type(0),
        sd_sample_otu,
        true
      );

      // Sample x OTU effects:
      //
      // b_sample_otu[j]
      // ~ Normal(0, sd_sample_otu^2).
    }
  }



  // ==========================================================
  // OCCUPANCY LINEAR PREDICTOR
  // ==========================================================


  vector<Type> eta_occ =
    X_occ * beta_occ;

  // Start with fixed effects:
  //
  //     eta_occ(g)
  //     =
  //     X_occ(g,.) beta_occ
  //
  // g indexes SITE x OTU combinations.
  //
  // Therefore occupancy is explicitly represented at
  // the SITE x OTU level.


  if (use_occ_otu == 1)
  {

    for (int g = 0;
         g < n_site_groups;
         g++)
    {

      eta_occ(g) +=
        b_occ_otu(
          occ_otu(g)
        );

      // Add the random effect of the OTU associated with
      // site x OTU group g.
      //
      // Final predictor:
      //
      // eta_occ_im
      // =
      // X_occ_im beta_occ
      // +
      // b_occ_m.
    }
  }



  // ==========================================================
  // CAPTURE LINEAR PREDICTOR
  // ==========================================================


  vector<Type> eta_cap =
    X_cap * beta_cap;

  // Fixed-effect capture predictor:
  //
  // eta_cap_ijm
  // =
  // X_cap_ijm beta_cap.


  if (use_cap_otu == 1)
  {

    for (int s = 0;
         s < n_sample_groups;
         s++)
    {

      eta_cap(s) +=
        b_cap_otu(
          sample_otu(s)
        );

      // Add OTU-specific capture effect:
      //
      // eta_cap_ijm
      // =
      // X_cap beta_cap
      // +
      // b_cap_m.
    }
  }



  // ==========================================================
  // ABUNDANCE LINEAR PREDICTOR
  // ==========================================================


  vector<Type> eta_abund =
    X_abund * beta_abund;

  // Fixed-effect abundance predictor.
  //
  // One value for each read/PCR observation r.


  for (int r = 0;
       r < y.size();
       r++)
  {

    eta_abund(r) +=
      offset_abund(r);

    // Add sequencing-depth or other offset.


    if (use_abund_otu == 1)
    {

      eta_abund(r) +=
        b_abund_otu(
          row_otu(r)
        );

      // Add abundance OTU effect.
    }


    if (use_sample_re == 1)
    {

      eta_abund(r) +=
        b_sample(
          row_sample_id(r)
        );

      // Add random intercept belonging to biological sample
      // containing observation r.
    }


    if (use_sample_otu_re == 1)
    {

      eta_abund(r) +=
        b_sample_otu(
          row_sample_otu_re(r)
        );

      // Add interaction-specific sample x OTU effect.
    }
  }



  // ==========================================================
  // READ-LEVEL COUNT LOG LIKELIHOOD
  // ==========================================================


  vector<Type> log_count(
    y.size()
  );

  // Storage vector:
  //
  // log_count(r)
  //
  // will contain the log probability of read/count observation r
  // conditional on capture A = 1.



  for (int r = 0;
       r < y.size();
       r++)
  {

    Type mu =
      exp(
        eta_abund(r)
      );

    // Abundance model uses a log link:
    //
    //     lambda_r = mu_r
    //              = exp(eta_abund_r)
    //
    // mu must therefore always be positive.


    Type lp =
      Type(0);

    // Temporary variable containing the log probability
    // for observation r.



    // --------------------------------------------------------
    // Poisson
    // --------------------------------------------------------


    if (family_code == 0)
    {

      lp =
        dpois(
          y(r),
          mu,
          true
        );

      // Evaluates:
      //
      //     log P(Y_r = y_r | mu_r)
      //
      // for a Poisson distribution.
      //
      // true means return the LOG probability.
    }



    // --------------------------------------------------------
    // Negative Binomial 2
    //
    // E(Y)   = mu
    //
    // Var(Y) = mu + mu^2/theta
    // --------------------------------------------------------


    else if (family_code == 1)
    {

      Type variance =
        mu +
        mu * mu /
        theta;

      // Compute NB2 variance.


      lp =
        dnbinom2(
          y(r),
          mu,
          variance,
          true
        );

      // NB2 log probability of observed count y(r).
    }



    // --------------------------------------------------------
    // Zero-Inflated Poisson
    // --------------------------------------------------------


    else if (family_code == 2)
    {

      Type log_count_component =
        dpois(
          y(r),
          mu,
          true
        );

      // Standard Poisson log probability.


      if (y(r) == Type(0))
      {

        // For a zero count there are TWO ways to observe zero:
        //
        // 1. Structural zero:
        //
        //        probability = pi_zi
        //
        // 2. Observation came from the Poisson component
        //    and Poisson generated zero:
        //
        //        probability =
        //        (1-pi_zi) * P_Pois(Y=0)


        lp =
          logspace_add(

            log(pi_zi),

            log(
              Type(1) -
              pi_zi
            ) +
            log_count_component
          );

        // This calculates:
        //
        // log[
        //     pi_zi
        //     +
        //     (1-pi_zi) P_Pois(0)
        // ].
      }


      else
      {

        // Positive observations cannot originate from the
        // structural-zero component.


        lp =
          log(
            Type(1) -
            pi_zi
          ) +
          log_count_component;

        // Therefore:
        //
        // P(Y=y>0)
        // =
        // (1-pi_zi) P_Pois(y).
      }
    }



    // --------------------------------------------------------
    // Zero-Inflated Negative Binomial 2
    // --------------------------------------------------------


    else if (family_code == 3)
    {

      Type variance =
        mu +
        mu * mu /
        theta;

      // NB2 variance.


      Type log_count_component =
        dnbinom2(
          y(r),
          mu,
          variance,
          true
        );

      // NB2 component likelihood.


      if (y(r) == Type(0))
      {

        lp =
          logspace_add(

            log(pi_zi),

            log(
              Type(1) -
              pi_zi
            ) +
            log_count_component
          );

        // For a zero:
        //
        // P(Y=0)
        // =
        // pi_zi
        // +
        // (1-pi_zi) P_NB(0).
      }


      else
      {

        lp =
          log(
            Type(1) -
            pi_zi
          ) +
          log_count_component;

        // For positive counts:
        //
        // P(Y=y)
        // =
        // (1-pi_zi) P_NB(y).
      }
    }


    log_count(r) =
      lp;

    // Store count-level log likelihood.
  }



  // ==========================================================
  // SAMPLE x OTU CONDITIONAL LIKELIHOOD
  // ==========================================================
  //
  // Here we move from PCR/read observations to biological
  // samples.
  //
  // First compute:
  //
  //     sum_k log f(y_ijkm | A_ijm = 1)
  //
  // which is equivalent to:
  //
  //     log product_k f(y_ijkm | A_ijm = 1)
  //
  // Then analytically marginalize the latent capture state A.
  // ==========================================================


  vector<Type> log_sample_count(
    n_sample_groups
  );

  // Stores the count likelihood conditional on A = 1
  // for each sample x OTU.


  log_sample_count.setZero();

  // Initialize all sample-level log likelihoods to zero.



  for (int r = 0;
       r < y.size();
       r++)
  {

    int s =
      row_sample_group(r);

    // Find sample x OTU group s that read observation r
    // belongs to.


    log_sample_count(s) +=
      log_count(r);

    // Add PCR log likelihoods within the same biological sample.
    //
    // At the end:
    //
    // log_sample_count(s)
    // =
    // sum_k log f(y_ijkm)
    //
    // =
    // log product_k f(y_ijkm).
  }



  vector<Type> log_sample(
    n_sample_groups
  );

  // Final sample-level likelihood AFTER accounting for
  // uncertainty in capture A.



  for (int s = 0;
       s < n_sample_groups;
       s++)
  {

    Type log_p =
      log_invlogit(
        eta_cap(s)
      );

    // Capture model:
    //
    // p_ijm = invlogit(eta_cap_ijm)
    //
    // This stores:
    //
    // log(p_ijm).


    Type log_1mp =
      log1m_invlogit(
        eta_cap(s)
      );

    // This stores:
    //
    // log(1-p_ijm).



    // --------------------------------------------------------
    // SAMPLE HAS AT LEAST ONE POSITIVE READ
    //
    // If any y > 0, A = 1 is required.
    //
    // L =
    //
    // p * product_k f(y_k)
    // --------------------------------------------------------


    if (sample_positive(s) == 1)
    {

      log_sample(s) =
        log_p +
        log_sample_count(s);

      // On original probability scale:
      //
      // L_ijm
      // =
      // p_ijm
      // *
      // product_k f(y_ijkm | A_ijm=1).
      //
      // There is NO contribution from A=0 because A=0 would
      // force all counts to equal zero.
    }



    // --------------------------------------------------------
    // ENTIRE BIOLOGICAL SAMPLE IS ZERO
    //
    // A is unknown and therefore marginalized.
    //
    // Two possibilities:
    //
    // A = 0:
    //
    //     probability = 1-p
    //
    // A = 1 but every PCR produced zero:
    //
    //     probability =
    //     p * product_k f(0)
    //
    // Therefore:
    //
    // L =
    //
    // (1-p)
    // +
    // p * product_k f(0)
    // --------------------------------------------------------


    else
    {

      log_sample(s) =
        logspace_add(

          log_1mp,

          log_p +
          log_sample_count(s)

        );

      // This is the CAPTURE MARGINALIZATION.
      //
      // Mathematically:
      //
      // L_ijm
      //
      // =
      //
      // SUM over A_ijm in {0,1}
      //
      // P(A_ijm | Z_im=1)
      // P(Y_ijm | A_ijm)
      //
      // =
      //
      // (1-p_ijm)
      //
      // +
      //
      // p_ijm
      // product_k f(0 | lambda_ijkm).
      //
      // A is therefore not sampled and not explicitly
      // estimated by the optimizer.
    }
  }



  // ==========================================================
  // SITE x OTU CONDITIONAL LIKELIHOOD GIVEN Z = 1
  // ==========================================================
  //
  // We now move from biological samples to SITE x OTU.
  //
  // For site i and OTU m:
  //
  // L(data_im | Z_im=1)
  //
  // =
  //
  // product_j L_ijm
  //
  // where j indexes biological samples.
  // ==========================================================


  vector<Type> log_site_conditional(
    n_site_groups
  );

  // Stores:
  //
  // log L(data_im | Z_im = 1).


  log_site_conditional.setZero();



  for (int s = 0;
       s < n_sample_groups;
       s++)
  {

    int g =
      sample_site_group(s);

    // Find site x OTU group g associated with sample group s.


    log_site_conditional(g) +=
      log_sample(s);

    // Sum sample log likelihoods:
    //
    // sum_j log L_ijm
    //
    // =
    //
    // log product_j L_ijm.
  }



  // ==========================================================
  // OCCUPANCY MARGINALIZATION
  // ==========================================================
  //
  // We now marginalize the latent site-level state:
  //
  //     Z_im in {0,1}.
  //
  // This occurs ONCE per SITE x OTU combination.
  //
  // It does NOT occur independently for each PCR observation.
  // ==========================================================


  vector<Type> psi(
    n_site_groups
  );

  // Natural-scale occupancy probabilities.


  vector<Type> log_site(
    n_site_groups
  );

  // Final marginal log likelihood contribution for every
  // site x OTU combination.



  for (int g = 0;
       g < n_site_groups;
       g++)
  {

    psi(g) =
      Type(1) /
      (
        Type(1) +
        exp(
          -eta_occ(g)
        )
      );

    // Logistic inverse link:
    //
    // psi_im
    // =
    // P(Z_im = 1)
    // =
    // invlogit(eta_occ_im).


    Type log_psi =
      log_invlogit(
        eta_occ(g)
      );

    // log(psi_im).


    Type log_1mpsi =
      log1m_invlogit(
        eta_occ(g)
      );

    // log(1-psi_im).



    // --------------------------------------------------------
    // POSITIVE SITE HISTORY
    //
    // If at least one observation is positive,
    //
    // Z_im = 1 is required.
    //
    // L =
    //
    // psi *
    // L(data | Z=1)
    // --------------------------------------------------------


    if (site_positive(g) == 1)
    {

      log_site(g) =
        log_psi +
        log_site_conditional(g);

      // Original scale:
      //
      // L_im
      // =
      //
      // psi_im
      // *
      // L(data_im | Z_im=1).
      //
      // Z=0 is impossible because a positive observation
      // cannot occur when the taxon is absent.
    }



    // --------------------------------------------------------
    // ENTIRE SITE x OTU HISTORY IS ZERO
    //
    // Z is unknown.
    //
    // Two possible explanations:
    //
    // Z = 0:
    //
    //     probability = 1-psi
    //
    // OR
    //
    // Z = 1, but all samples/PCR observations are zero:
    //
    //     probability =
    //     psi * L(all-zero data | Z=1)
    //
    // Therefore:
    //
    // L =
    //
    // (1-psi)
    // +
    // psi * L(all-zero data | Z=1)
    // --------------------------------------------------------


    else
    {

      log_site(g) =
        logspace_add(

          log_1mpsi,

          log_psi +
          log_site_conditional(g)

        );

      // This is the OCCUPANCY MARGINALIZATION.
      //
      // Mathematically:
      //
      // L_im
      //
      // =
      //
      // SUM over Z_im in {0,1}
      //
      // P(Z_im)
      // P(Y_im | Z_im)
      //
      // =
      //
      // (1-psi_im)
      //
      // +
      //
      // psi_im
      // L(Y_im=0 | Z_im=1).
      //
      // Thus Z_im is not sampled or optimized directly.
    }



    nll -=
      log_site(g);

    // Add SITE x OTU likelihood contribution to the
    // negative log likelihood.
    //
    // Since:
    //
    // log L =
    // sum_g log_site(g),
    //
    // the objective minimized by TMB is:
    //
    // nll =
    // - sum_g log_site(g).
  }



  // ==========================================================
  // DERIVED QUANTITIES
  // ==========================================================


  vector<Type> capture_prob(
    n_sample_groups
  );

  // Natural-scale capture probabilities p_ijm.



  for (int s = 0;
       s < n_sample_groups;
       s++)
  {

    capture_prob(s) =
      Type(1) /
      (
        Type(1) +
        exp(
          -eta_cap(s)
        )
      );

    // Transform:
    //
    // p_ijm
    // =
    // invlogit(eta_cap_ijm).
  }



  vector<Type> lambda(
    y.size()
  );

  // Natural-scale expected sequencing abundance.



  for (int r = 0;
       r < y.size();
       r++)
  {

    lambda(r) =
      exp(
        eta_abund(r)
      );

    // Because abundance uses a log link:
    //
    // lambda_ijkm
    // =
    // exp(eta_abund_ijkm).
  }



  // ==========================================================
  // REPORT
  // ==========================================================
  //
  // REPORT() makes quantities available in R through:
  //
  //     obj$report()
  //
  // It does not automatically request delta-method SEs.
  // ==========================================================


  REPORT(
    eta_occ
  );

  // Occupancy linear predictor.


  REPORT(
    eta_cap
  );

  // Capture linear predictor.


  REPORT(
    eta_abund
  );

  // Abundance linear predictor.


  REPORT(
    psi
  );

  // Natural-scale occupancy probabilities.


  REPORT(
    capture_prob
  );

  // Natural-scale capture probabilities.


  REPORT(
    lambda
  );

  // Natural-scale abundance means.


  REPORT(
    log_sample
  );

  // Marginal sample x OTU log-likelihood contributions,
  // after integrating over capture A.


  REPORT(
    log_site
  );

  // Final marginal site x OTU log-likelihood contributions,
  // after integrating over occupancy Z.



  // ==========================================================
  // DELTA-METHOD DERIVED-PARAMETER STANDARD ERRORS
  // ==========================================================
  //
  // ADREPORT tells TMB to propagate parameter uncertainty
  // through a derived quantity using automatic differentiation
  // and the delta method.
  // ==========================================================


  ADREPORT(
    psi
  );

  // Allows sdreport() to calculate approximate SEs for
  // occupancy probabilities psi.


  ADREPORT(
    capture_prob
  );

  // Allows sdreport() to calculate approximate SEs for
  // capture probabilities p.


  return nll;

  // Return the total negative log likelihood.
  //
  // nlminb()/TMB minimizes this quantity to estimate the
  // fixed effects and variance/dispersion parameters.
}
