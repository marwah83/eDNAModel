#include <TMB.hpp>


// ============================================================
// Numerically stable logistic helper functions
// ============================================================
//
// For a logistic probability
//
//     p = invlogit(x) = 1 / (1 + exp(-x))
//
// we frequently need:
//
//     log(p)
//     log(1 - p)
//
// Computing p first and then taking log(p) or log(1-p) can
// become numerically unstable when x is very large or very
// negative.
//
// These functions calculate the quantities directly on the
// log scale.
// ============================================================


template<class Type>
Type log_invlogit(Type x)
{
  // log(invlogit(x))
  //
  // = -log(1 + exp(-x))
  //
  // logspace_add() provides a numerically stable evaluation.

  return -logspace_add(
    Type(0),
    -x
  );
}


template<class Type>
Type log1m_invlogit(Type x)
{
  // log(1 - invlogit(x))
  //
  // = -log(1 + exp(x))

  return -logspace_add(
    Type(0),
    x
  );
}



// ============================================================
// Joint observed-data eDNA likelihood
// ============================================================
//
// HIERARCHICAL MODEL
//
// ------------------------------------------------------------
// LEVEL 1: SITE x OTU OCCUPANCY
// ------------------------------------------------------------
//
//     Z_im ~ Bernoulli(psi_im)
//
// where:
//
//     i = site
//     m = OTU
//
// There is ONE latent occupancy state for each site x OTU.
//
// The occupancy linear predictor is:
//
//     logit(psi_im)
//       = X_occ beta_occ
//         + b_occ_m
//
// when the OTU random effect is enabled.
//
// ------------------------------------------------------------
// LEVEL 2: BIOLOGICAL SAMPLE x OTU CAPTURE
// ------------------------------------------------------------
//
// Conditional on occupancy:
//
//     A_ijm | Z_im = 1
//       ~ Bernoulli(p_ijm)
//
// where:
//
//     j = biological sample.
//
// The capture linear predictor is:
//
//     logit(p_ijm)
//       = X_cap beta_cap
//         + b_cap_m
//
// when the OTU random effect is enabled.
//
// ------------------------------------------------------------
// LEVEL 3: PCR / READ COUNTS
// ------------------------------------------------------------
//
// Conditional on successful capture:
//
//     Y_ijkm | A_ijm = 1
//       ~ Count(lambda_ijkm)
//
// where:
//
//     k = PCR replicate.
//
// The abundance linear predictor is:
//
//     log(lambda_ijkm)
//       = X_abund beta_abund
//         + offset
//         + b_abund_m
//         + b_sample_j
//         + b_sampleOTU_jm
//
// depending on which random effects are enabled.
//
// Supported distributions:
//
//     family_code = 0 : Poisson
//     family_code = 1 : NB2
//     family_code = 2 : ZIP
//     family_code = 3 : ZINB
//
// ------------------------------------------------------------
// ZERO-INFLATION MODEL
// ------------------------------------------------------------
//
// For ZIP/ZINB:
//
//     logit(pi_m)
//       = zi_intercept
//         + b_zi_m
//
// where:
//
//     b_zi_m ~ Normal(0, sigma_zi^2)
//
// Therefore the structural-zero probability is allowed to
// differ among OTUs.
//
// If use_zi_otu = 0:
//
//     logit(pi) = zi_intercept
//
// and all OTUs share the same structural-zero probability.
//
// ------------------------------------------------------------
// MARGINALIZATION
// ------------------------------------------------------------
//
// The discrete latent states Z and A are NOT optimized as
// parameters.
//
// Instead:
//
//     A is analytically marginalized at the sample x OTU level.
//
// Then:
//
//     Z is analytically marginalized at the site x OTU level.
//
// Gaussian random effects are integrated by TMB through the
// Laplace approximation when supplied in the `random` argument
// of MakeADFun().
// ============================================================


template<class Type>
Type objective_function<Type>::operator() ()
{

  // ==========================================================
  // 1. DATA
  // ==========================================================


  DATA_VECTOR(y);

  // PCR/read-level observed counts.
  //
  // y(r) is the count associated with observation row r.



  // ----------------------------------------------------------
  // Design matrices
  // ----------------------------------------------------------


  DATA_MATRIX(X_occ);

  // Occupancy design matrix.
  //
  // IMPORTANT:
  // One row per SITE x OTU combination.
  //
  // Occupancy therefore does NOT live at the PCR/read level.


  DATA_MATRIX(X_cap);

  // Capture design matrix.
  //
  // One row per BIOLOGICAL SAMPLE x OTU combination.


  DATA_MATRIX(X_abund);

  // Abundance design matrix.
  //
  // One row per PCR/read observation.



  // ----------------------------------------------------------
  // Abundance offset
  // ----------------------------------------------------------


  DATA_VECTOR(offset_abund);

  // Optional abundance offset.
  //
  // For example:
  //
  //     log(total_reads)
  //
  // for sequencing-depth adjustment.



  // ----------------------------------------------------------
  // Hierarchical mapping indices
  // ----------------------------------------------------------


  DATA_IVECTOR(occ_otu);

  // Maps SITE x OTU occupancy group g
  // to OTU m.
  //
  // Used for:
  //
  //     b_occ_otu(m)


  DATA_IVECTOR(sample_otu);

  // Maps SAMPLE x OTU group s
  // to OTU m.
  //
  // Used for:
  //
  //     b_cap_otu(m)


  DATA_IVECTOR(row_sample_group);

  // Maps PCR/read row r
  // to SAMPLE x OTU group s.


  DATA_IVECTOR(sample_site_group);

  // Maps SAMPLE x OTU group s
  // to SITE x OTU group g.
  //
  // This is what makes all biological samples from the
  // same site x OTU share the same occupancy state.


  DATA_IVECTOR(row_otu);

  // Maps PCR/read row r
  // to OTU m.
  //
  // Used by:
  //
  //     b_abund_otu(m)
  //     b_zi_otu(m)


  DATA_IVECTOR(row_sample_id);

  // Maps PCR/read row r
  // to biological sample j.


  DATA_IVECTOR(row_sample_otu_re);

  // Maps PCR/read row r
  // to the sample x OTU random-effect level.



  // ----------------------------------------------------------
  // Positive-history indicators
  // ----------------------------------------------------------


  DATA_IVECTOR(sample_positive);

  // sample_positive(s) = 1 when at least one PCR replicate
  // within sample x OTU group s has y > 0.
  //
  // In that case:
  //
  //     A_ijm = 1
  //
  // must hold.


  DATA_IVECTOR(site_positive);

  // site_positive(g) = 1 when at least one count anywhere
  // within site x OTU group g is positive.
  //
  // In that case:
  //
  //     Z_im = 1
  //
  // must hold.



  // ----------------------------------------------------------
  // Number of hierarchical groups
  // ----------------------------------------------------------


  DATA_INTEGER(n_site_groups);

  // Number of SITE x OTU groups.


  DATA_INTEGER(n_sample_groups);

  // Number of SAMPLE x OTU groups.



  // ----------------------------------------------------------
  // Abundance distribution
  // ----------------------------------------------------------


  DATA_INTEGER(family_code);

  // family_code:
  //
  // 0 = Poisson
  // 1 = NB2
  // 2 = ZIP
  // 3 = ZINB



  // ==========================================================
  // 2. RANDOM-EFFECT SWITCHES
  // ==========================================================


  DATA_INTEGER(use_occ_otu);

  // 1 = use OTU occupancy random intercept.


  DATA_INTEGER(use_cap_otu);

  // 1 = use OTU capture random intercept.


  DATA_INTEGER(use_abund_otu);

  // 1 = use OTU abundance random intercept.


  DATA_INTEGER(use_sample_re);

  // 1 = use biological-sample abundance random intercept.


  DATA_INTEGER(use_sample_otu_re);

  // 1 = use sample x OTU abundance random intercept.


  DATA_INTEGER(use_zi_otu);

  // 1 = use OTU-specific zero-inflation random intercept.
  //
  // Relevant only for ZIP/ZINB.



  // ==========================================================
  // 3. VALIDATE FAMILY CODE
  // ==========================================================


  if (
    family_code < 0 ||
    family_code > 3
  )
  {
    error(
      "family_code must be 0, 1, 2, or 3."
    );
  }



  // ==========================================================
  // 4. FIXED-EFFECT PARAMETERS
  // ==========================================================


  PARAMETER_VECTOR(beta_occ);

  // Occupancy fixed effects.


  PARAMETER_VECTOR(beta_cap);

  // Capture fixed effects.


  PARAMETER_VECTOR(beta_abund);

  // Abundance fixed effects.



  // ==========================================================
  // 5. RANDOM-EFFECT VECTORS
  // ==========================================================


  PARAMETER_VECTOR(b_occ_otu);

  // OTU-specific occupancy effects:
  //
  //     b_occ_m ~ N(0, sigma_occ^2)


  PARAMETER_VECTOR(b_cap_otu);

  // OTU-specific capture effects:
  //
  //     b_cap_m ~ N(0, sigma_cap^2)


  PARAMETER_VECTOR(b_abund_otu);

  // OTU-specific abundance effects:
  //
  //     b_abund_m ~ N(0, sigma_abund^2)


  PARAMETER_VECTOR(b_sample);

  // Biological-sample abundance effects.


  PARAMETER_VECTOR(b_sample_otu);

  // Sample x OTU abundance effects.


  PARAMETER_VECTOR(b_zi_otu);

  // OTU-specific zero-inflation effects:
  //
  //     b_zi_m ~ N(0, sigma_zi^2)



  // ==========================================================
  // 6. RANDOM-EFFECT STANDARD DEVIATIONS
  // ==========================================================


  PARAMETER(log_sd_occ_otu);

  PARAMETER(log_sd_cap_otu);

  PARAMETER(log_sd_abund_otu);

  PARAMETER(log_sd_sample);

  PARAMETER(log_sd_sample_otu);

  PARAMETER(log_sd_zi_otu);



  // Transform log standard deviations to positive scale.

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



  // ==========================================================
  // 7. DISTRIBUTION PARAMETERS
  // ==========================================================


  PARAMETER(log_theta);

  // NB2 dispersion parameter.
  //
  // Currently ONE theta is shared across all OTUs.


  Type theta =
    exp(
      log_theta
    );

  // NB2 parameterization:
  //
  //     E(Y) = mu
  //
  //     Var(Y)
  //       = mu + mu^2 / theta



  PARAMETER(zi_intercept);

  // Population-level zero-inflation intercept.
  //
  // For ZIP/ZINB:
  //
  //     logit(pi_m)
  //       =
  //       zi_intercept
  //       +
  //       b_zi_otu(m)
  //
  // when use_zi_otu = 1.



  // ==========================================================
  // 8. NEGATIVE LOG-LIKELIHOOD
  // ==========================================================


  Type nll =
    Type(0);



  // ==========================================================
  // 9. RANDOM-EFFECT DISTRIBUTIONS
  // ==========================================================
  //
  // All random-effect blocks are currently assumed independent.
  //
  // There are NO correlations between:
  //
  //   occupancy OTU effects
  //   capture OTU effects
  //   abundance OTU effects
  //   zero-inflation OTU effects
  //
  // or the sample-level abundance effects.
  // ==========================================================


  // ----------------------------------------------------------
  // Occupancy OTU random effects
  // ----------------------------------------------------------


  if (
    use_occ_otu == 1
  )
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
    }
  }



  // ----------------------------------------------------------
  // Capture OTU random effects
  // ----------------------------------------------------------


  if (
    use_cap_otu == 1
  )
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



  // ----------------------------------------------------------
  // Abundance OTU random effects
  // ----------------------------------------------------------


  if (
    use_abund_otu == 1
  )
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



  // ----------------------------------------------------------
  // Biological-sample abundance random effects
  // ----------------------------------------------------------


  if (
    use_sample_re == 1
  )
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



  // ----------------------------------------------------------
  // Sample x OTU abundance random effects
  // ----------------------------------------------------------


  if (
    use_sample_otu_re == 1
  )
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
  // OTU-specific zero-inflation random effects
  // ----------------------------------------------------------
  //
  // These are included only for ZIP/ZINB.
  //
  //     b_zi_m ~ N(0, sigma_zi^2)
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
    }
  }



  // ==========================================================
  // 10. OCCUPANCY LINEAR PREDICTOR
  // ==========================================================


  vector<Type> eta_occ =
    X_occ *
    beta_occ;


  if (
    use_occ_otu == 1
  )
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
    }
  }


  // Model:
  //
  //     logit(psi_im)
  //
  //       =
  //
  //     X_occ_im beta_occ
  //
  //       +
  //
  //     b_occ_m



  // ==========================================================
  // 11. CAPTURE LINEAR PREDICTOR
  // ==========================================================


  vector<Type> eta_cap =
    X_cap *
    beta_cap;


  if (
    use_cap_otu == 1
  )
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
    }
  }


  // Model:
  //
  //     logit(p_ijm)
  //
  //       =
  //
  //     X_cap_ijm beta_cap
  //
  //       +
  //
  //     b_cap_m



  // ==========================================================
  // 12. ABUNDANCE LINEAR PREDICTOR
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
    // Sequencing-depth or other abundance offset
    // --------------------------------------------------------

    eta_abund(r) +=
      offset_abund(r);


    // --------------------------------------------------------
    // OTU abundance effect
    // --------------------------------------------------------

    if (
      use_abund_otu == 1
    )
    {
      eta_abund(r) +=
        b_abund_otu(
          row_otu(r)
        );
    }


    // --------------------------------------------------------
    // Biological-sample abundance effect
    // --------------------------------------------------------

    if (
      use_sample_re == 1
    )
    {
      eta_abund(r) +=
        b_sample(
          row_sample_id(r)
        );
    }


    // --------------------------------------------------------
    // Sample x OTU abundance effect
    // --------------------------------------------------------

    if (
      use_sample_otu_re == 1
    )
    {
      eta_abund(r) +=
        b_sample_otu(
          row_sample_otu_re(r)
        );
    }
  }



  // ==========================================================
  // 13. ZERO-INFLATION LINEAR PREDICTOR BY OTU
  // ==========================================================
  //
  // This vector is mainly useful for reporting.
  //
  // For OTU m:
  //
  //     eta_zi_m
  //
  //       =
  //
  //     zi_intercept
  //
  //       +
  //
  //     b_zi_otu(m)
  //
  // and:
  //
  //     pi_m = invlogit(eta_zi_m)
  //
  // When random_zi_otu is disabled:
  //
  //     eta_zi_m = zi_intercept
  //
  // for every OTU.
  // ==========================================================


  vector<Type> eta_zi_otu(
    b_zi_otu.size()
  );


  vector<Type> pi_zi_otu(
    b_zi_otu.size()
  );


  for (
    int m = 0;
    m < b_zi_otu.size();
    m++
  )
  {

    eta_zi_otu(m) =
      zi_intercept;


    if (
      use_zi_otu == 1 &&
      (
        family_code == 2 ||
        family_code == 3
      )
    )
    {
      eta_zi_otu(m) +=
        b_zi_otu(m);
    }


    // Probability is calculated here for reporting.
    //
    // The likelihood below does NOT use log(pi_zi_otu)
    // directly. Instead it calculates log(pi) from eta_zi
    // using the numerically stable helper functions.

    pi_zi_otu(m) =
      Type(1) /
      (
        Type(1) +
        exp(
          -eta_zi_otu(m)
        )
      );
  }



  // ==========================================================
  // 14. PCR / READ-LEVEL COUNT LOG-LIKELIHOOD
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
    // Expected count
    // --------------------------------------------------------
    //
    // Log link:
    //
    //     lambda_ijkm
    //       =
    //       exp(eta_abund_ijkm)
    // --------------------------------------------------------


    Type mu =
      exp(
        eta_abund(r)
      );


    Type lp =
      Type(0);



    // --------------------------------------------------------
    // Zero-inflation linear predictor for this observation
    // --------------------------------------------------------
    //
    // All observations belonging to OTU m share the same
    // zero-inflation random intercept b_zi_m.
    // --------------------------------------------------------


    Type eta_zi =
      zi_intercept;


    if (
      use_zi_otu == 1 &&
      (
        family_code == 2 ||
        family_code == 3
      )
    {
      eta_zi +=
        b_zi_otu(
          row_otu(r)
        );
    }



    // --------------------------------------------------------
    // Numerically stable log probabilities
    // --------------------------------------------------------
    //
    //     log_pi
    //       = log(pi_m)
    //
    //     log_1m_pi
    //       = log(1 - pi_m)
    //
    // These are preferable to:
    //
    //     log(pi)
    //     log(1-pi)
    //
    // after explicitly constructing pi.
    // --------------------------------------------------------


    Type log_pi =
      log_invlogit(
        eta_zi
      );


    Type log_1m_pi =
      log1m_invlogit(
        eta_zi
      );



    // ========================================================
    // POISSON
    // ========================================================
    //
    //     Y ~ Poisson(mu)
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
    }



    // ========================================================
    // NEGATIVE BINOMIAL 2
    // ========================================================
    //
    //     E(Y) = mu
    //
    //     Var(Y)
    //       =
    //       mu + mu^2/theta
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
    }



    // ========================================================
    // ZERO-INFLATED POISSON
    // ========================================================
    //
    // For Y = 0:
    //
    //     P(Y=0)
    //
    //       =
    //
    //     pi_m
    //
    //       +
    //
    //     (1-pi_m) P_Pois(Y=0 | mu)
    //
    //
    // For Y > 0:
    //
    //     P(Y=y)
    //
    //       =
    //
    //     (1-pi_m) P_Pois(Y=y | mu)
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


      if (
        y(r) == Type(0)
      )
      {

        lp =
          logspace_add(

            log_pi,

            log_1m_pi +
            log_count_component

          );
      }

      else
      {

        lp =
          log_1m_pi +
          log_count_component;
      }
    }



    // ========================================================
    // ZERO-INFLATED NEGATIVE BINOMIAL 2
    // ========================================================
    //
    // For Y = 0:
    //
    //     P(Y=0)
    //
    //       =
    //
    //     pi_m
    //
    //       +
    //
    //     (1-pi_m)
    //     P_NB(Y=0 | mu, theta)
    //
    //
    // For Y > 0:
    //
    //     P(Y=y)
    //
    //       =
    //
    //     (1-pi_m)
    //     P_NB(Y=y | mu, theta)
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
        y(r) == Type(0)
      )
      {

        lp =
          logspace_add(

            log_pi,

            log_1m_pi +
            log_count_component

          );
      }

      else
      {

        lp =
          log_1m_pi +
          log_count_component;
      }
    }



    // Store PCR/read-level conditional log likelihood.

    log_count(r) =
      lp;
  }



  // ==========================================================
  // 15. COMBINE PCR REPLICATES WITHIN SAMPLE x OTU
  // ==========================================================
  //
  // Conditional on successful capture:
  //
  //     A_ijm = 1
  //
  // the PCR/read observations are conditionally independent.
  //
  // Therefore:
  //
  //     L_count_ijm
  //
  //       =
  //
  //     product_k
  //     f(y_ijkm | A_ijm = 1)
  //
  //
  // On the log scale:
  //
  //     log L_count_ijm
  //
  //       =
  //
  //     sum_k log f(y_ijkm | A_ijm = 1)
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
  // 16. ANALYTICALLY MARGINALIZE CAPTURE STATE A
  // ==========================================================
  //
  // Capture lives at:
  //
  //     biological sample x OTU
  //
  // NOT at the PCR/read level.
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
    // CASE 1:
    // At least one positive PCR replicate
    // --------------------------------------------------------
    //
    // If any y_ijkm > 0:
    //
    //     A_ijm = 1
    //
    // necessarily.
    //
    // Therefore:
    //
    //     L_ijm
    //
    //       =
    //
    //     p_ijm
    //
    //       *
    //
    //     product_k f(y_ijkm)
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
    // CASE 2:
    // All PCR replicates are zero
    // --------------------------------------------------------
    //
    // Two latent explanations are possible.
    //
    // A = 0:
    //
    //     capture failed
    //
    //     probability = 1-p
    //
    // OR
    //
    // A = 1:
    //
    //     capture succeeded
    //
    //     but all PCR/read observations happened to be zero.
    //
    //
    // Therefore:
    //
    //     L_ijm
    //
    //       =
    //
    //     (1-p_ijm)
    //
    //       +
    //
    //     p_ijm
    //     product_k f(0)
    //
    //
    // For ZIP/ZINB, f(0) already includes the structural-zero
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
  // 17. COMBINE BIOLOGICAL SAMPLES WITHIN SITE x OTU
  // ==========================================================
  //
  // Conditional on:
  //
  //     Z_im = 1
  //
  // the sample likelihood is:
  //
  //     L(data_im | Z_im=1)
  //
  //       =
  //
  //     product_j L_ijm
  //
  //
  // On the log scale:
  //
  //     log L(data_im | Z_im=1)
  //
  //       =
  //
  //     sum_j log L_ijm
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
  // 18. ANALYTICALLY MARGINALIZE OCCUPANCY STATE Z
  // ==========================================================
  //
  // Occupancy lives at:
  //
  //     SITE x OTU
  //
  // This is the critical hierarchical distinction.
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

    // Occupancy probability for this site x OTU.

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
    // CASE 1:
    // At least one positive observation at site x OTU
    // --------------------------------------------------------
    //
    // If anything is positive:
    //
    //     Z_im = 1
    //
    // necessarily.
    //
    // Therefore:
    //
    //     L_im
    //
    //       =
    //
    //     psi_im
    //
    //       *
    //
    //     L(data_im | Z_im=1)
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
    // CASE 2:
    // Entire site x OTU history is zero
    // --------------------------------------------------------
    //
    // Two possibilities:
    //
    // Z = 0:
    //
    //     species absent
    //
    //     probability = 1-psi
    //
    // OR
    //
    // Z = 1:
    //
    //     species present
    //
    //     but all samples/PCR observations are zero.
    //
    //
    // Therefore:
    //
    //     L_im
    //
    //       =
    //
    //     (1-psi_im)
    //
    //       +
    //
    //     psi_im
    //     L(all-zero data | Z_im=1)
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


    // Add this SITE x OTU contribution to total
    // negative log likelihood.

    nll -=
      log_site(g);
  }



  // ==========================================================
  // 19. DERIVED CAPTURE PROBABILITIES
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
  }



  // ==========================================================
  // 20. DERIVED ABUNDANCE
  // ==========================================================


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
  }



  // ==========================================================
  // 21. REPORT MODEL QUANTITIES TO R
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
    eta_zi_otu
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


  REPORT(
    log_sample
  );


  REPORT(
    log_site
  );



  // ==========================================================
  // 22. ADREPORT
  // ==========================================================
  //
  // ADREPORT allows TMB/sdreport to calculate approximate
  // standard errors for transformed/derived quantities using
  // automatic differentiation and the delta method.
  // ==========================================================


  ADREPORT(
    psi
  );


  ADREPORT(
    capture_prob
  );


  // Only scientifically meaningful for ZIP/ZINB.
  //
  // The R wrapper maps the ZI parameters out for Poisson/NB.

  if (
    family_code == 2 ||
    family_code == 3
  )
  {
    ADREPORT(
      pi_zi_otu
    );
  }



  // ==========================================================
  // 23. RETURN NEGATIVE LOG-LIKELIHOOD
  // ==========================================================


  return nll;
}
