#include <TMB.hpp>


// ============================================================
// Numerically stable logistic helper functions
// ============================================================

// log(invlogit(x))
template<class Type>
Type log_invlogit(Type x)
{
  return -logspace_add(Type(0), -x);
}


// log(1 - invlogit(x))
template<class Type>
Type log1m_invlogit(Type x)
{
  return -logspace_add(Type(0), x);
}


// ============================================================
// Joint hierarchical eDNA model
// ============================================================
//
// Hierarchy:
//
// SITE x OTU:
//
//   Z_im ~ Bernoulli(psi_im)
//
// SAMPLE x OTU:
//
//   A_ijm | Z_im = 1 ~ Bernoulli(p_ijm)
//
// PCR replicate:
//
//   Y_ijkm | A_ijm = 1 ~ Count(lambda_ijkm)
//
// where:
//
//   i = site
//   j = biological sample
//   k = PCR replicate
//   m = OTU
//
// Occupancy:
//
//   logit(psi_im)
//     = X_occ beta_occ
//       + b_occ_m
//
// Capture:
//
//   logit(p_ijm)
//     = X_cap beta_cap
//       + b_cap_m
//
// Abundance:
//
//   log(lambda_ijkm)
//     = X_abund beta_abund
//       + offset
//       + b_abund_m
//       + b_sample_j
//       + b_sampleOTU_jm
//
// Zero inflation:
//
//   logit(pi_m)
//     = zi_intercept
//       + b_zi_m
//
// The latent capture state A is marginalized analytically.
//
// The latent occupancy state Z is marginalized analytically.
//
// Gaussian random effects are integrated by TMB using the
// Laplace approximation when supplied to MakeADFun(random=...).
// ============================================================


template<class Type>
Type objective_function<Type>::operator() ()
{

  // ==========================================================
  // 1. DATA
  // ==========================================================

  DATA_VECTOR(y);

  // y:
  // PCR/read-level observed counts.


  // ----------------------------------------------------------
  // Design matrices
  // ----------------------------------------------------------

  DATA_MATRIX(X_occ);

  // One row per site x OTU.


  DATA_MATRIX(X_cap);

  // One row per sample x OTU.


  DATA_MATRIX(X_abund);

  // One row per PCR/read observation.


  // ----------------------------------------------------------
  // Abundance offset
  // ----------------------------------------------------------

  DATA_VECTOR(offset_abund);


  // ----------------------------------------------------------
  // Hierarchical mapping indices
  // ----------------------------------------------------------

  DATA_IVECTOR(occ_otu);

  // site x OTU group -> OTU


  DATA_IVECTOR(sample_otu);

  // sample x OTU group -> OTU


  DATA_IVECTOR(row_sample_group);

  // PCR/read row -> sample x OTU group


  DATA_IVECTOR(sample_site_group);

  // sample x OTU group -> site x OTU group


  DATA_IVECTOR(row_otu);

  // PCR/read row -> OTU


  DATA_IVECTOR(row_sample_id);

  // PCR/read row -> biological sample


  DATA_IVECTOR(row_sample_otu_re);

  // PCR/read row -> sample x OTU random-effect level


  // ----------------------------------------------------------
  // Positive-history indicators
  // ----------------------------------------------------------

  DATA_IVECTOR(sample_positive);

  // sample_positive(s) = 1 if at least one PCR replicate
  // within sample x OTU group s is positive.


  DATA_IVECTOR(site_positive);

  // site_positive(g) = 1 if at least one observation
  // within site x OTU group g is positive.


  // ----------------------------------------------------------
  // Number of groups
  // ----------------------------------------------------------

  DATA_INTEGER(n_site_groups);

  DATA_INTEGER(n_sample_groups);


  // ----------------------------------------------------------
  // Distribution
  // ----------------------------------------------------------

  DATA_INTEGER(family_code);

  // 0 = Poisson
  // 1 = NB2
  // 2 = ZIP
  // 3 = ZINB


  // ==========================================================
  // 2. RANDOM-EFFECT SWITCHES
  // ==========================================================

  DATA_INTEGER(use_occ_otu);

  DATA_INTEGER(use_cap_otu);

  DATA_INTEGER(use_abund_otu);

  DATA_INTEGER(use_sample_re);

  DATA_INTEGER(use_sample_otu_re);

  DATA_INTEGER(use_zi_otu);


  // ==========================================================
  // 3. VALIDATE FAMILY
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
  // 4. FIXED EFFECTS
  // ==========================================================

  PARAMETER_VECTOR(beta_occ);

  PARAMETER_VECTOR(beta_cap);

  PARAMETER_VECTOR(beta_abund);


  // ==========================================================
  // 5. RANDOM EFFECTS
  // ==========================================================

  PARAMETER_VECTOR(b_occ_otu);

  // OTU-specific occupancy random effects.


  PARAMETER_VECTOR(b_cap_otu);

  // OTU-specific capture random effects.


  PARAMETER_VECTOR(b_abund_otu);

  // OTU-specific abundance random effects.


  PARAMETER_VECTOR(b_sample);

  // Biological-sample abundance random effects.


  PARAMETER_VECTOR(b_sample_otu);

  // Sample x OTU abundance random effects.


  PARAMETER_VECTOR(b_zi_otu);

  // OTU-specific zero-inflation random effects.


  // ==========================================================
  // 6. RANDOM-EFFECT STANDARD DEVIATIONS
  // ==========================================================

  PARAMETER(log_sd_occ_otu);

  PARAMETER(log_sd_cap_otu);

  PARAMETER(log_sd_abund_otu);

  PARAMETER(log_sd_sample);

  PARAMETER(log_sd_sample_otu);

  PARAMETER(log_sd_zi_otu);


  Type sd_occ_otu =
    exp(log_sd_occ_otu);


  Type sd_cap_otu =
    exp(log_sd_cap_otu);


  Type sd_abund_otu =
    exp(log_sd_abund_otu);


  Type sd_sample =
    exp(log_sd_sample);


  Type sd_sample_otu =
    exp(log_sd_sample_otu);


  Type sd_zi_otu =
    exp(log_sd_zi_otu);


  // ==========================================================
  // 7. NB DISPERSION
  // ==========================================================

  PARAMETER(log_theta);


  Type theta =
    exp(log_theta);

  // NB2 parameterization:
  //
  // E(Y) = mu
  //
  // Var(Y) = mu + mu^2 / theta
  //
  // Currently theta is shared across OTUs.


  // ==========================================================
  // 8. ZERO-INFLATION INTERCEPT
  // ==========================================================

  PARAMETER(zi_intercept);

  // For ZIP/ZINB:
  //
  // logit(pi_m)
  //   =
  // zi_intercept + b_zi_otu(m)


  // ==========================================================
  // 9. NEGATIVE LOG-LIKELIHOOD
  // ==========================================================

  Type nll =
    Type(0);


  // ==========================================================
  // 10. RANDOM-EFFECT DISTRIBUTIONS
  // ==========================================================
  //
  // All random-effect blocks are independent.
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
      nll -=
        dnorm(
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
      nll -=
        dnorm(
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
      nll -=
        dnorm(
          b_abund_otu(m),
          Type(0),
          sd_abund_otu,
          true
        );
    }
  }


  // ----------------------------------------------------------
  // Sample random effects
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
      nll -=
        dnorm(
          b_sample(j),
          Type(0),
          sd_sample,
          true
        );
    }
  }


  // ----------------------------------------------------------
  // Sample x OTU random effects
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
      nll -=
        dnorm(
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
  // Only active for ZIP/ZINB.
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
      nll -=
        dnorm(
          b_zi_otu(m),
          Type(0),
          sd_zi_otu,
          true
        );
    }
  }


  // ==========================================================
  // 11. OCCUPANCY LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_occ =
    X_occ * beta_occ;


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


  // ==========================================================
  // 12. CAPTURE LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_cap =
    X_cap * beta_cap;


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


  // ==========================================================
  // 13. ABUNDANCE LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_abund =
    X_abund * beta_abund;


  for (
    int r = 0;
    r < y.size();
    r++
  )
  {

    // Add abundance offset.

    eta_abund(r) +=
      offset_abund(r);


    // OTU abundance effect.

    if (
      use_abund_otu == 1
    )
    {
      eta_abund(r) +=
        b_abund_otu(
          row_otu(r)
        );
    }


    // Biological-sample abundance effect.

    if (
      use_sample_re == 1
    )
    {
      eta_abund(r) +=
        b_sample(
          row_sample_id(r)
        );
    }


    // Sample x OTU abundance effect.

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
  // 14. ZERO-INFLATION PREDICTOR BY OTU
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


    pi_zi_otu(m) =
      exp(
        log_invlogit(
          eta_zi_otu(m)
        )
      );
  }


  // ==========================================================
  // 15. PCR/READ-LEVEL COUNT LIKELIHOOD
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


    Type lp =
      Type(0);


    // --------------------------------------------------------
    // Zero-inflation predictor
    // --------------------------------------------------------

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
        b_zi_otu(
          row_otu(r)
        );
    }


    // --------------------------------------------------------
    // Stable log zero-inflation probabilities
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
    // Poisson
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
    // Negative Binomial 2
    // ========================================================

    else if (
      family_code == 1
    )
    {
      Type variance =
        mu +
        mu * mu / theta;


      lp =
        dnbinom2(
          y(r),
          mu,
          variance,
          true
        );
    }


    // ========================================================
    // Zero-Inflated Poisson
    // ========================================================
    //
    // If y = 0:
    //
    // P(Y=0)
    //   =
    // pi
    // +
    // (1-pi) Pois(0 | mu)
    //
    // If y > 0:
    //
    // P(Y=y)
    //   =
    // (1-pi) Pois(y | mu)
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
    // Zero-Inflated Negative Binomial 2
    // ========================================================
    //
    // If y = 0:
    //
    // P(Y=0)
    //   =
    // pi
    // +
    // (1-pi) NB(0 | mu, theta)
    //
    // If y > 0:
    //
    // P(Y=y)
    //   =
    // (1-pi) NB(y | mu, theta)
    // ========================================================

    else if (
      family_code == 3
    )
    {
      Type variance =
        mu +
        mu * mu / theta;


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


    log_count(r) =
      lp;
  }


  // ==========================================================
  // 16. COMBINE PCR REPLICATES WITHIN SAMPLE x OTU
  // ==========================================================
  //
  // Conditional on A_ijm = 1:
  //
  // L_count_ijm
  //   =
  // product_k f(y_ijkm)
  //
  // Therefore:
  //
  // log L_count_ijm
  //   =
  // sum_k log f(y_ijkm)
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
  // 17. MARGINALIZE SAMPLE-LEVEL CAPTURE A
  // ==========================================================
  //
  // A_ijm exists once per biological sample x OTU.
  //
  // If any PCR replicate is positive:
  //
  //   A_ijm = 1
  //
  // and:
  //
  //   L_ijm
  //     =
  //   p_ijm * product_k f(y_ijkm)
  //
  //
  // If all PCR replicates are zero:
  //
  //   L_ijm
  //     =
  //   (1-p_ijm)
  //     +
  //   p_ijm * product_k f(0)
  //
  // This sums over A_ijm = 0 and A_ijm = 1.
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
    // At least one positive PCR replicate
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
    // All PCR replicates are zero
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
  // 18. COMBINE SAMPLES WITHIN SITE x OTU
  // ==========================================================
  //
  // Conditional on occupancy Z_im = 1:
  //
  // L(data_im | Z_im = 1)
  //
  //   =
  //
  // product_j L_ijm
  //
  // therefore:
  //
  // log L(data_im | Z_im = 1)
  //
  //   =
  //
  // sum_j log L_ijm
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
  // 19. MARGINALIZE SITE-LEVEL OCCUPANCY Z
  // ==========================================================
  //
  // Z_im exists ONCE per site x OTU.
  //
  // This is the important distinction:
  //
  // occupancy is NOT at the PCR/read observation level.
  //
  //
  // If site x OTU has any positive observation:
  //
  //   Z_im = 1
  //
  // therefore:
  //
  //   L_im
  //     =
  //   psi_im *
  //   L(data_im | Z_im=1)
  //
  //
  // If the entire site x OTU history is zero:
  //
  //   L_im
  //     =
  //   (1-psi_im)
  //     +
  //   psi_im *
  //   L(all-zero data | Z_im=1)
  //
  // This analytically sums over Z_im.
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

    // Occupancy probability on natural scale.

    psi(g) =
      exp(
        log_invlogit(
          eta_occ(g)
        )
      );


    // Stable log occupancy probability.

    Type log_psi =
      log_invlogit(
        eta_occ(g)
      );


    // Stable log non-occupancy probability.

    Type log_1mpsi =
      log1m_invlogit(
        eta_occ(g)
      );


    // --------------------------------------------------------
    // Positive site x OTU history
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

    else
    {
      log_site(g) =
        logspace_add(
          log_1mpsi,
          log_psi +
            log_site_conditional(g)
        );
    }


    // Add contribution to negative log-likelihood.

    nll -=
      log_site(g);
  }


  // ==========================================================
  // 20. DERIVED CAPTURE PROBABILITIES
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
      exp(
        log_invlogit(
          eta_cap(s)
        )
      );
  }


  // ==========================================================
  // 21. DERIVED ABUNDANCE
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
  // 22. REPORT
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
    log_sample_count
  );


  REPORT(
    log_sample
  );


  REPORT(
    log_site_conditional
  );


  REPORT(
    log_site
  );


  // ==========================================================
  // 23. ADREPORT
  // ==========================================================
  //
  // ADREPORT allows sdreport() to obtain approximate SEs
  // for these transformed parameters using the delta method.
  // ==========================================================

  ADREPORT(
    psi
  );


  ADREPORT(
    capture_prob
  );


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
  // 24. RETURN NEGATIVE LOG-LIKELIHOOD
  // ==========================================================

  return nll;
}
