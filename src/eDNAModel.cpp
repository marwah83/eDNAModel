#include <TMB.hpp>


// ------------------------------------------------------------
// log(invlogit(x))
// ------------------------------------------------------------

template<class Type>
Type log_invlogit(Type x)
{
  return -logspace_add(Type(0), -x);
}


// ------------------------------------------------------------
// log(1 - invlogit(x))
// ------------------------------------------------------------

template<class Type>
Type log1m_invlogit(Type x)
{
  return -logspace_add(Type(0), x);
}


// ============================================================
// Joint eDNA likelihood
// ============================================================

template<class Type>
Type objective_function<Type>::operator() ()
{

  // ==========================================================
  // DATA
  // ==========================================================

  DATA_VECTOR(y);

  // Fixed-effect design matrices
  DATA_MATRIX(X_occ);
  DATA_MATRIX(X_cap);
  DATA_MATRIX(X_abund);

  // Offset for abundance
  DATA_VECTOR(offset_abund);

  // Taxon index for each site-taxon occupancy row
  DATA_IVECTOR(occ_otu);

  // Taxon index for each sample-taxon group
  DATA_IVECTOR(sample_otu);

  // Biological sample index
  DATA_IVECTOR(sample_id);

  // Sample x OTU interaction index
  DATA_IVECTOR(sample_otu_re);

  // Mapping PCR rows -> sample-taxon group
  DATA_IVECTOR(row_sample_group);

  // Mapping sample-taxon groups -> site-taxon group
  DATA_IVECTOR(sample_site_group);

  // Mapping PCR rows -> OTU
  DATA_IVECTOR(row_otu);

  // Mapping PCR rows -> biological sample
  DATA_IVECTOR(row_sample_id);

  // Mapping PCR rows -> sample x OTU random effect
  DATA_IVECTOR(row_sample_otu_re);

  // Does each sample-taxon group have any positive read?
  DATA_IVECTOR(sample_positive);

  // Does each site-taxon group have any positive read?
  DATA_IVECTOR(site_positive);

  // Number of site x OTU combinations
  DATA_INTEGER(n_site_groups);

  // Number of sample x OTU combinations
  DATA_INTEGER(n_sample_groups);

  // Count distribution:
  // 0 = Poisson
  // 1 = NB2
  // 2 = ZIP
  // 3 = ZINB
  DATA_INTEGER(family_code);

  // Random-effect switches
  DATA_INTEGER(use_occ_otu);
  DATA_INTEGER(use_cap_otu);
  DATA_INTEGER(use_abund_otu);
  DATA_INTEGER(use_sample_re);
  DATA_INTEGER(use_sample_otu_re);


  // ==========================================================
  // FIXED-EFFECT PARAMETERS
  // ==========================================================

  PARAMETER_VECTOR(beta_occ);

  PARAMETER_VECTOR(beta_cap);

  PARAMETER_VECTOR(beta_abund);


  // ==========================================================
  // RANDOM EFFECTS
  // ==========================================================

  PARAMETER_VECTOR(b_occ_otu);
  PARAMETER_VECTOR(b_cap_otu);
  PARAMETER_VECTOR(b_abund_otu);

  PARAMETER_VECTOR(b_sample);
  PARAMETER_VECTOR(b_sample_otu);


  // ==========================================================
  // RANDOM-EFFECT STANDARD DEVIATIONS
  // ==========================================================

  PARAMETER(log_sd_occ_otu);
  PARAMETER(log_sd_cap_otu);
  PARAMETER(log_sd_abund_otu);

  PARAMETER(log_sd_sample);
  PARAMETER(log_sd_sample_otu);


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


  // ==========================================================
  // ABUNDANCE PARAMETERS
  // ==========================================================

  // NB2 dispersion
  PARAMETER(log_theta);

  Type theta =
    exp(log_theta);

  // Zero inflation
  PARAMETER(zi_intercept);

  Type pi_zi =
    Type(1) /
    (
      Type(1) +
      exp(-zi_intercept)
    );


  // ==========================================================
  // NEGATIVE LOG-LIKELIHOOD
  // ==========================================================

  Type nll = 0;


  // ==========================================================
  // RANDOM-EFFECT PRIORS
  // ==========================================================

  if (use_occ_otu == 1)
  {
    for (int m = 0; m < b_occ_otu.size(); m++)
    {
      nll -= dnorm(
        b_occ_otu(m),
        Type(0),
        sd_occ_otu,
        true
      );
    }
  }


  if (use_cap_otu == 1)
  {
    for (int m = 0; m < b_cap_otu.size(); m++)
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
    for (int m = 0; m < b_abund_otu.size(); m++)
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
    for (int s = 0; s < b_sample.size(); s++)
    {
      nll -= dnorm(
        b_sample(s),
        Type(0),
        sd_sample,
        true
      );
    }
  }


  if (use_sample_otu_re == 1)
  {
    for (int q = 0; q < b_sample_otu.size(); q++)
    {
      nll -= dnorm(
        b_sample_otu(q),
        Type(0),
        sd_sample_otu,
        true
      );
    }
  }


  // ==========================================================
  // OCCUPANCY LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_occ =
    X_occ * beta_occ;

  for (int g = 0; g < n_site_groups; g++)
  {
    if (use_occ_otu == 1)
    {
      eta_occ(g) +=
        b_occ_otu(
          occ_otu(g)
        );
    }
  }


  // ==========================================================
  // CAPTURE LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_cap =
    X_cap * beta_cap;

  for (int s = 0; s < n_sample_groups; s++)
  {
    if (use_cap_otu == 1)
    {
      eta_cap(s) +=
        b_cap_otu(
          sample_otu(s)
        );
    }
  }


  // ==========================================================
  // ABUNDANCE LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_abund =
    X_abund * beta_abund;

  for (int r = 0; r < y.size(); r++)
  {

    eta_abund(r) +=
      offset_abund(r);

    if (use_abund_otu == 1)
    {
      eta_abund(r) +=
        b_abund_otu(
          row_otu(r)
        );
    }

    if (use_sample_re == 1)
    {
      eta_abund(r) +=
        b_sample(
          row_sample_id(r)
        );
    }

    if (use_sample_otu_re == 1)
    {
      eta_abund(r) +=
        b_sample_otu(
          row_sample_otu_re(r)
        );
    }
  }


  // ==========================================================
  // COUNT LOG LIKELIHOOD
  // ==========================================================

  vector<Type> log_count(
    y.size()
  );

  for (int r = 0; r < y.size(); r++)
  {

    Type mu =
      exp(
        eta_abund(r)
      );

    Type lp;


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
    }


    // --------------------------------------------------------
    // NB2
    //
    // Var(Y) = mu + mu^2 / theta
    // --------------------------------------------------------

    else if (family_code == 1)
    {

      lp =
        dnbinom2(
          y(r),
          mu,
          mu +
            mu * mu /
            theta,
          true
        );
    }


    // --------------------------------------------------------
    // ZIP
    // --------------------------------------------------------

    else if (family_code == 2)
    {

      Type log_count_component =
        dpois(
          y(r),
          mu,
          true
        );

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
    // ZINB
    // --------------------------------------------------------

    else
    {

      Type log_count_component =
        dnbinom2(
          y(r),
          mu,
          mu +
            mu * mu /
            theta,
          true
        );

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


    log_count(r) =
      lp;
  }


  // ==========================================================
  // SAMPLE-LEVEL LIKELIHOOD
  //
  // A is analytically marginalized here
  // ==========================================================

  vector<Type> log_sample(
    n_sample_groups
  );

  log_sample.setZero();


  // Sum replicate count likelihoods within sample x OTU
  for (int r = 0; r < y.size(); r++)
  {

    int s =
      row_sample_group(r);

    log_sample(s) +=
      log_count(r);
  }


  for (int s = 0; s < n_sample_groups; s++)
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
    //
    // A must equal 1
    // --------------------------------------------------------

    if (sample_positive(s) == 1)
    {

      log_sample(s) =
        log_p +
        log_sample(s);
    }


    // --------------------------------------------------------
    // All-zero biological sample
    //
    // Marginalize:
    //
    // (1-p) + p P(Y=0 | A=1)
    // --------------------------------------------------------

    else
    {

      log_sample(s) =
        logspace_add(
          log_1mp,

          log_p +
          log_sample(s)
        );
    }
  }


  // ==========================================================
  // SITE x OTU LIKELIHOOD
  //
  // Z is analytically marginalized here
  // ==========================================================

  vector<Type> log_site(
    n_site_groups
  );

  log_site.setZero();


  for (int s = 0; s < n_sample_groups; s++)
  {

    int g =
      sample_site_group(s);

    log_site(g) +=
      log_sample(s);
  }


  // ==========================================================
  // Marginalize occupancy Z
  // ==========================================================

  vector<Type> psi(
    n_site_groups
  );

  for (int g = 0; g < n_site_groups; g++)
  {

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
    // Positive site x OTU
    //
    // Z must equal 1
    // --------------------------------------------------------

    if (site_positive(g) == 1)
    {

      log_site(g) =
        log_psi +
        log_site(g);
    }


    // --------------------------------------------------------
    // Entire site x OTU history is zero
    //
    // Marginalize:
    //
    // (1-psi) + psi P(all data = 0 | Z=1)
    // --------------------------------------------------------

    else
    {

      log_site(g) =
        logspace_add(
          log_1mpsi,

          log_psi +
          log_site(g)
        );
    }


    nll -=
      log_site(g);
  }


  // ==========================================================
  // Derived parameters
  // ==========================================================

  vector<Type> capture_prob(
    n_sample_groups
  );

  for (int s = 0; s < n_sample_groups; s++)
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


  vector<Type> lambda(
    y.size()
  );

  for (int r = 0; r < y.size(); r++)
  {

    lambda(r) =
      exp(
        eta_abund(r)
      );
  }


  // ==========================================================
  // REPORT
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
    log_site
  );


  // Formal delta-method SEs from joint model
  ADREPORT(
    psi
  );

  ADREPORT(
    capture_prob
  );


  return nll;
}
