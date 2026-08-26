#include <TMB.hpp>

// ============================================================
// Numerically stable logistic helpers
// TMB already supplies logspace_add()
// ============================================================

template<class Type>
Type log_invlogit(Type x)
{
  return -logspace_add(Type(0), -x);
}

template<class Type>
Type log1m_invlogit(Type x)
{
  return -logspace_add(Type(0), x);
}


// ============================================================
// Joint observed-data eDNA likelihood
//
// Hierarchy:
//
// Z_im ~ Bernoulli(psi_im)
//
// A_ijm | Z_im = 1 ~ Bernoulli(p_ijm)
//
// Y_ijkm | A_ijm = 1 ~ Count(lambda_ijkm)
//
// Z and A are analytically marginalized.
// ============================================================

template<class Type>
Type objective_function<Type>::operator() ()
{

  // ==========================================================
  // DATA
  // ==========================================================

  DATA_VECTOR(y);

  DATA_MATRIX(X_occ);
  DATA_MATRIX(X_cap);
  DATA_MATRIX(X_abund);

  DATA_VECTOR(offset_abund);

  // Site x OTU -> OTU
  DATA_IVECTOR(occ_otu);

  // Sample x OTU -> OTU
  DATA_IVECTOR(sample_otu);

  // PCR/read row -> sample x OTU
  DATA_IVECTOR(row_sample_group);

  // Sample x OTU -> site x OTU
  DATA_IVECTOR(sample_site_group);

  // PCR/read row -> OTU
  DATA_IVECTOR(row_otu);

  // PCR/read row -> biological sample
  DATA_IVECTOR(row_sample_id);

  // PCR/read row -> sample x OTU RE
  DATA_IVECTOR(row_sample_otu_re);

  // Whether sample x OTU has >= 1 positive count
  DATA_IVECTOR(sample_positive);

  // Whether site x OTU has >= 1 positive count
  DATA_IVECTOR(site_positive);

  DATA_INTEGER(n_site_groups);
  DATA_INTEGER(n_sample_groups);

  // ----------------------------------------------------------
  // Count family
  //
  // 0 = Poisson
  // 1 = NB2
  // 2 = ZIP
  // 3 = ZINB
  // ----------------------------------------------------------

  DATA_INTEGER(family_code);

  // ----------------------------------------------------------
  // Random-effect switches
  // ----------------------------------------------------------

  DATA_INTEGER(use_occ_otu);
  DATA_INTEGER(use_cap_otu);
  DATA_INTEGER(use_abund_otu);
  DATA_INTEGER(use_sample_re);
  DATA_INTEGER(use_sample_otu_re);


  // ==========================================================
  // FIXED EFFECTS
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

  // ZIP / ZINB structural-zero probability
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

  Type nll = Type(0);


  // ==========================================================
  // RANDOM-EFFECT GAUSSIAN DENSITIES
  // ==========================================================

  if (use_occ_otu == 1)
  {

    for (int j = 0;
         j < b_occ_otu.size();
         j++)
    {

      nll -= dnorm(
        b_occ_otu(j),
        Type(0),
        sd_occ_otu,
        true
      );
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
    }
  }


  // ==========================================================
  // OCCUPANCY LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_occ =
    X_occ * beta_occ;

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
    }
  }


  // ==========================================================
  // CAPTURE LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_cap =
    X_cap * beta_cap;

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
    }
  }


  // ==========================================================
  // ABUNDANCE LINEAR PREDICTOR
  // ==========================================================

  vector<Type> eta_abund =
    X_abund * beta_abund;

  for (int r = 0;
       r < y.size();
       r++)
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
  // READ-LEVEL COUNT LOG LIKELIHOOD
  // ==========================================================

  vector<Type> log_count(
    y.size()
  );


  for (int r = 0;
       r < y.size();
       r++)
  {

    Type mu =
      exp(
        eta_abund(r)
      );


    Type lp =
      Type(0);


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
    // Negative Binomial 2
    //
    // Var(Y) =
    // mu + mu^2 / theta
    // --------------------------------------------------------

    else if (family_code == 1)
    {

      Type variance =
        mu +
        mu * mu /
        theta;

      lp =
        dnbinom2(
          y(r),
          mu,
          variance,
          true
        );
    }


    // --------------------------------------------------------
    // Zero-inflated Poisson
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
    // Zero-inflated NB2
    // --------------------------------------------------------

    else if (family_code == 3)
    {

      Type variance =
        mu +
        mu * mu /
        theta;


      Type log_count_component =
        dnbinom2(
          y(r),
          mu,
          variance,
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
  // SAMPLE x OTU CONDITIONAL LIKELIHOOD
  //
  // First:
  //
  //   sum_k log f(y_ijk | A = 1)
  //
  // Then marginalize A.
  // ==========================================================

  vector<Type> log_sample_count(
    n_sample_groups
  );

  log_sample_count.setZero();


  for (int r = 0;
       r < y.size();
       r++)
  {

    int s =
      row_sample_group(r);

    log_sample_count(s) +=
      log_count(r);
  }


  vector<Type> log_sample(
    n_sample_groups
  );


  for (int s = 0;
       s < n_sample_groups;
       s++)
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
    // At least one positive read
    //
    // A = 1 is required.
    //
    // L =
    // p * product_k f(y_k)
    // --------------------------------------------------------

    if (sample_positive(s) == 1)
    {

      log_sample(s) =
        log_p +
        log_sample_count(s);
    }


    // --------------------------------------------------------
    // Entire sample is zero
    //
    // A is marginalized:
    //
    // L =
    // (1-p) +
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
    }
  }


  // ==========================================================
  // SITE x OTU CONDITIONAL LIKELIHOOD GIVEN Z = 1
  // ==========================================================

  vector<Type> log_site_conditional(
    n_site_groups
  );

  log_site_conditional.setZero();


  for (int s = 0;
       s < n_sample_groups;
       s++)
  {

    int g =
      sample_site_group(s);

    log_site_conditional(g) +=
      log_sample(s);
  }


  // ==========================================================
  // OCCUPANCY MARGINALIZATION
  // ==========================================================

  vector<Type> psi(
    n_site_groups
  );

  vector<Type> log_site(
    n_site_groups
  );


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


    Type log_psi =
      log_invlogit(
        eta_occ(g)
      );


    Type log_1mpsi =
      log1m_invlogit(
        eta_occ(g)
      );


    // --------------------------------------------------------
    // Positive site history
    //
    // Z = 1 required:
    //
    // L =
    // psi * L(data | Z=1)
    // --------------------------------------------------------

    if (site_positive(g) == 1)
    {

      log_site(g) =
        log_psi +
        log_site_conditional(g);
    }


    // --------------------------------------------------------
    // Entire site history zero
    //
    // Marginalize Z:
    //
    // L =
    // (1-psi) +
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
    }


    nll -=
      log_site(g);
  }


  // ==========================================================
  // DERIVED QUANTITIES
  // ==========================================================

  vector<Type> capture_prob(
    n_sample_groups
  );


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
  }


  vector<Type> lambda(
    y.size()
  );


  for (int r = 0;
       r < y.size();
       r++)
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
    log_sample
  );

  REPORT(
    log_site
  );


  // ==========================================================
  // DELTA-METHOD DERIVED-PARAMETER SEs
  // ==========================================================

  ADREPORT(
    psi
  );

  ADREPORT(
    capture_prob
  );


  return nll;
}
