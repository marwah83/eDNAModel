#include <TMB.hpp>


// ============================================================
// NUMERICALLY STABLE LOGISTIC FUNCTIONS
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
// JOINT HIERARCHICAL eDNA MODEL
// ============================================================
//
// Hierarchy:
//
//   Site x OTU:
//       Z_im ~ Bernoulli(psi_im)
//
//   Sample x OTU:
//       A_ijm | Z_im = 1 ~ Bernoulli(p_ijm)
//
//   PCR replicate:
//       Y_ijkm | A_ijm = 1 ~ Count(lambda_ijkm)
//
// Indices:
//
//   i = site
//   j = biological sample
//   k = PCR replicate
//   m = OTU
//
// Occupancy:
//   logit(psi_im) = X_occ beta_occ + b_occ_m
//
// Capture:
//   logit(p_ijm) = X_cap beta_cap + b_cap_m
//
// Abundance:
//   log(lambda_ijkm)
//       = X_abund beta_abund
//       + offset
//       + b_abund_m
//       + b_sample_j
//       + b_sampleOTU_jm
//
// Zero inflation for ZIP/ZINB:
//   logit(pi_m) = zi_intercept + b_zi_m
//
// A and Z are analytically marginalized.
//
// Gaussian random effects can be integrated by TMB using
// the Laplace approximation when they are supplied in the
// random argument to MakeADFun().
// ============================================================


template<class Type>
Type objective_function<Type>::operator() ()
{

    // ========================================================
    // 1. DATA
    // ========================================================

    DATA_VECTOR(y);

    // Design matrices
    DATA_MATRIX(X_occ);
    DATA_MATRIX(X_cap);
    DATA_MATRIX(X_abund);

    // Abundance offset
    DATA_VECTOR(offset_abund);


    // --------------------------------------------------------
    // Hierarchical mappings
    // --------------------------------------------------------

    // site x OTU group -> OTU
    DATA_IVECTOR(occ_otu);

    // sample x OTU group -> OTU
    DATA_IVECTOR(sample_otu);

    // PCR row -> sample x OTU group
    DATA_IVECTOR(row_sample_group);

    // sample x OTU group -> site x OTU group
    DATA_IVECTOR(sample_site_group);

    // PCR row -> OTU
    DATA_IVECTOR(row_otu);

    // PCR row -> biological sample
    DATA_IVECTOR(row_sample_id);

    // PCR row -> sample x OTU RE level
    DATA_IVECTOR(row_sample_otu_re);


    // --------------------------------------------------------
    // Positive-history indicators
    // --------------------------------------------------------

    DATA_IVECTOR(sample_positive);
    DATA_IVECTOR(site_positive);


    // --------------------------------------------------------
    // Number of groups
    // --------------------------------------------------------

    DATA_INTEGER(n_site_groups);
    DATA_INTEGER(n_sample_groups);


    // --------------------------------------------------------
    // Abundance family
    // --------------------------------------------------------
    //
    // 0 = Poisson
    // 1 = NB2
    // 2 = ZIP
    // 3 = ZINB
    // --------------------------------------------------------

    DATA_INTEGER(family_code);


    // ========================================================
    // 2. RANDOM-EFFECT SWITCHES
    // ========================================================

    DATA_INTEGER(use_occ_otu);
    DATA_INTEGER(use_cap_otu);
    DATA_INTEGER(use_abund_otu);
    DATA_INTEGER(use_sample_re);
    DATA_INTEGER(use_sample_otu_re);
    DATA_INTEGER(use_zi_otu);


    // ========================================================
    // 3. VALIDATE FAMILY
    // ========================================================

    if (family_code < 0 || family_code > 3)
    {
        error("family_code must be 0, 1, 2, or 3.");
    }


    // ========================================================
    // 4. FIXED EFFECTS
    // ========================================================

    PARAMETER_VECTOR(beta_occ);
    PARAMETER_VECTOR(beta_cap);
    PARAMETER_VECTOR(beta_abund);


    // ========================================================
    // 5. RANDOM EFFECTS
    // ========================================================

    // OTU-specific occupancy RE
    PARAMETER_VECTOR(b_occ_otu);

    // OTU-specific capture RE
    PARAMETER_VECTOR(b_cap_otu);

    // OTU-specific abundance RE
    PARAMETER_VECTOR(b_abund_otu);

    // Sample-level abundance RE
    PARAMETER_VECTOR(b_sample);

    // Sample x OTU abundance RE
    PARAMETER_VECTOR(b_sample_otu);

    // OTU-specific zero-inflation RE
    PARAMETER_VECTOR(b_zi_otu);


    // ========================================================
    // 6. RANDOM-EFFECT STANDARD DEVIATIONS
    // ========================================================

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


    // ========================================================
    // 7. NB DISPERSION
    // ========================================================

    PARAMETER(log_theta);

    Type theta =
        exp(log_theta);

    // NB2:
    //
    // E(Y) = mu
    //
    // Var(Y) = mu + mu^2 / theta
    //
    // theta is currently shared among OTUs.


    // ========================================================
    // 8. ZERO-INFLATION INTERCEPT
    // ========================================================

    PARAMETER(zi_intercept);

    // ZIP/ZINB:
    //
    // logit(pi_m)
    //   =
    // zi_intercept + b_zi_otu(m)


    // ========================================================
    // 9. NEGATIVE LOG-LIKELIHOOD
    // ========================================================

    Type nll =
        Type(0);


    // ========================================================
    // 10. RANDOM-EFFECT DISTRIBUTIONS
    // ========================================================
    //
    // Random-effect blocks are currently independent.
    // ========================================================


    // --------------------------------------------------------
    // Occupancy OTU random effects
    // --------------------------------------------------------

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


    // --------------------------------------------------------
    // Capture OTU random effects
    // --------------------------------------------------------

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


    // --------------------------------------------------------
    // Abundance OTU random effects
    // --------------------------------------------------------

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


    // --------------------------------------------------------
    // Sample random effects
    // --------------------------------------------------------

    if (use_sample_re == 1)
    {
        for (int j = 0; j < b_sample.size(); j++)
        {
            nll -= dnorm(
                b_sample(j),
                Type(0),
                sd_sample,
                true
            );
        }
    }


    // --------------------------------------------------------
    // Sample x OTU random effects
    // --------------------------------------------------------

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


    // --------------------------------------------------------
    // OTU-specific zero-inflation random effects
    // --------------------------------------------------------
    //
    // Only relevant for ZIP and ZINB.
    // --------------------------------------------------------

    if (
        use_zi_otu == 1 &&
        (
            family_code == 2 ||
            family_code == 3
        )
    )
    {
        for (int m = 0; m < b_zi_otu.size(); m++)
        {
            nll -= dnorm(
                b_zi_otu(m),
                Type(0),
                sd_zi_otu,
                true
            );
        }
    }


    // ========================================================
    // 11. OCCUPANCY LINEAR PREDICTOR
    // ========================================================
    //
    // logit(psi_im)
    //   =
    // X_occ beta_occ
    // +
    // b_occ_m
    // ========================================================

    vector<Type> eta_occ =
        X_occ * beta_occ;


    if (use_occ_otu == 1)
    {
        for (int g = 0; g < n_site_groups; g++)
        {
            eta_occ(g) +=
                b_occ_otu(
                    occ_otu(g)
                );
        }
    }


    // ========================================================
    // 12. CAPTURE LINEAR PREDICTOR
    // ========================================================
    //
    // logit(p_ijm)
    //   =
    // X_cap beta_cap
    // +
    // b_cap_m
    // ========================================================

    vector<Type> eta_cap =
        X_cap * beta_cap;


    if (use_cap_otu == 1)
    {
        for (int s = 0; s < n_sample_groups; s++)
        {
            eta_cap(s) +=
                b_cap_otu(
                    sample_otu(s)
                );
        }
    }


    // ========================================================
    // 13. ABUNDANCE LINEAR PREDICTOR
    // ========================================================
    //
    // log(lambda_ijkm)
    //
    //   =
    //
    // X_abund beta_abund
    // + offset
    // + b_abund_m
    // + b_sample_j
    // + b_sampleOTU_jm
    // ========================================================

    vector<Type> eta_abund =
        X_abund * beta_abund;


    for (int r = 0; r < y.size(); r++)
    {
        // Offset
        eta_abund(r) +=
            offset_abund(r);


        // OTU abundance random effect
        if (use_abund_otu == 1)
        {
            eta_abund(r) +=
                b_abund_otu(
                    row_otu(r)
                );
        }


        // Biological-sample random effect
        if (use_sample_re == 1)
        {
            eta_abund(r) +=
                b_sample(
                    row_sample_id(r)
                );
        }


        // Sample x OTU random effect
        if (use_sample_otu_re == 1)
        {
            eta_abund(r) +=
                b_sample_otu(
                    row_sample_otu_re(r)
                );
        }
    }


    // ========================================================
    // 14. OTU-SPECIFIC ZERO-INFLATION PROBABILITIES
    // ========================================================
    //
    // logit(pi_m)
    //
    //   =
    //
    // zi_intercept
    //
    //   +
    //
    // b_zi_otu(m)
    //
    //
    // These quantities are constructed primarily for reporting.
    //
    // The likelihood itself uses log_pi and log_1m_pi directly
    // from eta_zi for numerical stability.
    // ========================================================

    vector<Type> eta_zi_otu(
        b_zi_otu.size()
    );

    vector<Type> pi_zi_otu(
        b_zi_otu.size()
    );


    for (int m = 0; m < b_zi_otu.size(); m++)
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


    // ========================================================
    // 15. PCR / READ-LEVEL COUNT LIKELIHOOD
    // ========================================================

    vector<Type> log_count(
        y.size()
    );


    for (int r = 0; r < y.size(); r++)
    {
        // ----------------------------------------------------
        // Expected abundance
        // ----------------------------------------------------

        Type mu =
            exp(
                eta_abund(r)
            );


        Type lp =
            Type(0);


        // ----------------------------------------------------
        // Zero-inflation linear predictor
        // ----------------------------------------------------

        Type eta_zi =
            zi_intercept;


        // IMPORTANT:
        //
        // The parentheses here are deliberately explicit.
        //
        // The OTU-specific ZI effect is added only when:
        //
        //   random_zi_otu is enabled
        //
        // AND
        //
        //   family is ZIP or ZINB.
        // ----------------------------------------------------

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


        // ----------------------------------------------------
        // Numerically stable ZI probabilities
        // ----------------------------------------------------
        //
        // log_pi
        //   =
        // log(pi_m)
        //
        // log_1m_pi
        //   =
        // log(1-pi_m)
        //
        // These are declared here so they remain in scope for
        // BOTH the ZIP and ZINB likelihood branches.
        // ----------------------------------------------------

        Type log_pi =
            log_invlogit(
                eta_zi
            );


        Type log_1m_pi =
            log1m_invlogit(
                eta_zi
            );


        // ====================================================
        // POISSON
        // ====================================================

        if (family_code == 0)
        {
            lp =
                dpois(
                    y(r),
                    mu,
                    true
                );
        }


        // ====================================================
        // NEGATIVE BINOMIAL 2
        // ====================================================
        //
        // E(Y) = mu
        //
        // Var(Y)
        //   =
        // mu + mu^2/theta
        // ====================================================

        else if (family_code == 1)
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


        // ====================================================
        // ZERO-INFLATED POISSON
        // ====================================================
        //
        // If y = 0:
        //
        // P(Y=0)
        //
        //   =
        //
        // pi
        //
        //   +
        //
        // (1-pi) Poisson(0 | mu)
        //
        //
        // If y > 0:
        //
        // P(Y=y)
        //
        //   =
        //
        // (1-pi) Poisson(y | mu)
        // ====================================================

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


        // ====================================================
        // ZERO-INFLATED NEGATIVE BINOMIAL 2
        // ====================================================
        //
        // If y = 0:
        //
        // P(Y=0)
        //
        //   =
        //
        // pi
        //
        //   +
        //
        // (1-pi) NB(0 | mu, theta)
        //
        //
        // If y > 0:
        //
        // P(Y=y)
        //
        //   =
        //
        // (1-pi) NB(y | mu, theta)
        // ====================================================

        else if (family_code == 3)
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


            if (y(r) == Type(0))
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


        // Store conditional PCR/read-level log likelihood.

        log_count(r) =
            lp;
    }


    // ========================================================
    // 16. COMBINE PCR REPLICATES WITHIN SAMPLE x OTU
    // ========================================================
    //
    // Conditional on A_ijm = 1:
    //
    // L_count_ijm
    //
    //   =
    //
    // product_k f(y_ijkm)
    //
    //
    // Therefore:
    //
    // log L_count_ijm
    //
    //   =
    //
    // sum_k log f(y_ijkm)
    // ========================================================

    vector<Type> log_sample_count(
        n_sample_groups
    );


    log_sample_count.setZero();


    for (int r = 0; r < y.size(); r++)
    {
        int s =
            row_sample_group(r);


        log_sample_count(s) +=
            log_count(r);
    }


    // ========================================================
    // 17. MARGINALIZE SAMPLE-LEVEL CAPTURE STATE A
    // ========================================================
    //
    // Capture lives at SAMPLE x OTU.
    //
    //
    // POSITIVE SAMPLE:
    //
    // If at least one PCR count is > 0:
    //
    //   A_ijm = 1
    //
    // and:
    //
    //   L_ijm
    //
    //     =
    //
    //   p_ijm *
    //   product_k f(y_ijkm)
    //
    //
    // ALL-ZERO SAMPLE:
    //
    // Two possibilities:
    //
    //   A_ijm = 0
    //
    // OR
    //
    //   A_ijm = 1 and all PCR counts are zero.
    //
    //
    // Therefore:
    //
    // L_ijm
    //
    //   =
    //
    // (1-p_ijm)
    //
    //   +
    //
    // p_ijm *
    // product_k f(0)
    // ========================================================

    vector<Type> log_sample(
        n_sample_groups
    );


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


        // ----------------------------------------------------
        // At least one positive PCR replicate
        // ----------------------------------------------------

        if (sample_positive(s) == 1)
        {
            log_sample(s) =
                log_p +
                log_sample_count(s);
        }


        // ----------------------------------------------------
        // All PCR replicates are zero
        // ----------------------------------------------------

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


    // ========================================================
    // 18. COMBINE SAMPLES WITHIN SITE x OTU
    // ========================================================
    //
    // Conditional on:
    //
    //   Z_im = 1
    //
    // we have:
    //
    // L(data_im | Z_im=1)
    //
    //   =
    //
    // product_j L_ijm
    //
    //
    // Hence:
    //
    // log L(data_im | Z_im=1)
    //
    //   =
    //
    // sum_j log L_ijm
    // ========================================================

    vector<Type> log_site_conditional(
        n_site_groups
    );


    log_site_conditional.setZero();


    for (int s = 0; s < n_sample_groups; s++)
    {
        int g =
            sample_site_group(s);


        log_site_conditional(g) +=
            log_sample(s);
    }


    // ========================================================
    // 19. MARGINALIZE SITE-LEVEL OCCUPANCY STATE Z
    // ========================================================
    //
    // IMPORTANT:
    //
    // There is ONE occupancy state:
    //
    //   Z_im
    //
    // per SITE x OTU.
    //
    // Occupancy therefore does NOT live at the PCR/read level.
    //
    //
    // POSITIVE SITE x OTU:
    //
    // If anything is positive:
    //
    //   Z_im = 1
    //
    // so:
    //
    // L_im
    //
    //   =
    //
    // psi_im *
    // L(data_im | Z_im=1)
    //
    //
    // ALL-ZERO SITE x OTU:
    //
    // Either:
    //
    //   Z_im = 0
    //
    // OR
    //
    //   Z_im = 1 but all lower-level observations are zero.
    //
    //
    // Therefore:
    //
    // L_im
    //
    //   =
    //
    // (1-psi_im)
    //
    //   +
    //
    // psi_im *
    // L(all-zero data | Z_im=1)
    // ========================================================

    vector<Type> psi(
        n_site_groups
    );


    vector<Type> log_site(
        n_site_groups
    );


    for (int g = 0; g < n_site_groups; g++)
    {
        // Occupancy probability on natural scale.

        psi(g) =
            exp(
                log_invlogit(
                    eta_occ(g)
                )
            );


        // Stable log(psi)

        Type log_psi =
            log_invlogit(
                eta_occ(g)
            );


        // Stable log(1-psi)

        Type log_1mpsi =
            log1m_invlogit(
                eta_occ(g)
            );


        // ----------------------------------------------------
        // Positive site x OTU
        // ----------------------------------------------------

        if (site_positive(g) == 1)
        {
            log_site(g) =
                log_psi +
                log_site_conditional(g);
        }


        // ----------------------------------------------------
        // All-zero site x OTU
        // ----------------------------------------------------

        else
        {
            log_site(g) =
                logspace_add(
                    log_1mpsi,
                    log_psi +
                        log_site_conditional(g)
                );
        }


        // Add contribution to negative log likelihood.

        nll -=
            log_site(g);
    }


    // ========================================================
    // 20. DERIVED CAPTURE PROBABILITIES
    // ========================================================

    vector<Type> capture_prob(
        n_sample_groups
    );


    for (int s = 0; s < n_sample_groups; s++)
    {
        capture_prob(s) =
            exp(
                log_invlogit(
                    eta_cap(s)
                )
            );
    }


    // ========================================================
    // 21. DERIVED ABUNDANCE
    // ========================================================

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


    // ========================================================
    // 22. REPORT
    // ========================================================

    REPORT(eta_occ);

    REPORT(eta_cap);

    REPORT(eta_abund);

    REPORT(eta_zi_otu);

    REPORT(psi);

    REPORT(capture_prob);

    REPORT(lambda);

    REPORT(pi_zi_otu);

    REPORT(log_count);

    REPORT(log_sample_count);

    REPORT(log_sample);

    REPORT(log_site_conditional);

    REPORT(log_site);


    // ========================================================
    // 23. ADREPORT
    // ========================================================
    //
    // Allows sdreport() to obtain approximate uncertainty for
    // derived quantities through automatic differentiation.
    // ========================================================

    ADREPORT(psi);

    ADREPORT(capture_prob);


    if (
        family_code == 2 ||
        family_code == 3
    )
    {
        ADREPORT(pi_zi_otu);
    }


    // ========================================================
    // 24. RETURN NEGATIVE LOG-LIKELIHOOD
    // ========================================================

    return nll;
}
