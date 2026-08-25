FitModel_joint <- function(
  phyloseq,
  site_col,
  sample_col = "Name",
  replicate_col = NULL,
  otu_col = "OTU",
  count_col = "y",

  occupancy_formula = ~ 1,
  capture_formula = ~ 1,
  abundance_formula = ~ 1,

  abundance_offset = NULL,

  abundance_family = c(
    "poisson",
    "nbinom",
    "zip",
    "zinb"
  ),

  min_species_sum = 10,
  min_detection_replicates = 1,

  random_occ_otu = TRUE,
  random_capture_otu = TRUE,
  random_abund_otu = TRUE,
  random_sample = TRUE,
  random_sample_otu = TRUE,

  verbose = TRUE
) {

  abundance_family <-
    match.arg(
      abundance_family
    )

  # ============================================================
  # 1. Prepare long data
  # ============================================================

  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    site_col = site_col,
    nested_cols = unique(
      stats::na.omit(
        c(
          sample_col,
          replicate_col
        )
      )
    )
  )

  dat <- prep$long_df


  # ============================================================
  # 2. Required columns
  # ============================================================

  required <- unique(
    c(
      site_col,
      sample_col,
      otu_col,
      count_col
    )
  )

  missing_cols <- setdiff(
    required,
    names(dat)
  )

  if (length(missing_cols) > 0) {

    stop(
      "Missing required columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }


  # ============================================================
  # 3. Basic types
  # ============================================================

  dat[[site_col]] <-
    as.character(
      dat[[site_col]]
    )

  dat[[sample_col]] <-
    as.character(
      dat[[sample_col]]
    )

  dat[[otu_col]] <-
    as.character(
      dat[[otu_col]]
    )

  dat[[count_col]] <-
    as.numeric(
      dat[[count_col]]
    )


  dat <- dat[
    is.finite(
      dat[[count_col]]
    ),
    ,
    drop = FALSE
  ]


  # ============================================================
  # IMPORTANT
  #
  # The joint likelihood uses the actual zero/non-zero count
  # process. Therefore zero is the detection boundary.
  # ============================================================

  if (
    any(
      dat[[count_col]] < 0,
      na.rm = TRUE
    )
  ) {

    stop(
      "Counts must be non-negative."
    )
  }


  # ============================================================
  # 4. OTU filtering
  # ============================================================

  otu_stats <- dat |>
    dplyr::group_by(
      .data[[otu_col]]
    ) |>
    dplyr::summarise(

      total_count =
        sum(
          .data[[count_col]],
          na.rm = TRUE
        ),

      detections =
        sum(
          .data[[count_col]] > 0,
          na.rm = TRUE
        ),

      .groups = "drop"
    )


  retained_otus <- otu_stats |>
    dplyr::filter(

      .data$total_count >=
        min_species_sum,

      .data$detections >=
        min_detection_replicates

    ) |>
    dplyr::pull(
      .data[[otu_col]]
    )


  dat <- dat |>
    dplyr::filter(
      .data[[otu_col]] %in%
        retained_otus
    )


  if (nrow(dat) == 0) {

    stop(
      "No taxa remain after filtering."
    )
  }


  # ============================================================
  # 5. Build index factors
  # ============================================================

  dat$.site <-
    factor(
      dat[[site_col]]
    )

  dat$.otu <-
    factor(
      dat[[otu_col]]
    )

  dat$.sample <-
    factor(
      dat[[sample_col]]
    )


  dat$.site_otu <-
    interaction(
      dat$.site,
      dat$.otu,
      drop = TRUE,
      lex.order = TRUE
    )


  dat$.sample_otu <-
    interaction(
      dat$.site,
      dat$.sample,
      dat$.otu,
      drop = TRUE,
      lex.order = TRUE
    )


  # ============================================================
  # 6. Sample x OTU data
  # ============================================================

  sample_df <- dat |>
    dplyr::group_by(
      .data$.sample_otu
    ) |>
    dplyr::summarise(

      .site_otu =
        dplyr::first(
          .data$.site_otu
        ),

      .otu =
        dplyr::first(
          .data$.otu
        ),

      .sample =
        dplyr::first(
          .data$.sample
        ),

      sample_positive =
        as.integer(
          any(
            .data[[count_col]] > 0
          )
        ),

      dplyr::across(
        dplyr::everything(),
        ~ dplyr::first(.x)
      ),

      .groups = "drop"
    )


  # ============================================================
  # 7. Site x OTU data
  # ============================================================

  site_df <- dat |>
    dplyr::group_by(
      .data$.site_otu
    ) |>
    dplyr::summarise(

      .otu =
        dplyr::first(
          .data$.otu
        ),

      site_positive =
        as.integer(
          any(
            .data[[count_col]] > 0
          )
        ),

      dplyr::across(
        dplyr::everything(),
        ~ dplyr::first(.x)
      ),

      .groups = "drop"
    )


  # ============================================================
  # 8. Design matrices
  # ============================================================

  X_occ <- stats::model.matrix(
    occupancy_formula,
    data = site_df
  )


  X_cap <- stats::model.matrix(
    capture_formula,
    data = sample_df
  )


  X_abund <- stats::model.matrix(
    abundance_formula,
    data = dat
  )


  # ============================================================
  # 9. Offset
  # ============================================================

  if (is.null(abundance_offset)) {

    offset_abund <-
      rep(
        0,
        nrow(dat)
      )

  } else {

    if (
      !(abundance_offset %in%
        names(dat))
    ) {

      stop(
        "abundance_offset '",
        abundance_offset,
        "' not found."
      )
    }


    offset_value <-
      as.numeric(
        dat[[abundance_offset]]
      )


    if (
      any(
        !is.finite(offset_value) |
          offset_value <= 0
      )
    ) {

      stop(
        "abundance_offset must be positive."
      )
    }


    offset_abund <-
      log(
        offset_value
      )
  }


  # ============================================================
  # 10. Integer mappings
  # ============================================================

  site_levels <-
    levels(
      dat$.site_otu
    )

  sample_levels <-
    levels(
      dat$.sample_otu
    )


  row_sample_group <-
    match(
      dat$.sample_otu,
      sample_levels
    ) - 1L


  sample_site_group <-
    match(
      sample_df$.site_otu,
      site_levels
    ) - 1L


  occ_otu <-
    as.integer(
      site_df$.otu
    ) - 1L


  sample_otu <-
    as.integer(
      sample_df$.otu
    ) - 1L


  row_otu <-
    as.integer(
      dat$.otu
    ) - 1L


  sample_id <-
    as.integer(
      sample_df$.sample
    ) - 1L


  row_sample_id <-
    as.integer(
      dat$.sample
    ) - 1L


  sample_otu_re <-
    seq_len(
      nrow(sample_df)
    ) - 1L


  row_sample_otu_re <-
    row_sample_group


  # ============================================================
  # 11. Family code
  # ============================================================

  family_code <- switch(

    abundance_family,

    poisson = 0L,

    nbinom = 1L,

    zip = 2L,

    zinb = 3L
  )


  # ============================================================
  # 12. TMB data
  # ============================================================

  data_tmb <- list(

    y =
      dat[[count_col]],

    X_occ =
      X_occ,

    X_cap =
      X_cap,

    X_abund =
      X_abund,

    offset_abund =
      offset_abund,

    occ_otu =
      occ_otu,

    sample_otu =
      sample_otu,

    sample_id =
      sample_id,

    sample_otu_re =
      sample_otu_re,

    row_sample_group =
      row_sample_group,

    sample_site_group =
      sample_site_group,

    row_otu =
      row_otu,

    row_sample_id =
      row_sample_id,

    row_sample_otu_re =
      row_sample_otu_re,

    sample_positive =
      as.integer(
        sample_df$sample_positive
      ),

    site_positive =
      as.integer(
        site_df$site_positive
      ),

    n_site_groups =
      nrow(site_df),

    n_sample_groups =
      nrow(sample_df),

    family_code =
      family_code,

    use_occ_otu =
      as.integer(
        random_occ_otu
      ),

    use_cap_otu =
      as.integer(
        random_capture_otu
      ),

    use_abund_otu =
      as.integer(
        random_abund_otu
      ),

    use_sample_re =
      as.integer(
        random_sample
      ),

    use_sample_otu_re =
      as.integer(
        random_sample_otu
      )
  )


  # ============================================================
  # 13. Starting parameters
  # ============================================================

  n_otu <-
    nlevels(
      dat$.otu
    )

  n_sample <-
    nlevels(
      dat$.sample
    )

  n_sample_otu <-
    nrow(
      sample_df
    )


  parameters <- list(

    beta_occ =
      rep(
        0,
        ncol(X_occ)
      ),

    beta_cap =
      rep(
        0,
        ncol(X_cap)
      ),

    beta_abund =
      rep(
        0,
        ncol(X_abund)
      ),

    b_occ_otu =
      rep(
        0,
        n_otu
      ),

    b_cap_otu =
      rep(
        0,
        n_otu
      ),

    b_abund_otu =
      rep(
        0,
        n_otu
      ),

    b_sample =
      rep(
        0,
        n_sample
      ),

    b_sample_otu =
      rep(
        0,
        n_sample_otu
      ),

    log_sd_occ_otu =
      log(0.5),

    log_sd_cap_otu =
      log(0.5),

    log_sd_abund_otu =
      log(0.5),

    log_sd_sample =
      log(0.5),

    log_sd_sample_otu =
      log(0.5),

    log_theta =
      log(10),

    zi_intercept =
      -3
  )


  # ============================================================
  # 14. Random-effect list
  # ============================================================

  random_effects <-
    character(0)


  if (random_occ_otu) {

    random_effects <-
      c(
        random_effects,
        "b_occ_otu"
      )
  }


  if (random_capture_otu) {

    random_effects <-
      c(
        random_effects,
        "b_cap_otu"
      )
  }


  if (random_abund_otu) {

    random_effects <-
      c(
        random_effects,
        "b_abund_otu"
      )
  }


  if (random_sample) {

    random_effects <-
      c(
        random_effects,
        "b_sample"
      )
  }


  if (random_sample_otu) {

    random_effects <-
      c(
        random_effects,
        "b_sample_otu"
      )
  }


  # ============================================================
  # 15. Build joint TMB model
  # ============================================================

  obj <- TMB::MakeADFun(

    data =
      data_tmb,

    parameters =
      parameters,

    random =
      random_effects,

    DLL =
      "edna_joint",

    silent =
      !verbose
  )


  # ============================================================
  # 16. Optimize
  # ============================================================

  opt <- stats::nlminb(

    start =
      obj$par,

    objective =
      obj$fn,

    gradient =
      obj$gr,

    control = list(

      iter.max =
        2000,

      eval.max =
        3000,

      rel.tol =
        1e-10
    )
  )


  # ============================================================
  # 17. Joint covariance / SE
  # ============================================================

  sdr <- TMB::sdreport(
    obj,
    par.fixed =
      opt$par,

    getJointPrecision =
      TRUE
  )


  # ============================================================
  # 18. Reports
  # ============================================================

  report <-
    obj$report(
      opt$par
    )


  # ============================================================
  # 19. Parameter table
  # ============================================================

  fixed_summary <-
    as.data.frame(
      summary(
        sdr,
        "fixed"
      )
    )


  fixed_summary$parameter <-
    rownames(
      fixed_summary
    )


  rownames(
    fixed_summary
  ) <- NULL


  # ============================================================
  # 20. Derived parameter table
  # ============================================================

  ad_summary <-
    tryCatch(

      as.data.frame(
        summary(
          sdr,
          "report"
        )
      ),

      error = function(e)
        data.frame()
    )


  # ============================================================
  # 21. Return
  # ============================================================

  list(

    fit =
      opt,

    tmb_object =
      obj,

    sdreport =
      sdr,

    report =
      report,

    fixed_effects =
      fixed_summary,

    derived =
      ad_summary,

    site_data =
      site_df,

    sample_data =
      sample_df,

    long_df =
      dat,

    otu_stats =
      otu_stats,

    retained_otus =
      retained_otus,

    formulas = list(

      occupancy =
        occupancy_formula,

      capture =
        capture_formula,

      abundance =
        abundance_formula
    ),

    abundance_family =
      abundance_family,

    convergence = list(

      code =
        opt$convergence,

      message =
        opt$message,

      objective =
        opt$objective
    ),

    note = paste(
      "Joint observed-data likelihood.",
      "Latent occupancy Z and biological-sample capture A",
      "are analytically marginalized rather than imputed."
    )
  )
}
