#' Fit Multi-Species Occupancy-Abundance Model from a phyloseq Object
#'
#' Prepares long-format data from a `phyloseq` object, builds species-interacted
#' GLMM formulas, and fits Poisson (abundance) and Binomial (occupancy) models
#' using an iterative burn-in approach.
#'
#' @param phyloseq A `phyloseq` object containing OTU, taxonomy, and sample metadata.
#' @param poisson_rhs A `formula` object defining the right-hand side of the Poisson model.
#' @param binomial_rhs A `formula` object defining the right-hand side of the Binomial model.
#' @param min_species_sum Minimum total count required to retain a taxon (default: 50).
#' @param sampletype_keep Optional character value. If specified, filters the data to keep only samples with this `sampletype`.
#' @param num_iterations Total number of iterations to run the model fitting (default: 100).
#' @param burn_in Number of iterations to discard as burn-in (default: 50).
#'
#' @return A list with fitted Poisson and Binomial model objects for each post-burn-in iteration.
#'
#' @examples
#' \dontrun{
#' data(physeq_example)
#' result <- FitModel(
#'   phyloseq = physeq_example,
#'   poisson_rhs = ~ (1 | Site) + (1 | Sample) + (1 | Replicate) + treatment,
#'   binomial_rhs = ~ (1 | Site) + treatment,
#'   sampletype_keep = "biologicalsample"
#' )
#' }
#'
#' @export
FitModel <- function(phyloseq,
                     poisson_rhs,
                     binomial_rhs,
                     min_species_sum = 50,
                     sampletype_keep = NULL,
                     num_iterations = 100,
                     burn_in = 50) {

  # Step 1: Prepare long-format data
  prep <- prepare_long_data(
    physeq_obj = phyloseq,
    min_species_sum = min_species_sum,
    sampletype_keep = sampletype_keep
  )
  long_df <- prep$long_df
  replicate_var <- "Replicate"

  # Step 2: Automatically detect most abundant OTU (based on the filtered object)
  message("🔍 Calculating most abundant OTU from filtered phyloseq object...")
  otu_abundances <- taxa_sums(prep$physeq_filtered)
  top_otu <- names(sort(otu_abundances, decreasing = TRUE))[1]
  message("📌 Using most abundant OTU as reference: ", top_otu)

  long_df$OTU <- relevel(factor(long_df$OTU), ref = top_otu)

  # Step 3: Optional - clean treatment variable
  if ("treatment" %in% names(long_df)) {
    long_df <- long_df %>%
      dplyr::filter(treatment != "0")
    long_df$treatment <- droplevels(factor(long_df$treatment))
  }

  # Step 4: Build model formulas with species interaction (except replicate)
  poisson_formula <- build_formula_with_species_interaction(
    rhs = poisson_rhs,
    response = "y",
    species_var = "OTU",
    lower_level = replicate_var
  )

  binomial_formula <- build_formula_with_species_interaction(
    rhs = binomial_rhs,
    response = "z_sim",
    species_var = "OTU",
    lower_level = replicate_var
  )

  message("📌 Poisson model formula: ", deparse(poisson_formula))
  message("📌 Binomial model formula: ", deparse(binomial_formula))

  # Step 5: Fit the model using iterative GLMMs
  result <- simulate_glm_burnin_iterations(
    data_glm = long_df,
    poisson_formula = poisson_formula,
    binomial_formula = binomial_formula,
    num_iterations = num_iterations,
    burn_in = burn_in
  )

  # Step 6: Return model output
  return(list(
    poisson_models = result$poisson_models,
    binomial_models = result$binomial_models,
    data_glm = long_df  # for inspection or reuse
  ))
}
