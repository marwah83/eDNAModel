#' Prepare Long Format Data from a Phyloseq Object
#'
#' Converts a `phyloseq` object into a tidy long-format data frame for modeling, with options for filtering taxa and creating nested replicate groupings.
#'
#' @param physeq_obj A `phyloseq` object containing OTU (or ASV) abundance table and sample metadata.
#' @param min_species_sum Integer. Minimum total count threshold across all samples for a taxon to be retained. Default is 50.
#' @param site_col Character. Column name in `sample_data` representing the sampling site variable. This will be used to group data for site-level modeling.
#' @param nested_cols Optional character vector. Column names in `sample_data` to be combined into a nested grouping factor (`SampleRep`). Used for hierarchical or repeated measures design. Default is `NULL`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{physeq_filtered}{A filtered `phyloseq` object with low-abundance and low-prevalence taxa removed.}
#'   \item{long_df}{A long-format `data.frame` with OTU counts and merged metadata. Contains one row per `SampleRep × OTU`. Columns include `SampleRep`, `OTU`, `y` (abundance), and all metadata variables.}
#' }
#'
#' @details
#' This function is designed to transform a `phyloseq` object into a format suitable for downstream generalized linear mixed modeling (GLMM). Key steps include:
#' \itemize{
#'   \item Taxa are filtered by total count (`min_species_sum`) and then by prevalence (present in more than 5 samples).
#'   \item A `SampleRep` variable is constructed from `nested_cols` if provided, otherwise defaults to sample names.
#'   \item The OTU table is melted into long format and joined with sample metadata using `SampleRep`.
#'   \item A new column `Site` is created using the `site_col`, for grouping and site-level modeling.
#' }
#'
#' The resulting data frame (`long_df`) can be used for modeling occupancy and abundance in ecological or microbiome analyses, including two-part models like those implemented in `FitModel()`.
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#' long_data <- prepare_long_data(GlobalPatterns, site_col = "SampleID", min_species_sum = 50)
#' }
#'
#' @seealso \code{\link{FitModel}} for modeling based on this long-format data.
#' @importFrom phyloseq filter_taxa sample_data taxa_are_rows otu_table sample_names
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @export
prepare_long_data <- function(physeq_obj,
                              min_species_sum = 50,
                              site_col,
                              nested_cols = NULL) {
  
  # Filter taxa by total abundance
  physeq_filtered <- filter_phyloseq_data(physeq_obj, min_species_sum = min_species_sum)
  
  # Extract sample metadata
  sample_meta <- as.data.frame(sample_data(physeq_filtered), stringsAsFactors = FALSE)

  # Create SampleRep based on nested_cols (if any), else use sample_names
  if (!is.null(nested_cols)) {
    if (!all(nested_cols %in% names(sample_meta))) {
      stop("Some nested_cols not found in sample_data.")
    }
    sample_data(physeq_filtered)$SampleRep <- do.call(interaction, sample_meta[, nested_cols, drop = FALSE])
  } else {
    sample_data(physeq_filtered)$SampleRep <- sample_names(physeq_filtered)
  }

  # Update metadata with SampleRep
  meta_df <- as.data.frame(sample_data(physeq_filtered), stringsAsFactors = FALSE)
  meta_df$SampleRep <- rownames(meta_df)

  # Filter low-prevalence taxa
  physeq_filtered <- filter_taxa(physeq_filtered, function(x) sum(x > 0) > 5, prune = TRUE)

  # OTU matrix to long format
  otu_mat <- as(otu_table(physeq_filtered), "matrix")
  if (taxa_are_rows(physeq_filtered)) {
    otu_mat <- t(otu_mat)
  }

  otu_long <- as.data.frame(otu_mat) %>%
    tibble::rownames_to_column("SampleRep") %>%
    tidyr::pivot_longer(-SampleRep, names_to = "OTU", values_to = "y")

  # Merge OTU and metadata
  long_df <- left_join(otu_long, meta_df, by = "SampleRep") %>%
    mutate(
      OTU = factor(OTU),
      y = as.integer(y),
      Site = .data[[site_col]]
    )

  return(list(
    physeq_filtered = physeq_filtered,
    long_df = long_df
  ))
}
