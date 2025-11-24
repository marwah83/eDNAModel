#' Prepare Long-format OTU Data from Phyloseq Object
#'
#' Filters rare taxa, auto-detects sample and replicate columns, and returns long-format data
#' with OTU abundances and metadata merged. Assumes presence of consistent metadata columns.
#'
#' @param physeq_obj A \code{phyloseq} object containing OTU table and sample metadata.
#' @param min_species_sum Integer. Minimum total abundance for a species (OTU) to be retained.
#' @param site_col Character string. Name of the column in metadata corresponding to site identity.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{physeq_filtered}}{Filtered phyloseq object with rare taxa removed}
#'   \item{\code{long_df}}{Long-format data frame with OTU abundances and metadata}
#'   \item{\code{sample_col}}{Auto-detected column used as sample ID}
#'   \item{\code{replicate_col}}{Auto-detected column used as replicate ID}
#' }
#' 
#'
#' @importFrom phyloseq sample_data otu_table taxa_are_rows filter_taxa
#' @importFrom dplyr mutate left_join across
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
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
