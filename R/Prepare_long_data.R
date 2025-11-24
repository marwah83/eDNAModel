#' Prepare long-format OTU data with inferred sample and replicate columns
#'
#' Converts a phyloseq object into a long-format data frame for modeling.
#'
#' @param physeq_obj A phyloseq object.
#' @param min_species_sum Minimum total abundance required to keep a species.
#' @param site_col Name of the column indicating the site or location.
#' @param nested_cols (Optional) Character vector of column names that are nested with OTU (e.g., sample, replicate).
#'
#' @return A list with:
#' \describe{
#'   \item{physeq_filtered}{Filtered phyloseq object.}
#'   \item{long_df}{Long-format data frame.}
#'   \item{sample_col}{Auto-detected sample column name.}
#'   \item{replicate_col}{Auto-detected replicate column name.}
#' }
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
