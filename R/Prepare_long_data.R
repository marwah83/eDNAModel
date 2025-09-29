#' Prepare Long-format OTU Data from a phyloseq Object
#'
#' Filters and reshapes a `phyloseq` object into a long-format data frame for
#' occupancy-abundance modeling. The function optionally filters on a sample type
#' (e.g., `"biologicalsample"`), parses replicate and sample identifiers, renames
#' location to `Site`, and transforms the OTU matrix to long format.
#'
#' @param physeq_obj A [`phyloseq::phyloseq`] object containing OTU table and sample metadata.
#' @param min_species_sum Minimum total abundance required to retain a taxon (default: 50).
#' @param sampletype_keep Optional value of `sampletype` column to filter by (e.g., `"biologicalsample"`).
#'
#' @return A list with:
#' \describe{
#'   \item{physeq_filtered}{Filtered `phyloseq` object after abundance and sample type filtering.}
#'   \item{long_df}{Long-format `data.frame` with OTU counts, metadata, and parsed identifiers.}
#' }
#' 
#' @importFrom phyloseq sample_data otu_table taxa_are_rows filter_taxa prune_samples sample_names
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr pivot_longer
#' 
#' @export
prepare_long_data <- function(physeq_obj,
                              min_species_sum = 50,
                              sampletype_keep = NULL) {
  
  # Filter by total abundance
  physeq_filtered <- filter_phyloseq_data(physeq_obj, min_species_sum = min_species_sum)
  
  # Extract sample metadata
  sample_meta <- as.data.frame(sample_data(physeq_filtered), stringsAsFactors = FALSE)
  
  # Optional filtering by sampletype
  if (!is.null(sampletype_keep)) {
    if (!"sampletype" %in% names(sample_meta)) {
      stop("sampletype column missing in sample_data")
    }
    keep_samples <- rownames(sample_meta[sample_meta$sampletype == sampletype_keep, ])
    physeq_filtered <- prune_samples(keep_samples, physeq_filtered)
  }
  
  # Parse replicate and sample
  sample_names_vec <- sample_names(physeq_filtered)
  parsed <- strcapture("^(.*)_r([0-9]+)$", sample_names_vec,
                       proto = list(Sample = character(), Replicate = integer()))
  parsed$Sample[is.na(parsed$Sample)] <- sample_names_vec[is.na(parsed$Sample)]
  parsed$Replicate[is.na(parsed$Replicate)] <- 1
  parsed$Sample <- factor(parsed$Sample)
  parsed$Replicate <- factor(parsed$Replicate)
  
  sample_data(physeq_filtered)$Sample <- parsed$Sample
  sample_data(physeq_filtered)$Replicate <- parsed$Replicate
  sample_data(physeq_filtered)$SampleRep <- interaction(parsed$Sample, parsed$Replicate)
  
  # Pull out metadata again, and rename "location" to "Site"
  meta_df <- as.data.frame(sample_data(physeq_filtered), stringsAsFactors = FALSE)
  if (!"location" %in% names(meta_df)) {
    stop("location column missing in sample_data")
  }
  meta_df$Site <- factor(trimws(meta_df$location))
  meta_df$SampleRep <- rownames(meta_df)
  
  # Filter low-prevalence taxa
  physeq_filtered <- filter_taxa(physeq_filtered, function(x) sum(x > 0) > 5, prune = TRUE)
  
  # Prepare OTU long format
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
      y = as.integer(y)
    )
  
  return(list(
    physeq_filtered = physeq_filtered,
    long_df = long_df
  ))
}
