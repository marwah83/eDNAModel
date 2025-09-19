#' Prepare Long-format OTU Data from a phyloseq Object
#'
#' This function filters and reshapes a phyloseq object into a long-format data frame
#' suitable for occupancy-abundance modeling. It also prepares metadata including site,
#' treatment, sample, and replicate identifiers.
#'
#' @param physeq_obj A `phyloseq` object containing OTU and sample metadata.
#' @param min_species_sum Minimum total abundance across all samples required to retain a taxon (default: 50).
#' @param sampletype_var Column name in sample metadata indicating the sample type (default: "sampletype").
#' @param sampletype_keep The value of `sampletype_var` to retain (default: "biologicalsample").
#' @param location_var Column name in sample metadata representing site or location (default: "location").
#' @param treatment_var Column name representing the treatment group (default: "treatment").
#' @param save_path Optional path to save filtered phyloseq object (default: NULL).
#'
#' @return A list with:
#' \describe{
#'   \item{physeq_filtered}{A filtered `phyloseq` object.}
#'   \item{long_df}{A long-format data frame with columns: i (OTU), Site, Sample, Replicate, treatment, y (count).}
#' }
#'
#' @importFrom phyloseq sample_data sample_names prune_samples filter_taxa otu_table taxa_are_rows
#' @importFrom dplyr left_join select all_of mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @export
prepare_long_data <- function(physeq_obj,
                              min_species_sum = 50,
                              sampletype_var = "sampletype",
                              sampletype_keep = "biologicalsample",
                              location_var = "location",
                              treatment_var = "treatment",
                              save_path = NULL) {

  # Filter by total abundance across taxa
  physeq_filtered <- filter_phyloseq_data(
    phyloseq_obj = physeq_obj,
    min_species_sum = min_species_sum,
    save_path = save_path
  )

  # Extract sample metadata
  sample_meta <- as.data.frame(sample_data(physeq_filtered))

  if (!sampletype_var %in% colnames(sample_meta)) {
    stop(paste("Column", sampletype_var, "not found in sample_data"))
  }

  # Subset samples using metadata
  keep_samples <- rownames(sample_meta[sample_meta[[sampletype_var]] == sampletype_keep, ])
  physeq_bio <- prune_samples(keep_samples, physeq_filtered)

  # Parse Sample and Replicate
  sample_names_vec <- sample_names(physeq_bio)
  parsed <- tryCatch({
    as.data.frame(strcapture(
      pattern = "^(.*)_r([0-9]+)$",
      x = sample_names_vec,
      proto = list(Sample = character(), Replicate = integer())
    ))
  }, error = function(e) {
    message("⚠️ Pattern matching failed: ", e$message)
    data.frame(Sample = sample_names_vec, Replicate = 1)
  })

  parsed$Sample[is.na(parsed$Sample)] <- sample_names_vec[is.na(parsed$Sample)]
  parsed$Replicate[is.na(parsed$Replicate)] <- 1
  parsed$Sample <- factor(parsed$Sample)
  parsed$Replicate <- factor(parsed$Replicate)

  sample_data(physeq_bio)$Sample <- parsed$Sample
  sample_data(physeq_bio)$Replicate <- parsed$Replicate
  sample_data(physeq_bio)$SampleRep <- interaction(parsed$Sample, parsed$Replicate)

  # Clean metadata
  meta_df <- as.data.frame(sample_data(physeq_bio), stringsAsFactors = FALSE)

  if (!location_var %in% colnames(meta_df)) {
    stop(paste("Column", location_var, "not found in sample_data"))
  }

  meta_df$Site <- factor(trimws(meta_df[[location_var]]))
  meta_df[[location_var]] <- NULL

  sample_data(physeq_bio) <- phyloseq::sample_data(meta_df)

  # Filter rare taxa
  physeq_bio_filtered <- filter_taxa(physeq_bio, function(x) sum(x > 0) > 5, prune = TRUE)

  # Convert OTU table to long format
  otu_mat <- as(otu_table(physeq_bio_filtered), "matrix")
  if (taxa_are_rows(physeq_bio_filtered)) {
    otu_mat <- t(otu_mat)
  }

  otu_long <- as.data.frame(otu_mat) %>%
    tibble::rownames_to_column(var = "SampleRep") %>%
    tidyr::pivot_longer(-SampleRep, names_to = "i", values_to = "y")

  # Clean metadata for joining
  meta_df <- as.data.frame(sample_data(physeq_bio_filtered), stringsAsFactors = FALSE)
  meta_df$SampleRep <- rownames(meta_df)

  # Join OTU and metadata
  long_df <- dplyr::left_join(otu_long, meta_df, by = "SampleRep") %>%
    mutate(
      i = factor(i),
      y = as.integer(y)
    )

  # Ensure treatment variable exists
  if (!treatment_var %in% colnames(long_df)) {
    stop(paste("Column", treatment_var, "not found in metadata"))
  }

  # Select only valid output columns
  output_cols <- c("i", "Site", "Sample", "Replicate", treatment_var, "y")
  output_cols <- output_cols[output_cols %in% colnames(long_df)]
  long_df <- long_df %>% dplyr::select(all_of(output_cols))

  return(list(
    physeq_filtered = physeq_bio_filtered,
    long_df = long_df
  ))
}
