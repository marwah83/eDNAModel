#' Prepare Long-format OTU Data from a phyloseq Object
#'
#' Filters and reshapes a phyloseq object into long format for occupancy-abundance modeling.
#'
#' @param physeq_obj A `phyloseq` object.
#' @param min_species_sum Minimum abundance to retain a taxon (default: 50).
#' @param sampletype_keep The value of `sampletype` to retain (default: "biologicalsample").
#' @param col_vars A named list mapping the following required fields:
#'   - "sampletype"
#'   - "location"
#'   - "treatment"
#'   - "sample"
#'   - "samplerep"
#'   Any other variables in the list will also be carried through to output.
#' @param save_path Optional path to save the filtered object.
#'
#' @return A list with:
#' \describe{
#'   \item{physeq_filtered}{Filtered `phyloseq` object}
#'   \item{long_df}{Long-format data frame with counts and metadata}
#' }
#' @export
prepare_long_data <- function(physeq_obj,
                              min_species_sum = 50,
                              sampletype_keep = "biologicalsample",
                              col_vars = list(
                                sampletype = "sampletype",
                                location = "location",
                                treatment = "treatment",
                                sample = "sample",
                                samplerep = "samplerep"
                              ),
                              save_path = NULL) {

  # ==== Validate col_vars ====
  required_keys <- c("sampletype", "location", "treatment", "sample", "samplerep")
  missing_keys <- setdiff(required_keys, names(col_vars))
  if (length(missing_keys) > 0) {
    stop("Missing required `col_vars`: ", paste(missing_keys, collapse = ", "))
  }

  # ==== Filter by abundance ====
  physeq_filtered <- filter_phyloseq_data(
    phyloseq_obj = physeq_obj,
    min_species_sum = min_species_sum,
    save_path = save_path
  )

  sample_meta <- as.data.frame(sample_data(physeq_filtered), stringsAsFactors = FALSE)

  # ==== Filter by sampletype ====
  sampletype_col <- col_vars$sampletype
  if (!sampletype_col %in% colnames(sample_meta)) {
    stop("Column '", sampletype_col, "' not found in sample_data.")
  }

  keep_samples <- rownames(sample_meta[sample_meta[[sampletype_col]] == sampletype_keep, ])
  physeq_bio <- prune_samples(keep_samples, physeq_filtered)

  # ==== Extract & assign sample + replicate ====
  sample_names_vec <- sample_names(physeq_bio)
  parsed <- tryCatch({
    strcapture(
      pattern = "^(.*)_r([0-9]+)$",
      x = sample_names_vec,
      proto = list(Sample = character(), Replicate = integer())
    )
  }, error = function(e) {
    warning("⚠️ Failed to parse replicates: ", e$message)
    data.frame(Sample = sample_names_vec, Replicate = 1)
  })

  parsed$Sample[is.na(parsed$Sample)] <- sample_names_vec[is.na(parsed$Sample)]
  parsed$Replicate[is.na(parsed$Replicate)] <- 1

  parsed$Sample <- factor(parsed$Sample)
  parsed$Replicate <- factor(parsed$Replicate)

  sample_data(physeq_bio)$Sample <- parsed$Sample
  sample_data(physeq_bio)$Replicate <- parsed$Replicate
  sample_data(physeq_bio)$SampleRep <- interaction(parsed$Sample, parsed$Replicate)

  # ==== Prepare metadata ====
  meta_df <- as.data.frame(sample_data(physeq_bio), stringsAsFactors = FALSE)
  meta_df$Site <- factor(trimws(meta_df[[col_vars$location]]))
  meta_df[[col_vars$location]] <- NULL
  sample_data(physeq_bio) <- phyloseq::sample_data(meta_df)

  # ==== Filter rare taxa ====
  physeq_bio_filtered <- filter_taxa(physeq_bio, function(x) sum(x > 0) > 5, prune = TRUE)

  # ==== Convert to long format ====
  otu_mat <- as(otu_table(physeq_bio_filtered), "matrix")
  if (taxa_are_rows(physeq_bio_filtered)) otu_mat <- t(otu_mat)

  otu_long <- as.data.frame(otu_mat) %>%
    tibble::rownames_to_column(var = "SampleRep") %>%
    tidyr::pivot_longer(-SampleRep, names_to = "i", values_to = "y")

  # ==== Clean metadata & join ====
  meta_df <- as.data.frame(sample_data(physeq_bio_filtered), stringsAsFactors = FALSE)
  meta_df$SampleRep <- rownames(meta_df)

  long_df <- dplyr::left_join(otu_long, meta_df, by = "SampleRep") %>%
    mutate(
      i = factor(i),
      y = as.integer(y)
    )

  # ==== Build column list ====
  output_cols <- c(
    "i", "Site", "Sample", "Replicate", col_vars$treatment,
    setdiff(names(col_vars), required_keys)  # include any extra columns
  )
  output_cols <- unique(output_cols[output_cols %in% colnames(long_df)])

  long_df <- long_df %>% dplyr::select(all_of(output_cols), y)

  # ==== Final Check ====
  if (nrow(long_df) == 0) {
    warning("No data remaining after filtering. Check sampletype_keep, min_species_sum, or col_vars.")
  }

  return(list(
    physeq_filtered = physeq_bio_filtered,
    long_df = long_df
  ))
}
