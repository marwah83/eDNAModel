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
                              site_col) {

  # --- Step 1: Filter low-abundance taxa ---
  physeq_filtered <- filter_phyloseq_data(physeq_obj, min_species_sum = min_species_sum)

  # --- Step 2: Extract metadata ---
  sample_meta <- as.data.frame(sample_data(physeq_filtered), stringsAsFactors = FALSE)

  # --- Step 3: Auto-detect replicate and sample columns ---
  col_names <- colnames(sample_meta)

  replicate_col_candidates <- col_names[grepl("rep", col_names, ignore.case = TRUE)]
  sample_col_candidates    <- col_names[grepl("name|id|sample", col_names, ignore.case = TRUE)]

  # Ensure sample != replicate
  sample_col_candidates <- setdiff(sample_col_candidates, replicate_col_candidates)

  if (length(sample_col_candidates) == 0 || length(replicate_col_candidates) == 0) {
    stop("❌ Could not detect 'sample_col' or 'replicate_col'. Check metadata column names.")
  }

  sample_col    <- sample_col_candidates[1]
  replicate_col <- replicate_col_candidates[1]

  message("🧠 Auto-detected sample_col: ", sample_col)
  message("🧠 Auto-detected replicate_col: ", replicate_col)

  # --- Step 4: Add interaction column SampleRep ---
  sample_data(physeq_filtered)$SampleRep <- interaction(
    sample_data(physeq_filtered)[[sample_col]],
    sample_data(physeq_filtered)[[replicate_col]]
  )

  # --- Step 5: Get updated metadata with SampleRep ---
  meta_df <- as.data.frame(sample_data(physeq_filtered), stringsAsFactors = FALSE)
  meta_df$SampleRep <- rownames(meta_df)

  # --- Step 6: Remove taxa detected in ≤5 samples ---
  physeq_filtered <- filter_taxa(physeq_filtered, function(x) sum(x > 0) > 5, prune = TRUE)

  # --- Step 7: Extract OTU table and convert to long format ---
  otu_mat <- as(otu_table(physeq_filtered), "matrix")
  if (taxa_are_rows(physeq_filtered)) {
    otu_mat <- t(otu_mat)
  }

  otu_long <- as.data.frame(otu_mat) %>%
    rownames_to_column("SampleRep") %>%
    pivot_longer(-SampleRep, names_to = "OTU", values_to = "y")

  # --- Step 8: Merge metadata with long OTU table ---
  long_df <- left_join(otu_long, meta_df, by = "SampleRep") %>%
    mutate(
      OTU = factor(OTU),
      y = as.integer(y),
      Site = .data[[site_col]],
      Sample = .data[[sample_col]],
      Replicate = .data[[replicate_col]]
    )

  # --- Return ---
  return(list(
    physeq_filtered = physeq_filtered,
    long_df = long_df,
    sample_col = sample_col,
    replicate_col = replicate_col
  ))
}
