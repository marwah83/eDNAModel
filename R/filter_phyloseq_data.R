#' Filter and Clean a Phyloseq Object Based on Abundance Thresholds
#'
#' This function filters a `phyloseq` object by removing all-zero taxa and samples, and excluding rare taxa based on a minimum abundance threshold. It performs validation checks and optionally saves the filtered object.
#'
#' @param phyloseq_obj A valid `phyloseq` object containing OTU (or ASV) tables, sample metadata, and (optionally) taxonomy.
#' @param min_species_sum Integer. Minimum total count threshold across all samples for a taxon (species/OTU) to be retained. Default is 30.
#' @param save_path Optional. File path (as character) to save the filtered object as an `.RDS` file. Default: `"phyloseq_filtered.RDS"`. Set to `NULL` to skip saving.
#'
#' @return A filtered `phyloseq` object with:
#' \itemize{
#'   \item All-zero species removed
#'   \item All-zero samples removed
#'   \item Species with total counts below `min_species_sum` removed
#' }
#'
#' @details
#' The function performs several steps:
#' \enumerate{
#'   \item Extracts and reorients the OTU matrix if needed.
#'   \item Removes all-zero taxa (columns) and all-zero samples (rows).
#'   \item Filters out taxa with total abundance < `min_species_sum`.
#'   \item Applies pruning to the original `phyloseq` object to maintain consistency.
#'   \item Performs internal checks to ensure all-zero entries are removed and threshold filtering was applied correctly.
#'   \item Optionally saves the result as an `.RDS` file for reproducibility.
#' }
#'
#' Warning messages will alert the user if any unexpected taxa or samples remain after filtering.
#'
#' @section Example:
#' \dontrun{
#' # Load your phyloseq object
#' data("GlobalPatterns")
#' 
#' # Filter rare species and save the result
#' filtered_ps <- filter_phyloseq_data(GlobalPatterns, min_species_sum = 50)
#' 
#' # Load later:
#' restored_ps <- readRDS("phyloseq_filtered.RDS")
#' }
#'
#' @seealso \code{\link[phyloseq]{prune_taxa}}, \code{\link[phyloseq]{prune_samples}}, \code{\link{prepare_long_data}}
#'
#' @importFrom phyloseq prune_taxa prune_samples otu_table taxa_are_rows nsamples ntaxa
#' @importFrom methods is
#' @importFrom utils saveRDS
#' @export
filter_phyloseq_data <- function(phyloseq_obj, min_species_sum = 30, save_path = "phyloseq_filtered.RDS") {

  # === Validate input ===
  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("Input must be a valid `phyloseq` object.")
  }

  message("Starting phyloseq data filtering process...")

  # === Step 1: Extract OTU matrix ===
  otu_mat <- if (phyloseq::taxa_are_rows(phyloseq_obj)) {
    t(phyloseq::otu_table(phyloseq_obj))
  } else {
    phyloseq::otu_table(phyloseq_obj)
  }
  otu_mat <- as.matrix(otu_mat)

  # === Step 2: Remove all-zero species ===
  species_to_keep <- colSums(otu_mat) > 0
  if (!any(species_to_keep)) stop("No species left after removing all-zero species.")
  otu_mat <- otu_mat[, species_to_keep, drop = FALSE]

  # === Step 3: Remove all-zero samples ===
  samples_to_keep <- rowSums(otu_mat) > 0
  if (!any(samples_to_keep)) stop("No samples left after removing all-zero samples.")
  otu_mat <- otu_mat[samples_to_keep, , drop = FALSE]
  message(" Removed all-zero samples. Remaining: ", nrow(otu_mat))

  # === Step 4: Filter rare species ===
  species_to_keep_final <- colSums(otu_mat) >= min_species_sum
if (!any(species_to_keep_final)) {
  stop(" No species meet the min_species_sum of ", min_species_sum)
}
otu_mat <- otu_mat[, species_to_keep_final, drop = FALSE]

  # === Step 5: Prune original phyloseq object ===
  physeq_filtered <- prune_taxa(colnames(otu_mat), phyloseq_obj)
  physeq_filtered <- prune_samples(rownames(otu_mat), physeq_filtered)

  # === Step 6: Validation checks ===
  message("\n  Running validation checks...")

  if (any(colSums(otu_table(physeq_filtered)) == 0)) {
    warning("Some species still have all zeros.")
  } else {
    message(" All-zero species removed.")
  }

  if (any(rowSums(otu_table(physeq_filtered)) == 0)) {
    warning(" Some samples still have all zeros.")
  } else {
    message(" All-zero samples removed.")
  }

  if (any(colSums(otu_table(physeq_filtered)) < min_species_sum)) {
    warning(" Some species with total counts <", min_species_sum, " still remain.")
  } else {
    message(" All species below threshold successfully removed.")
  }

  # === Step 7: Save object ===
  if (!is.null(save_path)) {
    tryCatch({
      saveRDS(physeq_filtered, file = save_path)
      message("Filtered phyloseq object saved as: ", save_path)
    }, error = function(e) {
      warning(" Failed to save filtered object: ", e$message)
    })
  }

  # === Final Summary ===
  message("\n Final filtered phyloseq object:")
  message("   Samples: ", nsamples(physeq_filtered))
  message("   Species: ", ntaxa(physeq_filtered))
  message("   Non-zero entries: ", sum(otu_table(physeq_filtered) > 0))

  return(physeq_filtered)
}
