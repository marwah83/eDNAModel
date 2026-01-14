#' Filter and Clean a Phyloseq Object
#'
#' This function filters a `phyloseq` object by removing:
#' - species (taxa) with all-zero counts,
#' - samples with all-zero counts,
#' - species whose total counts are below a specified threshold.
#' The resulting filtered object is optionally saved as an `.RDS` file.
#'
#' @param phyloseq_obj A `phyloseq` object to be filtered.
#' @param min_species_sum Minimum total abundance required for a species to be retained. Default is 30.
#' @param save_path Optional path to save the filtered `phyloseq` object as an `.RDS` file. Default is `"phyloseq_filtered.RDS"`. Use `NULL` to skip saving.
#'
#' @return A filtered `phyloseq` object with rare species and empty samples removed.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Removes OTUs (species) with zero total counts across all samples.
#'   \item Removes samples with zero total counts across all OTUs.
#'   \item Removes rare species whose total abundance is below `min_species_sum`.
#'   \item Prunes the original `phyloseq` object accordingly.
#'   \item Optionally saves the cleaned object.
#' }
#'
#' Messages are printed to summarize what was filtered and whether the output was saved.
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#' filtered <- filter_phyloseq_data(GlobalPatterns, min_species_sum = 50)
#' }
#'
#' @importFrom phyloseq taxa_are_rows otu_table prune_taxa prune_samples nsamples ntaxa
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
