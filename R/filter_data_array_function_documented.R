#' Filter Data Array Function
#'
#' Filters a 3D data array (Species x Sites x Replicates) by removing species and sites with zero sums,
#' and species with total counts less than a specified threshold.
#'
#' @param data_array A 3D array (Species x Sites x Replicates).
#' @param min_species_sum Minimum total count threshold to keep a species. Default is 30.
#' @param save_path Path to save the filtered data array as an .Rdata file. Default is "filtered_data.Rdata". Use NULL to avoid saving.
#'
#' @return A filtered 3D data array.
#'
#' @examples
#' @examples
#' # Example usage
#' # Load a phyloseq object and convert it to a 3D data array
#' load("Data/longdataexample.RDS")
#' data_array <- data_array_phyloseq("Data/longdataexample.RDS")
#' 
#' # Apply the filtering function
#' filtered_array <- filter_data_array(data_array, min_species_sum = 30, save_path = "filtered_data.Rdata")
#'
#' @export
filter_data_array <- function(data_array, min_species_sum = 30, save_path = "datanew_filtered.Rdata") {
  
  # ============================================
  # ✅ Validate Input
  # ============================================
  if (!is.array(data_array) || length(dim(data_array)) != 3) {
    stop("❌ Error: `data_array` must be a 3D array (Species x Sites x Replicates).")
  }

  message("✅ Starting data filtering process...")

  # Extract dimensions for reference
  species_count <- dim(data_array)[1]
  site_count <- dim(data_array)[2]
  replicate_count <- dim(data_array)[3]

  # Extract original names (if available)
  species_names <- dimnames(data_array)$Species
  site_names <- dimnames(data_array)$Sites
  replicate_names <- dimnames(data_array)$Replicates

  # ============================================
  # ✅ Step 1: Remove species with all zero counts
  # ============================================
  species_sums <- apply(data_array, 1, sum)
  species_to_keep <- species_sums > 0

  if (!any(species_to_keep)) stop("❌ No species left after removing all-zero species!")

  data_array <- data_array[species_to_keep, , , drop = FALSE]
  species_names_filtered <- species_names[species_to_keep]

  message("✅ Removed all-zero species. Remaining species: ", length(species_names_filtered))

  # ============================================
  # ✅ Step 2: Remove sites with all zero counts
  # ============================================
  site_sums <- apply(data_array, 2, sum)
  sites_to_keep <- site_sums > 0

  if (!any(sites_to_keep)) stop("❌ No sites left after removing all-zero sites!")

  data_array <- data_array[, sites_to_keep, , drop = FALSE]
  site_names_filtered <- site_names[sites_to_keep]

  message("✅ Removed all-zero sites. Remaining sites: ", length(site_names_filtered))

  # ============================================
  # ✅ Step 3: Remove rare species (below min_species_sum)
  # ============================================
  species_sums_after_site_filter <- apply(data_array, 1, sum)
  species_to_keep_final <- species_sums_after_site_filter >= min_species_sum

  if (!any(species_to_keep_final)) warning("⚠️ No species meet the minimum count threshold of ", min_species_sum, ". Array may be empty.")

  data_array_filtered <- data_array[species_to_keep_final, , , drop = FALSE]
  species_names_final <- species_names_filtered[species_to_keep_final]

  message("✅ Removed species with total counts <", min_species_sum, ". Remaining species: ", dim(data_array_filtered)[1])

  # ============================================
  # ✅ Preserve Dimension Names (if available)
  # ============================================
  dimnames(data_array_filtered) <- list(
    Species = species_names_final,
    Sites = site_names_filtered,
    Replicates = replicate_names
  )

  # ============================================
  # ✅ Validation Checks
  # ============================================
  message("\n✅ Performing validation checks...")

  # CHECK 1: All-zero species removed
  species_total_sums <- apply(data_array_filtered, 1, sum)
  if (any(species_total_sums == 0)) {
    warning("❌ Some species (rows) with all zeros still remain!")
  } else {
    message("✅ All-zero species (rows) successfully removed.")
  }

  # CHECK 2: All-zero sites removed
  site_total_sums <- apply(data_array_filtered, 2, sum)
  if (any(site_total_sums == 0)) {
    warning("❌ Some sites (columns) with all zeros still remain!")
  } else {
    message("✅ All-zero sites (columns) successfully removed.")
  }

  # CHECK 3: Rare species with total sum < threshold removed
  if (any(species_total_sums < min_species_sum)) {
    warning("❌ Some species (rows) with total counts <", min_species_sum, " still remain!")
  } else {
    message("✅ All species (rows) with total counts <", min_species_sum, " successfully removed.")
  }

  # ============================================
  # ✅ Final Summary
  # ============================================
  message("\n✅ Final 3D data array structure (Species x Sites x Replicates): ", paste(dim(data_array_filtered), collapse = " x "))
  nonzero_count_final <- sum(data_array_filtered != 0, na.rm = TRUE)
  message("✅ Total number of non-zero entries: ", nonzero_count_final)

  # ✅ Save filtered array if path provided
  if (!is.null(save_path)) {
    tryCatch({
      save(data_array_filtered, file = save_path)
      message("✅ Filtered data saved as '", save_path, "'.")
    }, error = function(e) {
      warning("⚠️ Could not save filtered data: ", e$message)
    })
  } else {
    message("ℹ️ No save path provided. Filtered data not saved to file.")
  }

  # ✅ Return filtered array
  return(data_array_filtered)
}
