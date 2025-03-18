#' Generate Data Array from Phyloseq Object
#'
#' This function reads a phyloseq object and converts it into a 3D array (Species x Sites x Replicates).
#'
#' @param phyloseq_path Path to the .RDS file containing a phyloseq object.
#' @param verbose Logical, if TRUE (default), prints progress messages.
#'
#' @return A 3D array with dimensions (Species x Sites x Replicates).
#'
#' @examples
#' \dontrun{
#' data_array_phyloseq("~/path/to/phyloseq_object.RDS")
#' }
#'
#' @export
data_array_phyloseq <- function(phyloseq_path, verbose = TRUE) {
  # Load phyloseq object using base R readRDS
  physeq <- readRDS(phyloseq_path)
  if (verbose) message("✅ Phyloseq object loaded successfully.")

  # Extract OTU Table (Species x Samples)
  otu_mat <- as.matrix(phyloseq::otu_table(physeq))
  otu_mat <- t(otu_mat)  # Transpose to ensure samples as rows, OTUs/species as columns

  # Extract Sample Metadata
  sample_data_df <- as.data.frame(phyloseq::sample_data(physeq))
  sample_data_df <- sample_data_df[sample_data_df$sampletype == "biologicalsample", ]
  if (verbose) message("✅ Biological samples retained: ", nrow(sample_data_df), " samples.")

  # Rename Columns to 'Sites' and 'Replicates'
  colnames(sample_data_df)[colnames(sample_data_df) == "sample"] <- "Sites"
  colnames(sample_data_df)[colnames(sample_data_df) == "samplerep"] <- "Replicates"

  # Extract real replicate name (assuming format "s116_r1")
  sample_data_df$Replicates <- sub(".*_", "", rownames(sample_data_df))  # Example: "r1"

  # Check for matching samples
  if (!all(rownames(sample_data_df) %in% rownames(otu_mat))) {
    stop("❌ Some samples in sample_data_df do not match samples in otu_mat!")
  } else if (verbose) {
    message("✅ All samples matched correctly between metadata and OTU table!")
  }

  # Convert Sites to numeric IDs
  sample_data_df$Sites <- as.factor(sample_data_df$Sites)
  sample_data_df$Site_ID <- as.numeric(sample_data_df$Sites)

  # Extract Names
  species_names <- colnames(otu_mat)             # OTUs
  site_names <- unique(sample_data_df$Site_ID)   # Sites
  replicate_names <- unique(sample_data_df$Replicates)  # Replicates

  if (verbose) {
    message("✅ Number of Species (OTUs): ", length(species_names))
    message("✅ Number of Sites: ", length(site_names))
    message("✅ Number of Replicates: ", length(replicate_names))  # Should be 4
  }

  # Create Empty 3D Array (Species x Sites x Replicates)
  data_array <- array(0, dim = c(length(species_names), length(site_names), length(replicate_names)),
                      dimnames = list(Species = species_names, Sites = site_names, Replicates = replicate_names))

  # Populate 3D Array
  for (i in seq_len(nrow(sample_data_df))) {
    site <- as.character(sample_data_df$Site_ID[i])
    rep <- as.character(sample_data_df$Replicates[i])
    sample_name <- rownames(sample_data_df)[i]

    if (site %in% site_names && rep %in% replicate_names && sample_name %in% rownames(otu_mat)) {
      data_array[, site, rep] <- otu_mat[sample_name, ]
    }
  }

  # Final Check of Data Structure
  if (verbose) {
    message("✅ Final 3D data array structure:")
    print(str(data_array))  # Print structure only if verbose
  }

  # Count non-zero entries
  nonzero_count <- sum(data_array != 0, na.rm = TRUE)
  if (nonzero_count == 0) {
    warning("⚠️ All entries in data_array are ZERO.")
  } else if (verbose) {
    message("✅ The data contains ", nonzero_count, " nonzero values.")
  }

  # Return the data array
  return(data_array)
}


