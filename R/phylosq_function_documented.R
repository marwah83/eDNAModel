#' Generate Data Array from Phyloseq Object
#'
#' This function reads a phyloseq object and converts it into a 3D array (Species x Sites x Replicates).
#'
#' @param phyloseq_path Path to the .RDS file containing a phyloseq object.
#' @param verbose Logical, if TRUE (default), prints progress messages.
#'
#' @return A 3D array with dimensions (Species x Sites x Replicates).
#'
@examples
#' \dontrun{
#' # Load a phyloseq object and convert it to a 3D data array
#' data_array <- data_array_phyloseq("~/Desktop/Diversea/longdataexample.RDS")
#'
#' # Now use this data_array with the filter function
#' filtered_array <- filter_data_array(data_array, min_species_sum = 30, save_path = "filtered_data.Rdata")
#' }
#'
#'
#' @export
data_array_phyloseq <- function(phyloseq_path, verbose = TRUE) {
  
  # ✅ Load phyloseq object
  physeq <- readRDS(phyloseq_path)
  if (verbose) message("✅ Phyloseq object loaded successfully.")

  # ✅ Extract OTU Table (Species x Samples)
  otu_mat <- as.matrix(phyloseq::otu_table(physeq))
  
  # Ensure correct orientation: Species (OTUs) as rows, Samples as columns
  if (phyloseq::taxa_are_rows(physeq)) {
    otu_mat <- t(otu_mat)
  }

  # ✅ Extract Sample Metadata
  sample_data_df <- as.data.frame(phyloseq::sample_data(physeq))

  # ✅ Keep only biological samples
  if ("sampletype" %in% colnames(sample_data_df)) {
    sample_data_df <- sample_data_df[sample_data_df$sampletype == "biologicalsample", ]
  }

  # ✅ Ensure column names match expected format
  if (!("sample" %in% colnames(sample_data_df)) || !("samplerep" %in% colnames(sample_data_df))) {
    stop("❌ Sample metadata must contain 'sample' (site) and 'samplerep' (replicate) columns!")
  }

  # ✅ Rename Columns to 'Sites' and 'Replicates'
  sample_data_df <- sample_data_df[, c("sample", "samplerep"), drop = FALSE]
  colnames(sample_data_df) <- c("Sites", "Replicates")

  # ✅ Extract replicate name (assuming format "s116_r1")
  sample_data_df$Replicates <- sub(".*_", "", rownames(sample_data_df))  # Extract replicate ID

  # ✅ Check for matching samples
  matching_samples <- intersect(rownames(sample_data_df), rownames(otu_mat))
  if (length(matching_samples) == 0) {
    stop("❌ No matching samples found between sample metadata and OTU table!")
  }
  if (verbose) message("✅ All samples matched correctly between metadata and OTU table!")

  # ✅ Convert Sites to numeric IDs
  sample_data_df$Sites <- as.factor(sample_data_df$Sites)  
  sample_data_df$Site_ID <- as.numeric(sample_data_df$Sites)

  # ✅ Extract Dimension Names
  species_names <- colnames(otu_mat)               # OTUs (species)
  site_names <- unique(sample_data_df$Site_ID)     # Unique site IDs
  replicate_names <- unique(sample_data_df$Replicates)  # Unique replicate IDs

  if (verbose) {
    message("✅ Number of Species (OTUs): ", length(species_names))
    message("✅ Number of Sites: ", length(site_names))
    message("✅ Number of Replicates: ", length(replicate_names))
  }

  # ✅ Create Empty 3D Array (Species x Sites x Replicates)
  data_array <- array(0, dim = c(length(species_names), length(site_names), length(replicate_names)),
                      dimnames = list(Species = species_names, Sites = as.character(site_names), Replicates = as.character(replicate_names)))

  # ✅ Populate 3D Array
  for (i in seq_len(nrow(sample_data_df))) {
    site <- as.character(sample_data_df$Site_ID[i])
    rep <- as.character(sample_data_df$Replicates[i])
    sample_name <- rownames(sample_data_df)[i]

    # Ensure sample exists in OTU matrix before assigning values
    if (sample_name %in% rownames(otu_mat) && site %in% site_names && rep %in% replicate_names) {
      data_array[, site, rep] <- otu_mat[sample_name, ]
    }
  }

  # ✅ Final Check of Data Structure
  if (verbose) {
    message("✅ Final 3D data array structure:")
    print(str(data_array))  # Print structure only if verbose
  }

  # ✅ Count non-zero entries
  nonzero_count <- sum(data_array != 0, na.rm = TRUE)
  if (nonzero_count == 0) {
    warning("⚠️ All entries in data_array are ZERO.")
  } else if (verbose) {
    message("✅ The data contains ", nonzero_count, " nonzero values.")
  }

  # ✅ Return the 3D data array
  return(data_array)
}
