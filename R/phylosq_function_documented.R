#' Generate Data Array from Phyloseq Object
#'
#' Converts a phyloseq object into a 3D data array (Species x Sites x Replicates).
#'
#' @param phyloseq_path Path to the .RDS file containing a phyloseq object.
#' @param verbose Logical, if TRUE (default), prints progress messages.
#'
#' @return A 3D array with dimensions (Species x Sites x Replicates).
#'
#' @examples
#' data_path <- system.file("extdata", "longdataexample.RDS", package = "eDNAModel")
#' data_array <- data_array_phyloseq(data_path)
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

  # ✅ Ensure required columns exist
  if (!all(c("sample", "samplerep") %in% colnames(sample_data_df))) {
    stop("❌ Sample metadata must contain 'sample' (site) and 'samplerep' (replicate) columns!")
  }

  # ✅ Rename Columns to 'Sites' and 'Replicates'
  sample_data_df <- sample_data_df[, c("sample", "samplerep"), drop = FALSE]
  colnames(sample_data_df) <- c("Sites", "Replicates")

  # ✅ Extract replicate name (assuming format "s116_r1")
  sample_data_df$Replicates <- sub(".*_", "", rownames(sample_data_df))  # Extract replicate ID

  # ✅ Convert Sites to numeric IDs
  sample_data_df$Sites <- as.factor(sample_data_df$Sites)  
  sample_data_df$Site_ID <- as.numeric(sample_data_df$Sites)

  # ✅ Create Empty 3D Array (Species x Sites x Replicates)
  species_names <- colnames(otu_mat)
  site_names <- unique(sample_data_df$Site_ID)
  replicate_names <- unique(sample_data_df$Replicates)

  data_array <- array(0, dim = c(length(species_names), length(site_names), length(replicate_names)),
                      dimnames = list(Species = species_names, Sites = as.character(site_names), Replicates = as.character(replicate_names)))

  # ✅ Populate 3D Array
  for (i in seq_len(nrow(sample_data_df))) {
    site <- as.character(sample_data_df$Site_ID[i])
    rep <- as.character(sample_data_df$Replicates[i])
    sample_name <- rownames(sample_data_df)[i]

    if (sample_name %in% rownames(otu_mat) && site %in% site_names && rep %in% replicate_names) {
      data_array[, site, rep] <- otu_mat[sample_name, ]
    }
  }

  # ✅ Return the 3D data array
  return(data_array)
}
