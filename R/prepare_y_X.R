prepare_y_X <- function(phyloseq_obj) {
  # Ensure OTUs (species) are columns
  if (phyloseq::taxa_are_rows(phyloseq_obj)) {
    otu <- t(phyloseq::otu_table(phyloseq_obj))
  } else {
    otu <- phyloseq::otu_table(phyloseq_obj)
  }
  otu_mat <- as.matrix(otu)
  
  # Extract sample metadata
  meta <- as.data.frame(phyloseq::sample_data(phyloseq_obj))
  
  # Filter to biological samples only
  if ("sampletype" %in% colnames(meta)) {
    meta <- meta[meta$sampletype == "biologicalsample", , drop = FALSE]
  } else {
    stop("❌ 'sampletype' column not found in metadata.")
  }
  
  if (nrow(meta) == 0) stop("❌ No biological samples found in metadata.")
  
  # Match OTU matrix and metadata by sample names
  common_samples <- intersect(rownames(meta), rownames(otu_mat))
  y <- otu_mat[common_samples, ]
  X <- meta[common_samples, , drop = FALSE]
  X <- data.frame(X)  # ensure it's a plain data frame
  
  # Create hierarchical IDs
  if ("location" %in% colnames(X)) {
    X$Site <- as.factor(X$location)
  } else {
    stop("❌ 'location' column not found in metadata.")
  }
  
  # Derive Sample and Replicate from rownames
  X$Sample <- sub("(_r[0-9]+)$", "", rownames(X))
  X$Replicate <- sub(".*_(r[0-9]+)$", "\\1", rownames(X))
  
  # Coerce hierarchical columns to factor
  X$Sample <- as.factor(X$Sample)
  X$Replicate <- as.factor(X$Replicate)
  X$Site <- as.factor(X$Site)
  
  # Drop unused factor levels
  X[] <- lapply(X, function(col) if (is.factor(col)) droplevels(col) else col)
  
  return(list(y = y, X = X))
}
