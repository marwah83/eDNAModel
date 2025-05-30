#' Prepare OTU Table and Metadata from a phyloseq Object
#'
#' This function extracts the OTU table and sample metadata from a `phyloseq` object,
#' filters for biological samples, and constructs a metadata frame including
#' `Site`, `Sample`, and `Replicate` columns needed for modeling.
#'
#' @param phyloseq_obj A `phyloseq` object containing OTU data and metadata.
#' @param verbose Logical; whether to print progress messages. Default is `TRUE`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{y}{A matrix of OTU counts (replicates × species).}
#'   \item{X}{A data frame of covariates including Site, Sample, and Replicate.}
#' }
#'
#' @examples
#' \dontrun{
#' physeq <- readRDS("path/to/your/phyloseq_object.rds")
#' data_list <- prepare_y_X(physeq)
#' y <- data_list$y
#' X <- data_list$X
#' }
#'
#' @export
prepare_y_X <- function(phyloseq_obj, verbose = TRUE) {
  if (verbose) message("🔧 Extracting OTU table and metadata...")

  # Ensure taxa are columns (species)
  if (phyloseq::taxa_are_rows(phyloseq_obj)) {
    otu <- t(phyloseq::otu_table(phyloseq_obj))
  } else {
    otu <- phyloseq::otu_table(phyloseq_obj)
  }

  otu_mat <- as.matrix(otu)

  # Extract and filter metadata
  meta <- as.data.frame(phyloseq::sample_data(phyloseq_obj))

  if (!"sampletype" %in% colnames(meta)) stop("❌ 'sampletype' column not found in sample_data.")
  meta <- meta[meta$sampletype == "biologicalsample", , drop = FALSE]
  if (nrow(meta) == 0) stop("❌ No biological samples found.")

  # Match OTU and metadata
  common_samples <- intersect(rownames(meta), rownames(otu_mat))
  if (length(common_samples) == 0) stop("❌ No overlapping samples between OTU table and metadata.")

  y <- otu_mat[common_samples, ]
  X <- meta[common_samples, , drop = FALSE]

  # Create Site, Replicate, and Sample columns
  if (!"location" %in% colnames(X)) stop("❌ 'location' column not found in metadata.")
  X$Site <- factor(X$location)
  X$Replicate <- sub(".*_(r[0-9]+)$", "\\1", rownames(X))
  X$Sample <- sub("(_r[0-9]+)$", "", rownames(X))

  # Clean factors
  X[] <- lapply(X, function(col) if (is.factor(col)) droplevels(col) else col)

  if (verbose) {
    message("✅ y dimensions: ", paste(dim(y), collapse = " × "))
    message("✅ X dimensions: ", paste(dim(X), collapse = " × "))
  }

  return(list(y = y, X = X))
}
