#' Prepare Long Format Data from a Phyloseq Object
#'
#' Converts a `phyloseq` object into a tidy long-format data frame for modeling.
#' This function performs **data restructuring only** and does not apply any
#' filtering of taxa. All filtering should be handled explicitly in downstream
#' modeling functions such as \code{\link{FitModel}}.
#'
#' @param physeq_obj A `phyloseq` object containing OTU (or ASV) abundance table and sample metadata.
#' @param site_col Character. Column name in `sample_data` representing the sampling site variable.
#' This variable is preserved in the output for site-level modeling.
#' @param nested_cols Optional character vector. Column names in `sample_data` to be combined into a
#' nested grouping factor (`SampleRep`). Used for hierarchical or repeated-measures designs.
#' Default is `NULL`, in which case original sample names are used.
#' @param otu_col Character. Name of the OTU (taxon) column in the output long data.
#' Default is `"OTU"`.
#' @param count_col Character. Name of the abundance/count column in the output long data.
#' Default is `"y"`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{physeq}{The original `phyloseq` object (unchanged).}
#'   \item{long_df}{A long-format `data.frame` with one row per `SampleRep × OTU`.
#'   Columns include:
#'   \itemize{
#'     \item `SampleRep`: sample or nested replicate identifier
#'     \item `otu_col`: taxon identifier (user-defined name)
#'     \item `count_col`: observed counts/abundance (user-defined name)
#'     \item all original sample metadata variables (including `site_col`)
#'   }}
#' }
#'
#' @details
#' This function reshapes a `phyloseq` object into a long (tidy) format suitable
#' for hierarchical modeling frameworks such as occupancy–detection or
#' zero-inflated models.
#'
#' Key steps:
#' \itemize{
#'   \item Constructs a `SampleRep` variable:
#'     \itemize{
#'       \item If `nested_cols` are provided, they are combined using `interaction()`
#'       to define nested or repeated-measures structure.
#'       \item Otherwise, original sample names are used.
#'     }
#'   \item Converts the OTU table into long format using \code{pivot_longer()}.
#'   \item Merges OTU counts with sample metadata via `SampleRep`.
#'   \item Ensures consistent data types for modeling (factor OTUs, numeric counts).
#' }
#'
#' \strong{Important:}
#' This function does \emph{not} perform any filtering of taxa (e.g., by abundance
#' or prevalence). This design ensures transparency and allows users to control
#' filtering explicitly in modeling functions such as \code{\link{FitModel}}.
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#'
#' long_data <- prepare_long_data(
#'   physeq_obj = GlobalPatterns,
#'   site_col = "SampleType",
#'   nested_cols = c("SampleType"),
#'   otu_col = "OTU",
#'   count_col = "y"
#' )
#' }
#'
#' @seealso \code{\link{FitModel}} for applying filtering and hierarchical modeling.
#'
#' @importFrom phyloseq sample_data taxa_are_rows otu_table sample_names
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @export
prepare_long_data <- function(physeq_obj,
                              site_col,
                              nested_cols = NULL) {
  

  sample_meta <- as.data.frame(sample_data(physeq_obj), stringsAsFactors = FALSE)

  if (!(site_col %in% names(sample_meta))) {
    stop("site_col not found in sample_data.")
  }

  # -------------------------------
  # Create SampleRep
  # -------------------------------
  if (!is.null(nested_cols)) {
    if (!all(nested_cols %in% names(sample_meta))) {
      stop("Some nested_cols not found in sample_data.")
    }

    sample_data(physeq_obj)$SampleRep <-
      do.call(interaction, sample_meta[, nested_cols, drop = FALSE])
  } else {
    sample_data(physeq_obj)$SampleRep <- sample_names(physeq_obj)
  }

  # -------------------------------
  # Metadata
  # -------------------------------
  meta_df <- as.data.frame(sample_data(physeq_obj), stringsAsFactors = FALSE)
  meta_df$SampleRep <- rownames(meta_df)

  # -------------------------------
  # OTU matrix → long
  # -------------------------------
  otu_mat <- as(otu_table(physeq_obj), "matrix")

  # OTU matrix to long format
  
  if (taxa_are_rows(physeq_obj)) {
    otu_mat <- t(otu_mat)
  }

  
  otu_long <- as.data.frame(otu_mat) |>
    tibble::rownames_to_column("SampleRep") |>
    tidyr::pivot_longer(-SampleRep, names_to = "OTU", values_to = "y")

  # -------------------------------
  # Merge
  # -------------------------------
  long_df <- dplyr::left_join(otu_long, meta_df, by = "SampleRep") |>
    dplyr::mutate(
      OTU = factor(OTU),
      y   = as.numeric(y),
      Site = .data[[site_col]]
    )

  return(list(
    physeq = physeq_obj,
    long_df = long_df
  ))
}
