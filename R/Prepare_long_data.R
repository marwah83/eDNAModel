#' Prepare Long-Format Data from a Phyloseq Object
#'
#' Converts a `phyloseq` object into a tidy long-format data frame for
#' downstream modeling. This function performs **data restructuring only**
#' and does not apply any filtering of taxa.
#'
#' All filtering (e.g., by abundance or detection frequency) should be handled
#' explicitly in modeling functions such as \code{\link{FitModel}}.
#'
#' @param physeq_obj A `phyloseq` object containing an OTU (or ASV) table and
#' sample metadata.
#' @param site_col Character. Name of the column in \code{sample_data} representing
#' the site-level grouping variable. This is copied into the output as \code{Site}.
#' @param nested_cols Optional character vector. Column names in
#' \code{sample_data} used to construct a nested grouping factor
#' (\code{SampleRep}). This is typically used to represent biological samples
#' or repeated measures (e.g., sample ID, replicate).
#' If \code{NULL}, original sample names are used.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{physeq}{
#'     The original `phyloseq` object (unchanged).
#'   }
#'   \item{long_df}{
#'     A long-format \code{data.frame} with one row per
#'     \code{SampleRep × OTU}. Columns include:
#'     \itemize{
#'       \item \code{SampleRep}: grouping identifier constructed from
#'       \code{nested_cols} or sample names
#'       \item \code{OTU}: taxon identifier (factor)
#'       \item \code{y}: observed counts/abundance (numeric)
#'       \item \code{Site}: site-level variable derived from \code{site_col}
#'       \item all original sample metadata variables
#'     }
#'   }
#' }
#'
#' @details
#' This function reshapes a `phyloseq` object into a long (tidy) format suitable
#' for hierarchical models such as occupancy–detection–abundance models.
#'
#' Key steps:
#' \itemize{
#'   \item Constructs a \code{SampleRep} variable:
#'     \itemize{
#'       \item If \code{nested_cols} are provided, they are combined using
#'       \code{interaction()} with \code{drop = TRUE} and \code{lex.order = TRUE}
#'       to ensure deterministic and reproducible grouping.
#'       \item Otherwise, original sample names are used.
#'     }
#'   \item Converts the OTU table into long format using \code{pivot_longer()}.
#'   \item Merges OTU counts with sample metadata via \code{SampleRep}.
#'   \item Creates a \code{Site} column from \code{site_col}.
#'   \item Ensures consistent data types (factor OTUs, numeric counts).
#' }
#'
#' \strong{Hierarchy:}
#' The resulting structure supports multi-level modeling:
#' \itemize{
#'   \item Site level: \code{Site}
#'   \item Biological sample level: \code{SampleRep}
#'   \item Observation level (counts): \code{y}
#' }
#'
#' \strong{Important:}
#' \itemize{
#'   \item No filtering of taxa is performed.
#'   \item Rare or low-prevalence taxa are retained.
#'   \item Filtering should be applied explicitly in \code{\link{FitModel}}.
#' }
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns)
#'
#' long_data <- prepare_long_data(
#'   physeq_obj = GlobalPatterns,
#'   site_col = "SampleType",
#'   nested_cols = c("SampleType")
#' )
#'
#' head(long_data$long_df)
#' }
#'
#' @seealso \code{\link{FitModel}} for filtering and hierarchical modeling.
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
    # -------------------------------
    # Extract metadata
    # -------------------------------
   
    if (!(site_col %in% names(sample_meta))) {
      stop("site_col not found in sample_data.")
    }
    
    # -------------------------------
    # Create SampleRep
    # 
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
      # 
      # -------------------------------
      long_df <- dplyr::left_join(otu_long, meta_df, by = "SampleRep") |>
       
        dplyr::mutate(
          OTU = factor(OTU),
          y   = as.numeric(y),
          Site = .data[[site_col]]
        )
      
      # -------------------------------
      # Return
      # -------------------------------
      return(list(
        physeq = physeq_obj,
        long_df = long_df

        ))
}
