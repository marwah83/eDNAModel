#' Prepare Long Format Data from a Phyloseq Object
#'
#' Converts a phyloseq object into a tidy long-format data frame.
#' No taxon filtering is performed.
#'
#' @param physeq_obj A phyloseq object containing an OTU table
#'   and sample metadata.
#' @param site_col Character. Name of the site-level variable
#'   in sample_data.
#' @param nested_cols Optional character vector of metadata
#'   columns used to construct a nested grouping identifier.
#'
#' @return A list containing:
#' \describe{
#'   \item{physeq}{The phyloseq object with SampleRep added.}
#'   \item{long_df}{Long-format data with SampleRep, OTU,
#'   y, Site, and sample metadata.}
#' }
#'
#' @importFrom phyloseq sample_data sample_names
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join mutate
#' @importFrom rlang .data
#'
#' @export
prepare_long_data <- function(
    physeq_obj,
    site_col,
    nested_cols = NULL
) {

  # ----------------------------------------------------------
  # 1. Extract sample metadata
  # ----------------------------------------------------------

  sample_meta <- as.data.frame(
    phyloseq::sample_data(physeq_obj),
    stringsAsFactors = FALSE
  )

  if (!(site_col %in% names(sample_meta))) {
    stop("site_col not found in sample_data.")
  }

  # Preserve original sample identifiers
  sample_meta$SampleID <- rownames(sample_meta)

  # ----------------------------------------------------------
  # 2. Construct SampleRep
  # ----------------------------------------------------------

  if (!is.null(nested_cols)) {

    if (!all(nested_cols %in% names(sample_meta))) {
      stop("Some nested_cols not found in sample_data.")
    }

    sample_meta$SampleRep <- as.character(
      do.call(
        interaction,
        c(
          sample_meta[, nested_cols, drop = FALSE],
          list(drop = TRUE, sep = "_")
        )
      )
    )

  } else {

    sample_meta$SampleRep <- sample_meta$SampleID

  }

  # ----------------------------------------------------------
  # 3. Update phyloseq metadata
  # ----------------------------------------------------------

  phyloseq::sample_data(physeq_obj)$SampleRep <-
    sample_meta$SampleRep

  # ----------------------------------------------------------
  # 4. Extract OTU matrix
  # ----------------------------------------------------------

  otu_mat <- as(
    phyloseq::otu_table(physeq_obj),
    "matrix"
  )

  if (phyloseq::taxa_are_rows(physeq_obj)) {
    otu_mat <- t(otu_mat)
  }

  # ----------------------------------------------------------
  # 5. Convert OTU matrix to long format
  # ----------------------------------------------------------

  otu_long <- as.data.frame(otu_mat) |>
    tibble::rownames_to_column("SampleID") |>
    tidyr::pivot_longer(
      cols = -tidyselect::all_of("SampleID"),
      names_to = "OTU",
      values_to = "y"
    )

  # ----------------------------------------------------------
  # 6. Merge with sample metadata
  # ----------------------------------------------------------

  long_df <- dplyr::left_join(
    otu_long,
    sample_meta,
    by = "SampleID"
  ) |>
    dplyr::mutate(
      OTU = factor(.data$OTU),
      y = as.numeric(.data$y),
      Site = .data[[site_col]]
    )

  # ----------------------------------------------------------
  # 7. Return results
  # ----------------------------------------------------------

  return(
    list(
      physeq = physeq_obj,
      long_df = long_df
    )
  )
}
