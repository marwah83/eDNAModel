#' Subset phyloseq object by site
#'
#' @param physeq A phyloseq object
#' @param site_col Column name in sample_data defining sites
#' @param keep_sites Vector of site names to retain
#'
#' @return A subsetted phyloseq object
#' @export
subset_phyloseq_by_sites <- function(physeq, site_col, keep_sites) {

  meta <- data.frame(phyloseq::sample_data(physeq), check.names = FALSE)

  if (!(site_col %in% names(meta))) {
    stop("site_col not found in sample_data.")
  }

  keep_samples <- rownames(meta)[meta[[site_col]] %in% keep_sites]

  phyloseq::prune_samples(keep_samples, physeq)
}
