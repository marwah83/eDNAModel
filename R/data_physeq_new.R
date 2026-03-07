#' Example phyloseq dataset for eDNAModel
#'
#' A phyloseq dataset containing eDNA sequencing data used to demonstrate
#' the occupancy–detection modeling workflow implemented in the
#' `eDNAModel` package.
#'
#' The object is a list of phyloseq objects corresponding to different
#' sampling locations. Each element contains OTU tables, taxonomic
#' information, and sample metadata used for model fitting and
#' visualization.
#'
#' @format A list of phyloseq objects.
#'
#' @details
#' The dataset includes marine invasive species monitoring data from
#' Thessaloniki. It is used in the package vignette to illustrate
#' the occupancy–detection modeling workflow.
#'
#' @source Environmental DNA monitoring study (example dataset)
#'
#' @examples
#' data(physeq_new)
#' names(physeq_new)
#'
"physeq_new"
