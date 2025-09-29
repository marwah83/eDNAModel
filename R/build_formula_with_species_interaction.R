#' Build Formula with Species Interaction
#'
#' Automatically constructs a model formula where each covariate term 
#' (except a designated lower-level random effect) is interacted with 
#' the species variable. This is used in occupancy-abundance modeling 
#' to fit multi-species GLMMs.
#'
#' @param rhs The right-hand side of the formula as a `formula` object (e.g., `~ (1 | Site) + treatment`).
#' @param response A character string specifying the response variable (default: `"y"`).
#' @param species_var The name of the species variable to interact with each term (default: `"OTU"`).
#' @param lower_level A variable name (usually a random effect) that should NOT be interacted with the species (default: `"Replicate"`).
#'
#' @return A `formula` object with the specified response and species interaction terms.
#'
#' @examples
#' build_formula_with_species_interaction(
#'   rhs = ~ (1 | Site) + (1 | Sample) + (1 | Replicate) + treatment,
#'   response = "y",
#'   species_var = "OTU",
#'   lower_level = "Replicate"
#' )
#'
#' @export
build_formula_with_species_interaction <- function(rhs, response = "y", species_var = "OTU", lower_level = "Replicate") {
  # Convert RHS formula to a string
  rhs_str <- paste(deparse(rhs), collapse = "")
  
  # Split terms while preserving proper formula components
  terms <- strsplit(rhs_str, "\\s*\\+\\s*")[[1]]
  terms <- trimws(terms)
  
  # Process each term to interact with species_var, unless it's the lower-level effect
  new_terms <- vapply(terms, function(term) {
    if (grepl("^offset\\s*\\(", term)) {
      return(term)  # Keep offset terms unchanged
    } else if (grepl(lower_level, term)) {
      return(term)  # Don't interact lower-level random effects with species
    } else {
      return(paste0("(", term, ") * ", species_var))  # Interact with species
    }
  }, character(1))
  
  # Construct and return the full formula
  full_formula_str <- paste(response, "~", paste(new_terms, collapse = " + "))
  as.formula(full_formula_str)
}
