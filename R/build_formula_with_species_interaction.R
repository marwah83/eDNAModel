#' Build model formula with species interaction (e.g., i)
#'
#' @param rhs A quoted right-hand side of the formula, e.g., quote((1 | Site) + treatment + offset(...))
#' @param response The response variable (default: "y" for Poisson, "z_sim" for Binomial)
#' @param species_var The species variable (default: "i")
#'
#' @return A full formula object: e.g., y ~ (1 | Site) * i + treatment * i + offset(log_reads)
build_formula_with_species_interaction <- function(rhs, response = "y", species_var = "i") {
  # Convert RHS to string
  rhs_str <- deparse(rhs)
  rhs_str <- paste(rhs_str, collapse = "")  # Flatten multiline
  
  # Split by "+" to get terms (safe for offsets now)
  terms <- strsplit(rhs_str, "\\s*\\+\\s*")[[1]]
  terms <- trimws(terms)
  
  # Separate offset terms from others
  offset_terms <- terms[grepl("^offset\\(", terms)]
  non_offset_terms <- terms[!grepl("^offset\\(", terms)]
  
  # Wrap and multiply non-offset terms by species_var
  interaction_terms <- paste0("(", non_offset_terms, ") * ", species_var)
  
  # Combine all terms (interaction + offset)
  all_terms <- c(interaction_terms, offset_terms)
  
  # Build the full formula string
  formula_str <- paste(response, "~", paste(all_terms, collapse = " + "))
  
  # Convert to formula
  as.formula(formula_str)
}
