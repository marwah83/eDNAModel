#' Fit an eDNA Model from a phyloseq Object
#'
#' @param phyloseq_obj A phyloseq object with OTU and sample data
#' @param a.formula Abundance formula (default: ~1)
#' @param o.formula Occupancy formula (default: ~1)
#' @param linko Occupancy link function (1 = logit, 2 = probit)
#' @param linka Abundance link function (0 = log, 1 = logit, 2 = probit, 3 = cloglog)
#' @param family Distribution family (0 = ZIP, 1 = ZINB, 2 = Binomial)
#' @param Ntrials Number of trials for binomial models (default: matrix(0))
#' @param offset Offset matrix (same dimensions as y)
#' @param control Control list for optimizer (default: list(maxit = 10000, trace = 1))
#' @param verbose Whether to print status messages (default: TRUE)
#' @param ... Additional arguments passed to `run_full_TMB`
#'
#' @return A fitted model of class `eDNAModel`
#' @export
fit.phyloseq <- function(phyloseq_obj,
                         a.formula = ~ 1,
                         o.formula = ~ 1,
                         linko = 1,
                         linka = 0,
                         family = 1,
                         Ntrials = matrix(0),
                         offset = NULL,
                         control = list(maxit = 10000, trace = 1),
                         verbose = TRUE,
                         ...) {
  
  # Internal helper to prepare y and X
  prepare_y_X <- function(phyloseq_obj, verbose = TRUE) {
    if (verbose) message("🔧 Extracting OTU table and metadata...")

    if (phyloseq::taxa_are_rows(phyloseq_obj)) {
      otu <- t(phyloseq::otu_table(phyloseq_obj))
    } else {
      otu <- phyloseq::otu_table(phyloseq_obj)
    }

    otu_mat <- as.matrix(otu)
    meta <- as.data.frame(phyloseq::sample_data(phyloseq_obj))

    if (!"sampletype" %in% colnames(meta)) stop("❌ 'sampletype' column not found.")
    meta <- meta[meta$sampletype == "biologicalsample", , drop = FALSE]
    if (nrow(meta) == 0) stop("❌ No biological samples found.")

    common_samples <- intersect(rownames(meta), rownames(otu_mat))
    if (length(common_samples) == 0) stop("❌ No overlapping samples between OTU table and metadata.")

    y <- otu_mat[common_samples, ]
    X <- meta[common_samples, , drop = FALSE]
    X <- data.frame(X)

    if (!"location" %in% colnames(X)) stop("❌ 'location' column not found.")
    X$Site <- factor(X$location)
    X$Sample <- sub("(_r[0-9]+)$", "", rownames(X))
    X$Replicate <- sub(".*_(r[0-9]+)$", "\\1", rownames(X))

    X$Site <- as.factor(X$Site)
    X$Sample <- as.factor(X$Sample)
    X$Replicate <- as.factor(X$Replicate)
    X[] <- lapply(X, function(col) if (is.factor(col)) droplevels(col) else col)

    if (verbose) {
      message("✅ y dimensions: ", paste(dim(y), collapse = " × "))
      message("✅ X dimensions: ", paste(dim(X), collapse = " × "))
    }

    return(list(y = y, X = X))
  }

  if (verbose) message("🚀 Preparing inputs...")
  data_list <- prepare_y_X(phyloseq_obj, verbose = verbose)
  y <- data_list$y
  X <- data_list$X

  if (verbose) message("📦 Running TMB model...")
  model <- run_full_TMB(
    y = y,
    X = X,
    a.formula = a.formula,
    o.formula = o.formula,
    linko = linko,
    linka = linka,
    family = family,
    Ntrials = Ntrials,
    offset = offset,
    control = control,
    ...
  )

  if (verbose) message("✅ Model fitting complete.")
  return(model)
}
