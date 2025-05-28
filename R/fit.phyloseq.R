fit.phyloseq <- function(phyloseq_obj,
                         a.formula = ~ 1,
                         o.formula = ~ 1,
                         linko = 1,
                         linka = 0,
                         family = 1,
                         Ntrials = matrix(0),
                         offset = NULL,
                         control = list(maxit = 10000, trace = 1),
                         ...) {
  
  # Internal helper to prepare y and X
  prepare_y_X <- function(phyloseq_obj) {
    if (phyloseq::taxa_are_rows(phyloseq_obj)) {
      otu <- t(phyloseq::otu_table(phyloseq_obj))
    } else {
      otu <- phyloseq::otu_table(phyloseq_obj)
    }
    otu_mat <- as.matrix(otu)
    
    meta <- as.data.frame(phyloseq::sample_data(phyloseq_obj))
    
    if ("sampletype" %in% colnames(meta)) {
      meta <- meta[meta$sampletype == "biologicalsample", , drop = FALSE]
    } else {
      stop("❌ 'sampletype' column not found.")
    }
    
    if (nrow(meta) == 0) stop("❌ No biological samples found.")
    
    common_samples <- intersect(rownames(meta), rownames(otu_mat))
    y <- otu_mat[common_samples, ]
    X <- meta[common_samples, , drop = FALSE]
    X <- data.frame(X)
    
    if ("location" %in% colnames(X)) {
      X$Site <- as.factor(X$location)
    } else {
      stop("❌ 'location' column not found.")
    }
    
    X$Sample <- sub("(_r[0-9]+)$", "", rownames(X))
    X$Replicate <- sub(".*_(r[0-9]+)$", "\\1", rownames(X))
    
    X$Site <- as.factor(X$Site)
    X$Sample <- as.factor(X$Sample)
    X$Replicate <- as.factor(X$Replicate)
    
    X[] <- lapply(X, function(col) if (is.factor(col)) droplevels(col) else col)
    
    return(list(y = y, X = X))
  }
  
  # Prepare data
  data_list <- prepare_y_X(phyloseq_obj)
  y <- data_list$y
  X <- data_list$X
  
  # Fit model
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
  
  return(model)
}
