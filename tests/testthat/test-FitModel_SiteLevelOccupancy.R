test_that("FitModel_SiteLevelOccupancy", {
  library(phyloseq)
  library(glmmTMB)
  library(dplyr)  # ✅ Add this line
  # -------------------------------
  # Simulate small eDNA-like dataset
  # -------------------------------
  
  # Simulate OTU table: 3 OTUs × 6 samples
  otu_mat <- matrix(
    c(5, 2, 3, 1, 2, 1,
      1, 2, 4, 1, 3, 1,
      3, 1, 1, 2, 1, 1),
    nrow = 3, byrow = TRUE
  )
  rownames(otu_mat) <- paste0("OTU", 1:3)
  colnames(otu_mat) <- paste0("S", 1:6)
  otu_tab <- otu_table(otu_mat, taxa_are_rows = TRUE)
  
  # Simulate sample metadata: 6 samples across 2 treatments
  sample_df <- data.frame(
    sampletype = rep("biologicalsample", 6),
    location   = rep(c("Loc1", "Loc2"), each = 3),
    treatment  = c("control", "control", "control", "rats", "rats", "rats"),
    Replicate  = rep(1:3, 2),
    row.names  = paste0("S", 1:6)
  )
  sample_tab <- sample_data(sample_df)
  
  # Create phyloseq object
  physeq <- phyloseq(otu_tab, sample_tab)
  
  # ----------------------------------------
  # Run the occupancy model (with no exclusion)
  # ----------------------------------------
  
  result <- suppressWarnings(
    FitModel_SiteLevelOccupancy(
      phyloseq = physeq,
      poisson_rhs = quote((1 | OTU) + (1 | Site/OTU) + (1 | Sample/OTU) + (1 | Replicate/OTU) + treatment * OTU),
      binomial_rhs = quote((1 | OTU) + (1 | Site/OTU)),
      sampletype_keep = "biologicalsample",
      min_species_sum = 1,
      abundance_threshold = 1,
      treatment_exclude = NULL,
      n_iter = 5,
      burn_in = 2
    )
  )
  
  # -----------------------------
  # Check model output structure
  # -----------------------------
  
  expect_type(result, "list")
  expect_named(result, c("summary", "psi_list", "lambda_list", "p_detect_list", "reduced_data"))
  
  # Check summary content
  expect_s3_class(result$summary, "data.frame")
  expect_gt(nrow(result$summary), 0)
  expect_true(all(c("psi_mean", "lambda_mean", "p_detect_mean") %in% names(result$summary)))
})
