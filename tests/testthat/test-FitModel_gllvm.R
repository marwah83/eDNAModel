test_that("FitModel_gllvm runs and returns expected structure with synthetic data (3 Sites)", {
     skip_if_not_installed("gllvm")
     skip_if_not_installed("glmmTMB")
     skip_if_not_installed("phyloseq")
     
     # Create synthetic OTU table: 3 OTUs × 6 Samples
     otu_mat <- matrix(
         c(5, 2, 3, 4, 6, 1,
           1, 4, 4, 3, 3, 2,
           3, 1, 2, 2, 5, 4),
         nrow = 3, byrow = TRUE
     )
     rownames(otu_mat) <- paste0("OTU", 1:3)
     colnames(otu_mat) <- paste0("S", 1:6)
     otu_tab <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
     
     # Sample metadata: 3 Sites, 6 Samples
     sample_df <- data.frame(
         Site           = rep(c("Loc1", "Loc2", "Loc3"), each = 2),
         Name           = paste0("Sample", 1:6),
         Samplingmonth  = factor(c("Jan", "Feb", "Jan", "Feb", "Mar", "Mar")),
         Replicate      = rep(1:2, 3),
         row.names      = paste0("S", 1:6)
     )
     sample_tab <- phyloseq::sample_data(sample_df)
     
     # Create phyloseq object
     physeq <- phyloseq::phyloseq(otu_tab, sample_tab)
     
     # Run FitModel_gllvm with low min_species_sum for testability
     out <- suppressWarnings(
         FitModel_gllvm(
             phyloseq = physeq,
             site_col = "Site",
             abundance_rhs = (1 | Samplingmonth / OTU) + (1 | OTU),
             occupancy_covars = "Samplingmonth",
             abundance_family = "poisson",
             min_species_sum = 1,
             n_iter = 3,
             burn_in = 1,
             num_lv_c=1
         )
     )
     
     # Test structure
     expect_type(out, "list")
     
     expected_components <- c(
         "summary", "psi_list", "lambda_list", "p_detect_list",
         "occupancy_models", "abundance_models", "reduced_data",
         "lv_sites", "lv_species", "mean_lv_sites", "mean_lv_species"
     )
     expect_true(all(expected_components %in% names(out)))
     expect_gt(nrow(out$summary), 0)
     expect_s3_class(out$lv_sites, "data.frame")
     expect_s3_class(out$mean_lv_species, "data.frame")
 })
