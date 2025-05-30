test_that("fit.phyloseq runs correctly on example data", {
  skip_on_cran()

  # Load example phyloseq object
  physeq_path <- system.file("extdata", "longdataexample.RDS", package = "eDNAModel")
  physeq <- readRDS(physeq_path)

  expect_s4_class(physeq, "phyloseq")
  expect_true(!is.null(otu_table(physeq)))

  # Fit the model with default formulas
  model <- fit.phyloseq(physeq)

  # Check structure of returned object
  expect_s3_class(model, "eDNAModel")
  expect_true("occupancy_probability" %in% names(model))
  expect_true("detection_probability" %in% names(model))
  expect_true(inherits(model$TMBobj, "ADFun"))
})

test_that("fit.phyloseq throws error with invalid phyloseq input", {
  expect_error(fit.phyloseq("not_a_phyloseq_object"), "phyloseq_obj must be a phyloseq object")
})
