# test-fit.phyloseq.R

test_that("fit.phyloseq() works on a small example", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("eDNAModel")

  library(phyloseq)
  library(eDNAModel)

  # Load test data
  data_path <- system.file("extdata", "longdataexample.RDS", package = "eDNAModel")
  skip_if(data_path == "", "Test RDS file not found.")

  physeq <- readRDS(data_path)

  # Run model
  expect_error({
    model <- fit.phyloseq(
      phyloseq_obj = physeq,
      a.formula = ~ 1,
      o.formula = ~ 1,
      linko = 1,
      linka = 0,
      family = 1,
      control = list(maxit = 50, trace = 0),
      verbose = FALSE
    )
  }, NA)  # Expect no error

  # Check output class
  expect_s3_class(model, "eDNAModel")

  # Basic content check
  expect_true("occupancy_probability" %in% names(model))
  expect_true("detection_probability" %in% names(model))
  expect_s3_class(model$TMBobj, "ADFun")
})
