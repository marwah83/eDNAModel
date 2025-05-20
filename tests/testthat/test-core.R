library(testthat)
library(eDNAModel)

test_that("summary returns valid data.frame", {
  array <- array(c(
    1, 2, 1, 2, 1, 2, 1, 2,
    3, 3, 3, 3, 3, 3, 3, 3,
    0, 0, 0, 0, 0, 0, 0, 0
  ), dim = c(3, 4, 2),
  dimnames = list(
    Species = paste0("OTU", 1:3),
    Sites = paste0("S", 1:4),
    Replicates = paste0("r", 1:2)
  ))

  X <- data.frame(
    Site = factor(rep(paste0("S", 1:4), each = 2), levels = paste0("S", 1:4)),
    Replicate = factor(rep(paste0("r", 1:2), 4), levels = paste0("r", 1:2))
  )

  model <- run_full_TMB(array, X, ~ Site, ~ Site, linko = 1, linka = 0, family = 1)
  smry <- summary(model)

  expect_s3_class(smry, "data.frame")
  expect_true(all(c("parameter", "estimate", "sd") %in% colnames(smry)))
})


test_that("fitted_TMB returns list of matrices", {
  array <- array(1:24, dim = c(3, 4, 2),
                 dimnames = list(Species = paste0("OTU", 1:3),
                                 Sites = paste0("S", 1:4),
                                 Replicates = paste0("r", 1:2)))

  X <- data.frame(
    Site = factor(rep(paste0("S", 1:4), each = 2), levels = paste0("S", 1:4)),
    Replicate = factor(rep(paste0("r", 1:2), 4), levels = paste0("r", 1:2))
  )

  model <- run_full_TMB(array, X, ~ Site, ~ Site, linko = 1, linka = 0, family = 1)
  fit <- fitted_TMB(model)

  expect_type(fit, "list")
  expect_true(all(c("fitted_occupancy", "fitted_abundance") %in% names(fit)))
  expect_true(is.matrix(fit$fitted_occupancy))
  expect_true(is.matrix(fit$fitted_abundance))
})


test_that("compute_residuals_TMB returns residuals of correct shape", {
  array <- array(1:24, dim = c(3, 4, 2),
                 dimnames = list(Species = paste0("OTU", 1:3),
                                 Sites = paste0("S", 1:4),
                                 Replicates = paste0("r", 1:2)))

  X <- data.frame(
    Site = factor(rep(paste0("S", 1:4), each = 2), levels = paste0("S", 1:4)),
    Replicate = factor(rep(paste0("r", 1:2), 4), levels = paste0("r", 1:2))
  )

  model <- run_full_TMB(array, X, ~ Site, ~ Site, linko = 1, linka = 0, family = 1)
  resids <- compute_residuals_TMB(model, array)

  expect_true(all(c("occupancy_residuals", "abundance_residuals") %in% names(resids)))
  expect_true(is.matrix(resids$occupancy_residuals))
  expect_true(is.matrix(resids$abundance_residuals))
  expect_equal(dim(resids$occupancy_residuals), dim(resids$abundance_residuals))
})


test_that("predict_TMB returns prediction data frame", {
  array <- array(1:24, dim = c(3, 4, 2),
                 dimnames = list(Species = paste0("OTU", 1:3),
                                 Sites = paste0("S", 1:4),
                                 Replicates = paste0("r", 1:2)))

  X <- data.frame(
    Site = factor(rep(paste0("S", 1:4), each = 2), levels = paste0("S", 1:4)),
    Replicate = factor(rep(paste0("r", 1:2), 4), levels = paste0("r", 1:2))
  )

  model <- run_full_TMB(array, X, ~ Site, ~ Site, linko = 1, linka = 0, family = 1)

  newX <- X[1:2, ]
  pred_abund <- predict_TMB(model, newX = newX, formula = ~ Site,
                            which = "abundance", type = "response")
  pred_occ <- predict_TMB(model, newX = newX, formula = ~ Site,
                          which = "occupancy", type = "response")

  expect_s3_class(pred_occ, "data.frame")
  expect_s3_class(pred_abund, "data.frame")
  expect_true(all(c("prediction", "se", "lower", "upper") %in% names(pred_occ)))
})
