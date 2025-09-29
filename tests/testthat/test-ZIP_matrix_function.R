# ✅ Test create_zero_inflation_matrix function
test_that("create_zero_inflation_matrix returns correct dimensions", {
  species <- 10
  sites <- 5
  ZIP <- 0.5
  matrix_result <- create_zero_inflation_matrix(species, sites, ZIP)
  expect_equal(dim(matrix_result), c(species, sites))
  expect_true(all(matrix_result == ZIP))
})



