test_that("gradient PASS/MARGINAL/FAIL thresholds behave as intended", {

  gradient_tol <- 1e-3
  gradient_marginal_factor <- 5

  classify_gradient <- function(g) {

    if (!is.finite(g)) {

      "FAIL"

    } else if (
      g <= gradient_tol
    ) {

      "PASS"

    } else if (
      g <=
        gradient_marginal_factor *
        gradient_tol
    ) {

      "MARGINAL"

    } else {

      "FAIL"
    }
  }


  # ----------------------------
  # PASS
  # ----------------------------

  expect_identical(
    classify_gradient(0),
    "PASS"
  )

  expect_identical(
    classify_gradient(0.0005),
    "PASS"
  )

  expect_identical(
    classify_gradient(0.001),
    "PASS"
  )


  # ----------------------------
  # MARGINAL
  # ----------------------------

  expect_identical(
    classify_gradient(0.0011),
    "MARGINAL"
  )

  # Reviewer's example
  expect_identical(
    classify_gradient(0.0038),
    "MARGINAL"
  )

  expect_identical(
    classify_gradient(0.005),
    "MARGINAL"
  )


  # ----------------------------
  # FAIL
  # ----------------------------

  expect_identical(
    classify_gradient(0.0051),
    "FAIL"
  )

  expect_identical(
    classify_gradient(0.018),
    "FAIL"
  )

  expect_identical(
    classify_gradient(Inf),
    "FAIL"
  )

  expect_identical(
    classify_gradient(NA_real_),
    "FAIL"
  )

})
