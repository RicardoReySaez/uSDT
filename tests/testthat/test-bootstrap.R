# test-bootstrap.R
# This script tests the bootstrap filter and its interval scales.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

test_that("the bootstrap drops only samples that failed to fit", {
  samples <- rbind(
    c(1, 0.2, 0.1, 0.3, 1, 0, 0),
    c(NA, 0.2, 0.1, 0.3, 1, 0, 0),
    c(1, 0.2, 0.1, 0.3, 0, 0, 0),
    c(1, 0.2, 0.1, 0.3, 1, 1, 0),
    c(1, 1.0, 0.1, 0.3, 1, 0, 1)
  )
  colnames(samples) <- c("diff", "rho", "intercept", "slope",
                         "converged", "singular", "boundary")
  result <- .boot_filter(samples)

  # A singular or boundary sample is a result and stays in.
  expect_identical(result$ok, c(TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(result$finite, c(TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_identical(result$converged, c(TRUE, TRUE, FALSE, TRUE, TRUE))
  expect_identical(result$singular, c(FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_identical(result$boundary, c(FALSE, FALSE, FALSE, FALSE, TRUE))
})

test_that("a correlation interval stays inside its range", {
  set.seed(4)
  v <- tanh(atanh(0.8) + stats::rnorm(2000, sd = 0.4))

  # The location intervals would leave [-1, 1] without the Fisher-z scale.
  for (type in c("perc", "basic", "norm")) {
    ci <- .boot_ci(v, 0.8, 0.95, type, link = .fisher_link)
    expect_gte(ci[1L], -1)
    expect_lte(ci[2L], 1)
  }

  # An unbounded quantity is left on its own scale.
  expect_equal(drop(.boot_ci(v, 0.8, 0.95, "perc")),
               stats::quantile(v, c(0.025, 0.975), names = FALSE))

  # One call covers a whole curve, one interval per column.
  m <- cbind(v, v + 0.5, v - 0.5)
  for (type in c("perc", "basic", "norm")) {
    expect_equal(.boot_ci(m, c(0.8, 1.3, 0.3), 0.95, type),
                 t(vapply(1:3, function(j)
                   drop(.boot_ci(m[, j], c(0.8, 1.3, 0.3)[j], 0.95, type)),
                   numeric(2L))))
  }
})
