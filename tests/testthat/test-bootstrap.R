# test-bootstrap.R
# This script tests the bootstrap arguments, filter, interval scales and
# hypothesis table.
# Author: Ricardo Rey-Sáez
# Last modified: 07-10-2026

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

test_that("the bootstrap rejects invalid arguments before refitting", {
  # The checks run before the model is touched, so an empty object of the
  # right class is enough to reach each of them.
  m <- structure(list(), class = "hsdt")

  expect_error(usdt_boot(list()), "must come from hsdt()", fixed = TRUE)
  expect_error(usdt_boot(m, nsim = 100), "`nsim`", fixed = TRUE)
  expect_error(usdt_boot(m, nsim = 500.5), "`nsim`", fixed = TRUE)
  expect_error(usdt_boot(m, ncores = 0), "`ncores`", fixed = TRUE)
  expect_error(usdt_boot(m, nsim = 500, max_attempts = 499),
               "`max_attempts`", fixed = TRUE)
  expect_error(usdt_boot(m, level = 1), "`level`", fixed = TRUE)
  expect_error(usdt_boot(m, seed = -1), "`seed`", fixed = TRUE)
  expect_error(usdt_boot(m, progress = NA), "`progress`", fixed = TRUE)
  expect_error(usdt_boot(m, type = "bca"), "should be one of")
})

test_that("the bootstrap fills the hypothesis table from the replicates", {
  tests <- rbind(
    .row("d'(indirect) - d'(direct)", 1, se = 0.3, statistic = 3.3,
         ci_method = "Wald"),
    .row("correlation", 0.6, ci_method = "Fisher-z"),
    .row_na("intercept", 5, "the intercept standard error is invalid"),
    .row("slope", 0.4, ci_method = "Wald")
  )

  # Ten replicates per term, placed at known distances from each estimate.
  shift <- c(-2, -1.5, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.5, 2)
  t <- cbind(diff      = 1 + shift,
             rho       = tanh(atanh(0.6) + shift / 4),
             intercept = 5 + shift / 10,
             slope     = 0.4 + shift / 2)
  result <- .boot_tests(tests, t, level = 0.9, type = "norm")

  # Each term reads its own column of replicates.
  expect_identical(result$term, tests$term)
  expect_identical(result$estimate, tests$estimate)
  expect_equal(result$se, unname(apply(t, 2L, stats::sd)))

  # Four centred replicates of the difference and the slope reach their
  # estimates, so p = (4 + 1) / (10 + 1). None reach the correlation or the
  # intercept, and their p-values stop at 1 / 11 rather than zero.
  expect_equal(result$p.value, c(5, 1, 1, 5) / 11)

  # The normal interval stays on the estimate's scale, except for the
  # correlation, which is built on the Fisher-z scale and transformed back.
  z <- stats::qnorm(0.95) * c(-1, 1)
  v <- t[, "diff"]
  expect_equal(c(result$conf.low[1L], result$conf.high[1L]),
               2 * 1 - mean(v) + z * stats::sd(v))
  v <- atanh(t[, "rho"])
  expect_equal(c(result$conf.low[2L], result$conf.high[2L]),
               tanh(2 * atanh(0.6) - mean(v) + z * stats::sd(v)))

  # The bootstrap replaces the Wald statistic and labels its own intervals.
  expect_true(all(is.na(result$statistic)))
  expect_identical(result$ci_method, rep("bootstrap (norm)", 4L))

  # A term the Wald test could not estimate is usable once replicates exist.
  expect_identical(result$status, rep("ok", 4L))
  expect_true(all(is.na(result$reason)))
})
