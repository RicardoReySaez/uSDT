# test-bootstrap.R
# This script tests the bootstrap arguments, filter, interval scales,
# hypothesis table and refitting loop.
# Author: Ricardo Rey-Sáez
# Last modified: 07-10-2026

# This function fits a small model with interior estimates.
boot_model <- function() {
  set.seed(21)
  trials <- make_trials(n_subj = 30L, n_trials = 80L)
  d <- usdt_data_long(trials, task_col = "task",
                      task_levels      = c(direct = "D", indirect = "I"),
                      subject_col      = "subj",
                      condition_col    = "cond",
                      condition_levels = c(signal = 1, noise = 0),
                      response_col     = "response",
                      response_levels  = c(signal = 1, noise = 0))
  hsdt(d)
}

# This function replaces .boot_mer(), the package's call to lme4::bootMer(),
# for the rest of the calling test. The fake lives in the uSDT namespace, so
# lme4 itself is never changed. Instead of refitting the model, it returns
# the statistic of the original fit with some noise in the difference.
# `fail` and `singular` mark replicates as non-converged or singular, counted
# across all calls, and `error_at` makes that call fail. The returned
# environment records the size of each batch and the arguments of the last
# call.
mock_boot_mer <- function(fail = integer(0), singular = integer(0),
                          error_at = 0L, env = parent.frame()) {
  calls <- new.env()
  calls$sizes <- integer(0)
  local_mocked_bindings(
    .boot_mer = function(x, FUN, nsim, ...) {
      done <- sum(calls$sizes)
      calls$sizes <- c(calls$sizes, nsim)
      calls$args <- list(...)
      if (length(calls$sizes) == error_at) stop("the refit ran out of memory")

      t0 <- FUN(x)
      t <- matrix(t0, nrow = nsim, ncol = length(t0), byrow = TRUE,
                  dimnames = list(NULL, names(t0)))
      t[, "diff"] <- t[, "diff"] + stats::rnorm(nsim, sd = 0.1)
      index <- done + seq_len(nsim)
      t[, "converged"] <- as.numeric(!index %in% fail)
      t[, "singular"]  <- as.numeric(index %in% singular)
      t[, "boundary"]  <- 0
      structure(list(t = t), boot.all.msgs = list())
    },
    .env = env
  )
  calls
}

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

test_that("the bootstrap refits in batches until it has enough replicates", {
  m <- boot_model()
  calls <- mock_boot_mer()
  expect_no_warning(b <- usdt_boot(m, nsim = 500))

  # Each replicate draws new random effects from the fitted model.
  expect_identical(calls$sizes, rep(100L, 5L))
  expect_identical(calls$args$type, "parametric")
  expect_false(calls$args$use.u)

  expect_true(b$boot$complete)
  expect_identical(b$boot$attempted, 500L)
  expect_identical(b$boot$usable, 500L)
  expect_identical(dim(b$boot$t), c(500L, 7L))
  expect_identical(dim(b$boot$subjects$t), c(500L, 30L, 4L))

  # The point estimates stay and the intervals come from the replicates.
  expect_identical(b$tests$estimate, m$tests$estimate)
  expect_identical(b$tests$ci_method, rep("bootstrap (perc)", 4L))
})

test_that("a seed makes the bootstrap reproducible", {
  m <- boot_model()
  mock_boot_mer()
  b1 <- usdt_boot(m, nsim = 500, seed = 1)
  b2 <- usdt_boot(m, nsim = 500, seed = 1)
  b3 <- usdt_boot(m, nsim = 500, seed = 2)

  expect_identical(b1$boot$t, b2$boot$t)
  expect_false(identical(b1$boot$t, b3$boot$t))
})

test_that("failed refits are replaced and singular ones are kept", {
  m <- boot_model()
  calls <- mock_boot_mer(fail = 1:50, singular = 101:110)
  expect_no_warning(b <- usdt_boot(m, nsim = 500))

  # The last batch asks only for the replicates still missing.
  expect_identical(calls$sizes, c(rep(100L, 5L), 50L))
  expect_identical(b$boot$attempted, 550L)
  expect_identical(b$boot$usable, 500L)
  expect_identical(b$boot$failures, c(non_finite = 0L, non_converged = 50L))
  expect_identical(b$boot$retained, c(singular = 10L, boundary = 0L))
  expect_identical(sum(b$boot$ok), 500L)
})

test_that("the bootstrap warns when more than 20% of refits fail", {
  m <- boot_model()
  mock_boot_mer(fail = 1:200)
  expect_warning(b <- usdt_boot(m, nsim = 500), "more than 20% were discarded")

  # The run still completes, so the summaries use the replicates.
  expect_true(b$boot$complete)
  expect_identical(b$boot$attempted, 700L)
  expect_identical(b$tests$ci_method, rep("bootstrap (perc)", 4L))
})

test_that("too few usable replicates keep the original results", {
  m <- boot_model()
  mock_boot_mer(fail = seq(2L, 600L, by = 2L))
  expect_warning(b <- usdt_boot(m, nsim = 500, max_attempts = 600),
                 "original hypothesis results remain")

  expect_false(b$boot$complete)
  expect_false(b$boot$summary_available)
  expect_identical(b$boot$attempted, 600L)
  expect_identical(b$boot$usable, 300L)
  expect_identical(b$tests, m$tests)
})

test_that("an incomplete run with enough replicates still summarises", {
  m <- boot_model()
  mock_boot_mer(fail = 1:50)
  expect_warning(b <- usdt_boot(m, nsim = 600, max_attempts = 600),
                 "summaries use the available replicates")

  expect_false(b$boot$complete)
  expect_true(b$boot$summary_available)
  expect_identical(b$boot$usable, 550L)
  expect_identical(b$tests$ci_method, rep("bootstrap (perc)", 4L))
})

test_that("a failed batch stops the bootstrap and keeps its message", {
  m <- boot_model()

  # A failure in the first batch leaves no replicates at all.
  mock_boot_mer(error_at = 1L)
  expect_warning(b <- usdt_boot(m, nsim = 500),
                 "The last batch failed: the refit ran out of memory")
  expect_identical(b$boot$error, "the refit ran out of memory")
  expect_identical(b$boot$attempted, 0L)
  expect_identical(dim(b$boot$t), c(0L, 7L))
  expect_identical(b$tests, m$tests)

  # A later failure keeps the replicates of the batches before it.
  mock_boot_mer(error_at = 3L)
  expect_warning(b <- usdt_boot(m, nsim = 500), "The last batch failed")
  expect_identical(b$boot$attempted, 200L)
  expect_identical(b$boot$usable, 200L)
})

test_that("the progress bar counts the usable replicates", {
  m <- boot_model()
  mock_boot_mer()
  expect_output(usdt_boot(m, nsim = 500, progress = TRUE), "100%")
})

test_that("several cores share one cluster that is stopped at the end", {
  skip_on_cran()
  m <- boot_model()
  calls <- mock_boot_mer()
  usdt_boot(m, nsim = 500, ncores = 2)

  expect_identical(calls$args$parallel, "snow")
  expect_identical(calls$args$ncpus, 2L)
  expect_s3_class(calls$args$cl, "cluster")
  expect_error(parallel::clusterEvalQ(calls$args$cl, 1))
})
