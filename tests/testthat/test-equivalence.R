# test-equivalence.R
# This script tests the equivalence, one-sided and minimal-effect tests: their
# arguments, their agreement with the reported intervals, the Fisher-z scale
# of the correlation, boundary rows and bootstrap replicates.
# Author: Ricardo Rey-Sáez
# Last modified: 27-10-2026

# This function fits a small model with interior estimates.
equivalence_model <- function() {
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

# This function gives bootstrap replicates of a known shape, in the layout
# that usdt_boot() stores. The last replicate failed to fit and is ignored.
fake_boot <- function(type = "perc") {
  shift <- c(-2, -1.5, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.5, 2, 50)
  t <- cbind(diff      = 0.2 + shift / 10,
             rho       = tanh(atanh(0.4) + shift / 4),
             intercept = 0.1 + shift / 10,
             slope     = 0.3 + shift / 10)
  list(t = t, t0 = c(diff = 0.2, rho = 0.4, intercept = 0.1, slope = 0.3),
       ok = c(rep(TRUE, 10L), FALSE), type = type, summary_available = TRUE)
}

test_that("every hypothesis must be stated, and only for frequentist fits", {
  m <- structure(list(estimation = "frequentist"), class = "hsdt")
  expect_error(usdt_equivalence(m), "no default")
  expect_error(usdt_equivalence(m, 0.2), "one or more strings")
  expect_error(usdt_equivalence(m, "diff in [-0.2, 0.2]", alpha = 0.5),
               "alpha")
  expect_error(usdt_equivalence(m, "beta = 0"), "quantities")

  bayes <- structure(list(estimation = "bayesian"), class = "hsdt")
  expect_error(usdt_equivalence(bayes, "diff > 0"), "usdt_bf()", fixed = TRUE)
  expect_error(usdt_equivalence(list(), "diff > 0"), "must come from")
})

test_that("the frequentist tests read the slope as it is written", {
  # Only a Bayes factor moves a comparison of the slope with zero to rho.
  h <- .parse_hypothesis("slope = 0")
  expect_identical(c(h$quantity, h$op), c("slope", "="))
  expect_null(h$swapped)
  expect_identical(.bf_parse("slope = 0")$quantity, "rho")
})

test_that("each test agrees with the interval it reports", {
  m <- equivalence_model()
  hypotheses <- c("diff in [-0.6, 0.6]", "diff out [-0.6, 0.6]",
                  "diff > -0.6", "diff < -0.6", "diff = -0.6",
                  "rho in [-0.3, 0.9]", "intercept in [-0.5, 0.5]",
                  "slope out [0.1, 0.2]")
  e <- usdt_equivalence(m, hypotheses)
  expect_s3_class(e, "usdt_equivalence")
  expect_identical(e$status, rep("ok", length(hypotheses)))

  # One-sided and equivalence tests report the 1 - 2 alpha interval, the
  # two-sided and minimal-effect tests the 1 - alpha one.
  expect_equal(e$conf.level, c(0.9, 0.95, 0.9, 0.9, 0.95, 0.9, 0.9, 0.95))
  ref90 <- usdt_tests(m, level = 0.9)
  ref95 <- usdt_tests(m, level = 0.95)
  expect_equal(c(e$conf.low[1L], e$conf.high[1L]),
               c(ref90$conf.low[1L], ref90$conf.high[1L]))
  expect_equal(c(e$conf.low[2L], e$conf.high[2L]),
               c(ref95$conf.low[1L], ref95$conf.high[1L]))
  expect_equal(c(e$conf.low[6L], e$conf.high[6L]),
               c(ref90$conf.low[2L], ref90$conf.high[2L]))

  # The decision rejects H0 exactly when the interval clears the bounds.
  inside <- function(i, a, b) e$conf.low[i] > a && e$conf.high[i] < b
  outside <- function(i, a, b) e$conf.high[i] < a || e$conf.low[i] > b
  clears <- c(inside(1L, -0.6, 0.6), outside(2L, -0.6, 0.6),
              e$conf.low[3L] > -0.6, e$conf.high[4L] < -0.6,
              outside(5L, -0.6, -0.6), inside(6L, -0.3, 0.9),
              inside(7L, -0.5, 0.5), outside(8L, 0.1, 0.2))
  expect_identical(e$decision == "reject H0", clears)

  # The difference is a Wald test on its own scale.
  est <- ref95$estimate[1L]
  se <- ref95$se[1L]
  expect_equal(e$p.lower[1L], stats::pnorm((-0.6 - est) / se))
  expect_equal(e$p.upper[1L], stats::pnorm((est - 0.6) / se))
  expect_equal(e$p.value[1L], max(e$p.lower[1L], e$p.upper[1L]))
  expect_equal(e$p.value[2L], min(1, 2 * min(e$p.lower[2L], e$p.upper[2L])))
  expect_equal(e$p.value[3L], stats::pnorm(-(est + 0.6) / se))
  expect_equal(e$p.value[5L], .wald(est, se, null = -0.6)$p.value)

  # The same model fitted directly with glmer gives the same tests.
  direct <- usdt_equivalence(m$fit, hypotheses[1:3])
  expect_equal(direct$p.value, e$p.value[1:3])
})

test_that("the correlation is bounded on its scale and tested on Fisher-z", {
  m <- equivalence_model()
  r <- latent_cor(m)
  se_z <- r$se / (1 - r$estimate^2)
  e <- usdt_equivalence(m, c("rho > -0.2", "rho in [-0.3, 0.9]"))

  expect_identical(e$ci_method, rep("Fisher-z", 2L))
  expect_equal(e$estimate, rep(r$estimate, 2L))
  expect_equal(e$se, rep(r$se, 2L))
  expect_equal(e$p.value[1L],
               stats::pnorm(-(atanh(r$estimate) - atanh(-0.2)) / se_z))
  expect_equal(e$p.lower[2L],
               stats::pnorm(-(atanh(r$estimate) - atanh(-0.3)) / se_z))
  expect_equal(e$p.upper[2L],
               stats::pnorm((atanh(r$estimate) - atanh(0.9)) / se_z))
  expect_true(all(e$conf.low > -1 & e$conf.high < 1))
})

test_that("a boundary fit keeps the difference and points to the bootstrap", {
  m <- equivalence_model()
  p <- .pars_of(m, "d_D", "d_I")
  p$joint_ok <- FALSE
  p$inference_reason <- "the sensitivity covariance is on the boundary"
  s <- .equiv_wald(p, bootable = TRUE)

  expect_true(s$diff$ok)
  rows <- do.call(rbind, lapply(c("rho", "intercept", "slope"), function(q) {
    .equiv_one(.parse_hypothesis(paste(q, "in [-0.3, 0.3]")), s[[q]], 0.05)
  }))
  expect_identical(rows$status, rep("not estimable", 3L))
  expect_true(all(is.na(rows$p.value) & is.na(rows$conf.low)))
  expect_match(rows$reason, "on the boundary; usdt_boot()", fixed = TRUE)

  e <- structure(rows, alpha = 0.05,
                 class = c("usdt_equivalence", "data.frame"))
  out <- utils::capture.output(print(e))
  expect_length(grep("Not estimable", out, fixed = TRUE), 3L)
})

test_that("a bootstrapped fit tests the usable replicates", {
  perc <- .equiv_boot(fake_boot("perc"))
  v <- fake_boot()$t[1:10, ]

  # Each p-value is the share of the distribution beyond the bound,
  # counted as (k + 1) / (B + 1).
  h <- .parse_hypothesis("diff in [0.1, 0.3]")
  row <- .equiv_one(h, perc$diff, 0.05)
  expect_equal(row$p.lower, (sum(v[, "diff"] <= 0.1) + 1) / 11)
  expect_equal(row$p.upper, (sum(v[, "diff"] >= 0.3) + 1) / 11)
  expect_identical(row$ci_method, "bootstrap (perc)")
  expect_equal(c(row$conf.low, row$conf.high),
               as.vector(.boot_ci(v[, "diff"], 0.2, 0.9, "perc")))
  expect_equal(row$se, stats::sd(v[, "diff"]))

  # The basic distribution of the correlation is reflected on Fisher-z.
  basic <- .equiv_boot(fake_boot("basic"))
  row <- .equiv_one(.parse_hypothesis("rho > 0.1"), basic$rho, 0.05)
  d <- 2 * atanh(0.4) - atanh(v[, "rho"])
  expect_equal(row$p.value, (sum(d <= atanh(0.1)) + 1) / 11)
  expect_equal(c(row$conf.low, row$conf.high),
               as.vector(.boot_ci(v[, "rho"], 0.4, 0.9, "basic",
                                  link = .fisher_link)))

  # The normal distribution has the bias and spread of the replicates.
  norm <- .equiv_boot(fake_boot("norm"))
  row <- .equiv_one(.parse_hypothesis("slope < 0.5"), norm$slope, 0.05)
  centre <- 2 * 0.3 - mean(v[, "slope"])
  expect_equal(row$p.value, stats::pnorm(0.5, centre, stats::sd(v[, "slope"]),
                                         lower.tail = FALSE))

  # A replicate at a correlation of 1 has no Fisher-z value, but the
  # percentile bootstrap still tests it.
  boot <- fake_boot("basic")
  boot$t[1L, "rho"] <- 1
  expect_false(.equiv_boot(boot)$rho$ok)
  boot$type <- "perc"
  expect_true(.equiv_boot(boot)$rho$ok)

  # usdt_equivalence() uses the replicates once usdt_boot() has run.
  m <- structure(list(estimation = "frequentist", boot = fake_boot()),
                 class = "hsdt")
  e <- usdt_equivalence(m, "intercept in [-0.5, 0.5]")
  expect_identical(e$ci_method, "bootstrap (perc)")
})

test_that("each hypothesis prints as a report", {
  m <- structure(list(estimation = "frequentist", boot = fake_boot()),
                 class = "hsdt")
  e <- usdt_equivalence(m, c("diff in [0.1, 0.3]", "rho > 0", "diff = 0"))

  out <- utils::capture.output(print(e))
  expect_length(grep("Equivalence test (TOST)", out, fixed = TRUE), 1L)
  expect_length(grep("p-value, lower bound", out, fixed = TRUE), 1L)
  expect_length(grep("90% CI, bootstrap (perc)", out, fixed = TRUE), 2L)
  expect_length(grep("95% CI, bootstrap (perc)", out, fixed = TRUE), 1L)
  expect_length(grep("Decision (alpha = 0.05)", out, fixed = TRUE), 3L)
  expect_output(print(e[, c("hypothesis", "p.value")]), "rho > 0",
                fixed = TRUE)
})
