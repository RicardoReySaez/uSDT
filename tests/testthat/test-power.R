# test-power.R
# This script tests the approximate power of usdt_power(): the error variance
# of d', the closed-form standard errors against those of the analysis, the
# decisions against usdt_equivalence(), the moment engine, the design grid,
# the reports and the argument checks.
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# This function gives a population with moderate reliability.
power_population <- function(...) {
  usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
                  sd_dprime = c(direct = 0.35, indirect = 0.30),
                  cor_dprime = 0.3, criterion = 0.1, sd_criterion = 0.3, ...)
}

# This function gives the estimates of Vadillo et al. (2025), whose direct
# task has a reliability near 0.13.
vadillo_population <- function() {
  usdt_population(dprime    = c(direct = 0.235, indirect = 0.128),
                  sd_dprime = c(direct = 0.122, indirect = 0.088),
                  cor_dprime = 0.49, criterion = 0.036, sd_criterion = 0.25)
}

test_that("the error variance is the d' variance of usdt_reliability()", {
  # Subjects that barely differ have the variance of the mean subject.
  pop <- as.list(as.data.frame(usdt_population(
    dprime = c(direct = 0.8, indirect = 0.4), sd_dprime = 1e-4,
    cor_dprime = 0, criterion = 0.2, sd_criterion = 0))[1L, ])
  cells <- list(signal = c(direct = 30, indirect = 50),
                noise = c(direct = 34, indirect = 50))
  v <- .error_variance(pop, cells, c(direct = FALSE, indirect = TRUE))
  direct <- .sensitivity_variance(data.frame(sig = c(FALSE, TRUE), n = c(34, 30)),
                                  cbind(-0.2 - 0.4, -0.2 + 0.4), free = TRUE)
  indirect <- .sensitivity_variance(data.frame(sig = c(FALSE, TRUE), n = c(50, 50)),
                                    cbind(-0.2, 0.2), free = FALSE)
  expect_equal(unname(v), c(direct, indirect), tolerance = 1e-6)

  # More trials, less error.
  more <- .error_variance(pop, list(signal = 2 * cells$signal,
                                    noise = 2 * cells$noise),
                          c(direct = FALSE, indirect = TRUE))
  expect_true(all(more < v))
})

test_that("the expected reliability matches a fitted model", {
  set.seed(5)
  trials <- make_trials(n_subj = 200L, n_trials = 100L)
  d <- usdt_data_long(trials, task_col = "task",
                      task_levels      = c(direct = "D", indirect = "I"),
                      subject_col      = "subj",
                      condition_col    = "cond",
                      condition_levels = c(signal = 1, noise = 0),
                      response_col     = "response",
                      response_levels  = c(signal = 1, noise = 0))
  fitted <- usdt_reliability(hsdt(d))$tasks$reliability
  pop <- usdt_population(dprime = c(direct = 0.8, indirect = 0.3),
                         sd_dprime = c(direct = 0.5, indirect = 0.3),
                         cor_dprime = 0.5, criterion = 0, sd_criterion = 0.3)
  p <- usdt_power(pop, subjects = 200, trials = 100)
  expect_equal(c(p$reliability_direct[1L], p$reliability_indirect[1L]), fitted,
               tolerance = 0.1)
})

test_that("the closed-form SEs are those the analysis computes", {
  # The covariance of the five values from the Wishart formula, written out.
  s2D <- 0.35^2; s2I <- 0.30^2; sDI <- 0.3 * 0.35 * 0.30
  v <- c(0.10, 0.03); N <- 80
  S <- matrix(c(s2D + v[1L], sDI, sDI, s2I + v[2L]), 2L)
  idx <- list(c(1, 1), c(2, 2), c(1, 2))
  V <- matrix(0, 5L, 5L, dimnames = rep(list(c("gamma_D", "gamma_I", "s2_D",
                                               "s2_I", "s_DI")), 2L))
  V[1:2, 1:2] <- S / N
  for (a in 1:3) for (b in 1:3) {
    i <- idx[[a]][1L]; j <- idx[[a]][2L]; k <- idx[[b]][1L]; l <- idx[[b]][2L]
    V[2 + a, 2 + b] <- (S[i, k] * S[j, l] + S[i, l] * S[j, k]) / N
  }
  est <- c(gamma_D = 0.15, gamma_I = 0.30, s2_D = s2D, s2_I = s2I, s_DI = sDI)
  p <- list(est = est, vcov = V, rho = 0.3, joint_ok = TRUE,
            fixed_vcov = matrix(V[1:2, 1:2], 2L,
                                dimnames = rep(list(c("d_D", "d_I")), 2L)),
            fixed_names = c("d_D", "d_I"))
  rows <- rbind(.diff_rows(p, 0.95), .cor_rows(p, 0.95), .reg_rows(p, 0.95))

  mine <- .moment_tests(t(est), S[1, 1], S[2, 2], S[1, 2], N)
  expect_equal(mine$diff$se, rows$se[1L])
  expect_equal(mine$rho$se, rows$se[2L])
  expect_equal(mine$intercept$se, rows$se[3L])
  expect_equal(mine$slope$se, rows$se[4L])
  expect_equal(mine$covariance$se, sqrt(V[["s_DI", "s_DI"]]))
})

test_that("the Vadillo design gives the standard errors of the real fit", {
  # hsdt() on vadillo_awareness and vadillo_cuing: 0.0354, 0.519 and 0.116.
  p <- usdt_power(vadillo_population(), subjects = 104,
                  trials = c(direct = 64, indirect = 384),
                  dichotomize = "indirect")
  expect_equal(p$se, c(0.0354, 0.519, 0.116), tolerance = 0.02)
  expect_equal(p$reliability_direct[1L], 0.129, tolerance = 0.02)
  expect_equal(p$reliability_indirect[1L], 0.323, tolerance = 0.02)
})

test_that("each decision agrees with usdt_equivalence()", {
  set.seed(3)
  statements <- c("q = 0.1", "q > 0.1", "q < 0.1", "q in [-0.1, 0.2]",
                  "q out [-0.1, 0.2]")
  for (text in statements) {
    h <- .parse_hypothesis(sub("^q", "diff", text))
    for (i in 1:40) {
      est <- stats::rnorm(1L, 0.05, 0.15)
      se <- stats::runif(1L, 0.02, 0.1)
      s <- list(ok = TRUE, estimate = est, se = se, method = "Wald",
                cdf = function(x) stats::pnorm(x, est, se),
                sf = function(x) stats::pnorm(x, est, se, lower.tail = FALSE),
                ci = function(level) est + c(-1, 1) *
                  stats::qnorm(1 - (1 - level) / 2) * se)
      row <- .equiv_one(h, s, 0.05)
      expect_identical(.power_reject(h$op, h$value, est, se, 0.05),
                       row$decision == "reject H0", info = text)
    }
  }
})

test_that("the normal power formulas match the decisions they summarise", {
  z <- .with_seed(9L, stats::rnorm(2e5))
  for (op in c("=", ">", "<", "in", "out")) {
    bounds <- if (op %in% c("in", "out")) c(-0.1, 0.2) else 0.1
    for (theta in c(-0.05, 0.1, 0.25)) {
      share <- mean(.power_reject(op, bounds, theta + 0.07 * z, 0.07, 0.05))
      expect_lt(abs(.power_normal(op, bounds, theta, 0.07, 0.05) - share),
                0.005, label = paste(op, theta))
    }
  }
})

test_that("the moment engine follows the closed form when nothing is on the boundary", {
  p <- usdt_power(power_population(), subjects = c(60, 120),
                  trials = c(direct = 128, indirect = 400),
                  dichotomize = "indirect")
  expect_true(all(p$estimable > 0.999))
  closed <- vapply(seq_len(nrow(p)), function(i) {
    .power_normal("=", 0, p$true[i], p$se[i], 0.05)
  }, 0)
  h1 <- p$hypothesis == paste0(.bf_name("diff"), " = 0")
  expect_equal(p$power[h1], closed[h1], tolerance = 0.02)
})

test_that("power rises with subjects and trials, and boundary fits fall", {
  p <- usdt_power(vadillo_population(), subjects = c(50, 100, 200, 400),
                  trials = list(direct = c(32, 64, 128), indirect = 384),
                  dichotomize = "indirect")
  for (h in unique(p$hypothesis)) {
    for (td in c(32, 64, 128)) {
      by_n <- p$power[p$hypothesis == h & p$trials_direct == td]
      expect_true(all(diff(by_n) >= -0.005), info = paste(h, td))
    }
    by_t <- p$estimable[p$hypothesis == h & p$subjects == 100]
    expect_true(all(diff(by_t) >= -0.005), info = h)
  }

  # With a direct reliability near 0.13, many fits are on the boundary.
  low <- p[p$subjects == 100 & p$trials_direct == 64, ]
  expect_true(all(low$estimable[low$quantity != "diff"] < 0.9))
})

test_that("a population in H0 is flagged and rejects near alpha", {
  p <- usdt_power(update(power_population(), cor_dprime = 0), subjects = 200,
                  trials = c(direct = 128, indirect = 400),
                  hypothesis = c("rho = 0", "rho > 0.2", "diff > 0"),
                  dichotomize = "indirect")
  expect_identical(p$in_null, c(TRUE, TRUE, FALSE))
  expect_equal(p$power[1L], 0.05, tolerance = 0.3)
  expect_lt(p$power[2L], 0.05)
})

test_that("the approximation is reproducible and leaves the session's random numbers alone", {
  set.seed(1)
  before <- stats::runif(1L)
  set.seed(1)
  a <- usdt_power(power_population(), subjects = 80, trials = 100)
  after <- stats::runif(1L)
  expect_identical(before, after)
  expect_identical(a$power,
                   usdt_power(power_population(), subjects = 80,
                              trials = 100)$power)
})

test_that("scenarios, subjects and trials make a grid", {
  pop <- update(power_population(), cor_dprime = c(0.1, 0.3))
  p <- usdt_power(pop, subjects = c(40, 80),
                  trials = list(direct = c(32, 64), indirect = 200),
                  combine = "all")
  expect_s3_class(p, "usdt_power")
  expect_equal(nrow(p), 2L * 2L * 2L * 4L)
  expect_setequal(unique(p$hypothesis),
                  c(paste0(.bf_name("diff"), " = 0"), "rho = 0",
                    "intercept = 0", "all"))
  joint <- p[p$hypothesis == "all", ]
  each <- p[p$hypothesis != "all", ]
  lowest <- tapply(each$power, paste(each$scenario, each$subjects,
                                     each$trials_direct), min)
  expect_true(all(joint$power <=
                    lowest[paste(joint$scenario, joint$subjects,
                                 joint$trials_direct)] + 1e-12))

  # A shared trial count keeps both tasks equal.
  shared <- usdt_power(power_population(), subjects = 40, trials = c(50, 100))
  expect_identical(shared$trials_direct, shared$trials_indirect)
})

test_that("the report is a card for one design and a table for several", {
  one <- capture.output(print(usdt_power(power_population(), subjects = 80,
                                         trials = 100, combine = "all")))
  expect_true(any(grepl("80 subjects", one)))
  expect_true(any(grepl("All at once", one)))

  many <- capture.output(print(usdt_power(update(power_population(),
                                                 cor_dprime = 0),
                                          subjects = c(40, 80), trials = 100)))
  expect_true(any(grepl("2 designs", many[1L])))
  expect_true(any(grepl("type I error", many)))
  expect_false(inherits(usdt_power(power_population(), subjects = 40,
                                   trials = 100)[1L, ], "usdt_power"))
})

test_that("invalid arguments stop with a message that names them", {
  pop <- power_population()
  expect_error(usdt_power(data.frame(), 40, 100), "usdt_population")
  expect_error(usdt_power(update(pop, sd_signal = 1.2), 40, 100),
               "equal-variance")
  expect_error(usdt_power(pop, subjects = 2, trials = 100), "at least 3")
  expect_error(usdt_power(pop, subjects = 40.5, trials = 100), "whole")
  expect_error(usdt_power(pop, subjects = 40, trials = 0), "positive whole")
  expect_error(usdt_power(pop, subjects = 40, trials = 101), "whole numbers")
  expect_error(usdt_power(pop, subjects = 40, trials = 100, signal = 1),
               "inside \\(0, 1\\)")
  expect_error(usdt_power(pop, subjects = 40, trials = 100,
                          hypothesis = "beta = 0"), "quantities")
  expect_error(usdt_power(pop, subjects = 40, trials = 100, alpha = 0.6),
               "alpha")
})
