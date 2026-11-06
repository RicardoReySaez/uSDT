# test-power.R
# This script tests the approximate power of usdt_power(): the error variance
# of d', the closed-form standard errors against those of the analysis, the
# decisions against usdt_equivalence(), the moment engine, the design grid,
# the reports and the argument checks.
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# This function gives a population with moderate reliability.
power_population <- function(...) {
  usdt_population(dprime    = c(awareness = 0.15, influence = 0.30),
                  sd_dprime = c(awareness = 0.35, influence = 0.30),
                  cor_dprime = 0.3, criterion = 0.1, sd_criterion = 0.3, ...)
}

# This function gives the estimates of Vadillo et al. (2025), whose awareness
# task has a reliability near 0.13.
vadillo_population <- function() {
  usdt_population(dprime    = c(awareness = 0.235, influence = 0.128),
                  sd_dprime = c(awareness = 0.122, influence = 0.088),
                  cor_dprime = 0.49, criterion = 0.036, sd_criterion = 0.25)
}

test_that("the error variance is the d' variance of usdt_reliability()", {
  # Subjects that barely differ have the variance of the mean subject.
  pop <- as.list(as.data.frame(usdt_population(
    dprime = c(awareness = 0.8, influence = 0.4), sd_dprime = 1e-4,
    cor_dprime = 0, criterion = 0.2, sd_criterion = 0))[1L, ])
  cells <- list(signal = c(awareness = 30, influence = 50),
                noise = c(awareness = 34, influence = 50))
  v <- .error_variance(pop, cells, c(awareness = FALSE, influence = TRUE))
  awareness <- .sensitivity_variance(
    data.frame(sig = c(FALSE, TRUE), n = c(34, 30)),
    cbind(-0.2 - 0.4, -0.2 + 0.4), free = TRUE)
  influence <- .sensitivity_variance(
    data.frame(sig = c(FALSE, TRUE), n = c(50, 50)),
    cbind(-0.2, 0.2), free = FALSE)
  expect_equal(unname(v), c(awareness, influence), tolerance = 1e-6)

  # More trials, less error.
  more <- .error_variance(pop, list(signal = 2 * cells$signal,
                                    noise = 2 * cells$noise),
                          c(awareness = FALSE, influence = TRUE))
  expect_true(all(more < v))
})

test_that("the expected reliability matches a fitted model", {
  set.seed(5)
  trials <- make_trials(n_subj = 200L, n_trials = 100L)
  d <- usdt_data_long(trials, task_col = "task",
                      task_levels      = c(awareness = "A", influence = "I"),
                      subject_col      = "subj",
                      condition_col    = "cond",
                      condition_levels = c(signal = 1, noise = 0),
                      response_col     = "response",
                      response_levels  = c(signal = 1, noise = 0))
  fitted <- usdt_reliability(hsdt(d))$tasks$reliability
  pop <- usdt_population(dprime = c(awareness = 0.8, influence = 0.3),
                         sd_dprime = c(awareness = 0.5, influence = 0.3),
                         cor_dprime = 0.5, criterion = 0, sd_criterion = 0.3)
  p <- usdt_power(pop, subjects = 200, trials = 100)
  expect_equal(c(p$reliability_awareness[1L], p$reliability_influence[1L]),
               fitted,
               tolerance = 0.1)
})

test_that("the closed-form SEs are those the analysis computes", {
  # The covariance of the five values from the Wishart formula, written out.
  s2A <- 0.35^2; s2I <- 0.30^2; sAI <- 0.3 * 0.35 * 0.30
  v <- c(0.10, 0.03); N <- 80
  S <- matrix(c(s2A + v[1L], sAI, sAI, s2I + v[2L]), 2L)
  idx <- list(c(1, 1), c(2, 2), c(1, 2))
  V <- matrix(0, 5L, 5L, dimnames = rep(list(c("gamma_A", "gamma_I", "s2_A",
                                               "s2_I", "s_AI")), 2L))
  V[1:2, 1:2] <- S / N
  for (a in 1:3) for (b in 1:3) {
    i <- idx[[a]][1L]; j <- idx[[a]][2L]; k <- idx[[b]][1L]; l <- idx[[b]][2L]
    V[2 + a, 2 + b] <- (S[i, k] * S[j, l] + S[i, l] * S[j, k]) / N
  }
  est <- c(gamma_A = 0.15, gamma_I = 0.30, s2_A = s2A, s2_I = s2I, s_AI = sAI)
  p <- list(est = est, vcov = V, rho = 0.3, joint_ok = TRUE,
            fixed_vcov = matrix(V[1:2, 1:2], 2L,
                                dimnames = rep(list(c("d_A", "d_I")), 2L)),
            fixed_names = c("d_A", "d_I"))
  rows <- rbind(.diff_rows(p, 0.95), .cor_rows(p, 0.95), .reg_rows(p, 0.95))

  mine <- .moment_tests(t(est), S[1, 1], S[2, 2], S[1, 2], N)
  expect_equal(mine$diff$se, rows$se[1L])
  expect_equal(mine$rho$se, rows$se[2L])
  expect_equal(mine$intercept$se, rows$se[3L])
  expect_equal(mine$slope$se, rows$se[4L])
  expect_equal(mine$covariance$se, sqrt(V[["s_AI", "s_AI"]]))
})

test_that("the Vadillo design gives the standard errors of the real fit", {
  # hsdt() on vadillo_awareness and vadillo_cuing: 0.0354, 0.519 and 0.116.
  p <- usdt_power(vadillo_population(), subjects = 104,
                  trials = c(awareness = 64, influence = 384),
                  dichotomize = "influence")
  expect_equal(p$se, c(0.0354, 0.519, 0.116), tolerance = 0.02)
  expect_equal(p$reliability_awareness[1L], 0.129, tolerance = 0.02)
  expect_equal(p$reliability_influence[1L], 0.323, tolerance = 0.02)
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
                  trials = c(awareness = 128, influence = 400),
                  dichotomize = "influence")
  expect_true(all(p$estimable > 0.999))
  closed <- vapply(seq_len(nrow(p)), function(i) {
    .power_normal("=", 0, p$true[i], p$se[i], 0.05)
  }, 0)
  h1 <- p$hypothesis == paste0(.bf_name("diff"), " = 0")
  expect_equal(p$power[h1], closed[h1], tolerance = 0.02)
})

test_that("power rises with subjects and trials, and boundary fits fall", {
  p <- usdt_power(vadillo_population(), subjects = c(50, 100, 200, 400),
                  trials = list(awareness = c(32, 64, 128), influence = 384),
                  dichotomize = "influence")
  for (h in unique(p$hypothesis)) {
    for (td in c(32, 64, 128)) {
      by_n <- p$power[p$hypothesis == h & p$trials_awareness == td]
      expect_true(all(diff(by_n) >= -0.005), info = paste(h, td))
    }
    by_t <- p$estimable[p$hypothesis == h & p$subjects == 100]
    expect_true(all(diff(by_t) >= -0.005), info = h)
  }

  # With an awareness reliability near 0.13, many fits are on the boundary.
  low <- p[p$subjects == 100 & p$trials_awareness == 64, ]
  expect_true(all(low$estimable[low$quantity != "diff"] < 0.9))
})

test_that("a population in H0 is flagged and rejects near alpha", {
  p <- usdt_power(update(power_population(), cor_dprime = 0), subjects = 200,
                  trials = c(awareness = 128, influence = 400),
                  hypothesis = c("rho = 0", "rho > 0.2", "diff > 0"),
                  dichotomize = "influence")
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

test_that("scoped seeds restore RNG settings even when a calculation fails", {
  set.seed(41)
  before <- .Random.seed
  kind <- RNGkind()
  expect_error(.with_seed(7, stop("calculation failed"), kind = "L'Ecuyer-CMRG"),
               "calculation failed")
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
})

test_that("scoped seeds leave a fresh session without a seed or a new generator", {
  scope <- .with_seed
  environment(scope) <- baseenv()
  result <- callr::r(function(scope) {
    stopifnot(!exists(".Random.seed", globalenv(), inherits = FALSE))
    before <- ls(globalenv(), all.names = TRUE)
    kind <- RNGkind()
    scope(7, stats::runif(2), kind = "L'Ecuyer-CMRG")
    c(objects = identical(ls(globalenv(), all.names = TRUE), before),
      generator = identical(RNGkind(), kind))
  }, args = list(scope = scope))
  expect_true(result[["objects"]])
  expect_true(result[["generator"]])
})

test_that("scenarios, subjects and trials make a grid", {
  pop <- update(power_population(), cor_dprime = c(0.1, 0.3))
  p <- usdt_power(pop, subjects = c(40, 80),
                  trials = list(awareness = c(32, 64), influence = 200),
                  combine = "all")
  expect_s3_class(p, "usdt_power")
  expect_equal(nrow(p), 2L * 2L * 2L * 4L)
  expect_setequal(unique(p$hypothesis),
                  c(paste0(.bf_name("diff"), " = 0"), "rho = 0",
                    "intercept = 0", "all"))
  joint <- p[p$hypothesis == "all", ]
  each <- p[p$hypothesis != "all", ]
  lowest <- tapply(each$power, paste(each$scenario, each$subjects,
                                     each$trials_awareness), min)
  expect_true(all(joint$power <=
                    lowest[paste(joint$scenario, joint$subjects,
                                 joint$trials_awareness)] + 1e-12))

  # A shared trial count keeps both tasks equal.
  shared <- usdt_power(power_population(), subjects = 40, trials = c(50, 100))
  expect_identical(shared$trials_awareness, shared$trials_influence)
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
