# test-population.R
# This script tests the population values of usdt_population(): both
# parameterizations, the tested quantities, the scenarios, the values taken
# from a fitted model and update().
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# This function gives the population used by most tests.
base_population <- function() {
  usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
                  sd_dprime = c(direct = 0.35, indirect = 0.30),
                  cor_dprime = 0.3)
}

test_that("the correlation parameterization derives the tested quantities", {
  pop <- base_population()
  expect_s3_class(pop, "usdt_population")
  expect_identical(attr(pop, "parameterization"), "correlation")
  expect_equal(nrow(pop), 1L)

  q <- .usdt_quantities(c(gamma_D = 0.15, gamma_I = 0.30, s2_D = 0.35^2,
                          s2_I = 0.30^2, s_DI = 0.3 * 0.35 * 0.30))
  expect_equal(pop$diff, q[, "diff"][[1L]])
  expect_equal(pop$intercept, q[, "intercept"][[1L]])
  expect_equal(pop$slope, q[, "slope"][[1L]])
  expect_equal(pop$rho, 0.3)

  # The criteria take their defaults for both tasks.
  expect_equal(c(pop$criterion_D, pop$criterion_I), c(0, 0))
  expect_equal(c(pop$sd_criterion_D, pop$sd_criterion_I), c(0.3, 0.3))
  expect_equal(c(pop$sd_signal_D, pop$sd_signal_I), c(1, 1))
})

test_that("the latent regression gives the intercept and derives the indirect mean", {
  pop <- usdt_population(dprime = 0.15, intercept = 0.25,
                         sd_dprime = c(direct = 0.35, indirect = 0.30),
                         cor_dprime = 0.3)
  expect_identical(attr(pop, "parameterization"), "regression")
  expect_equal(pop$mu_I, 0.25 + 0.3 * 0.30 / 0.35 * 0.15)
  expect_equal(pop$intercept, 0.25)
  expect_equal(pop$diff, pop$mu_I - 0.15)

  # Both parameterizations describe the same distribution.
  same <- usdt_population(dprime = c(direct = 0.15, indirect = pop$mu_I),
                          sd_dprime = c(direct = 0.35, indirect = 0.30),
                          cor_dprime = 0.3)
  expect_equal(same$intercept, 0.25)
  expect_equal(same$slope, pop$slope)

  # The direct mean may be named.
  named <- usdt_population(dprime = c(direct = 0.15), intercept = 0.25,
                           sd_dprime = 0.3, cor_dprime = 0.3)
  expect_equal(named$mu_D, 0.15)
})

test_that("shared values are linked across tasks and task values are crossed", {
  shared <- usdt_population(dprime = c(direct = 0.1, indirect = 0.2),
                            sd_dprime = c(0.2, 0.3), cor_dprime = 0.3)
  expect_equal(nrow(shared), 2L)
  expect_equal(shared$sigma_D, shared$sigma_I)

  crossed <- usdt_population(
    dprime = c(direct = 0.1, indirect = 0.2),
    sd_dprime = list(direct = c(0.2, 0.3), indirect = c(0.2, 0.3)),
    cor_dprime = 0.3)
  expect_equal(nrow(crossed), 4L)

  # Every argument with several values multiplies the scenarios.
  grid <- usdt_population(dprime = c(direct = 0.1, indirect = 0.2),
                          sd_dprime = c(0.2, 0.3),
                          cor_dprime = c(0.1, 0.3, 0.5))
  expect_equal(nrow(grid), 6L)
  expect_equal(grid$scenario, 1:6)
  expect_setequal(paste(grid$sigma_D, grid$rho),
                  paste(rep(c(0.2, 0.3), 3L), rep(c(0.1, 0.3, 0.5), each = 2L)))
})

test_that("invalid values stop with a message that names the argument", {
  expect_error(usdt_population(dprime = 0.1, cor_dprime = 0.3), "`sd_dprime`")
  expect_error(usdt_population(dprime = c(direct = 0.1), sd_dprime = 0.3,
                               cor_dprime = 0.3), "indirect mean")
  expect_error(usdt_population(dprime = 0.1, sd_dprime = 0, cor_dprime = 0.3),
               "must be positive")
  expect_error(usdt_population(dprime = 0.1, sd_dprime = 0.3, cor_dprime = 1),
               "inside \\(-1, 1\\)")
  expect_error(usdt_population(dprime = 0.1, sd_dprime = 0.3, cor_dprime = 0.3,
                               sd_criterion = -0.1), "cannot be negative")
  expect_error(usdt_population(dprime = c(direct = 0.1, indirect = 0.2),
                               intercept = 0.1, sd_dprime = 0.3,
                               cor_dprime = 0.3), "direct mean only")
  expect_error(usdt_population(dprime = 0.1, sd_dprime = 0.3,
                               cor_dprime = c(direct = 0.3, indirect = 0.3)),
               "single value")
  expect_error(usdt_population(dprime = NA, sd_dprime = 0.3, cor_dprime = 0.3),
               "finite numbers")
})

test_that("update() changes values, one task at a time if needed", {
  pop <- base_population()
  expect_equal(nrow(update(pop, cor_dprime = c(0.1, 0.3, 0.5))), 3L)

  one_task <- update(pop, sd_dprime = c(direct = 0.2))
  expect_equal(c(one_task$sigma_D, one_task$sigma_I), c(0.2, 0.30))

  # Under the correlation parameterization a new correlation keeps the means.
  weaker <- update(pop, cor_dprime = 0.1)
  expect_equal(weaker$mu_I, pop$mu_I)
  expect_false(isTRUE(all.equal(weaker$intercept, pop$intercept)))

  # An intercept switches to the latent regression, which keeps it fixed.
  reg <- update(pop, intercept = 0)
  expect_identical(attr(reg, "parameterization"), "regression")
  expect_equal(reg$mu_D, 0.15)
  expect_equal(update(reg, cor_dprime = 0.6)$intercept, 0)

  # Removing the intercept needs the indirect mean back.
  expect_error(update(reg, intercept = NULL), "indirect mean")
  back <- update(reg, intercept = NULL,
                 dprime = c(direct = 0.15, indirect = 0.30))
  expect_equal(back$intercept, pop$intercept)

  expect_error(update(pop, slope = 1), "not an argument")
})

test_that("a subset of the scenarios is a plain data frame", {
  pop <- update(base_population(), cor_dprime = c(0.1, 0.3))
  expect_false(inherits(pop[1L, ], "usdt_population"))
  expect_s3_class(pop[1L, ], "data.frame")
})

test_that("print() shows the parameters, the quantities and the scenarios", {
  out <- capture.output(print(base_population()))
  expect_true(any(grepl("correlation", out[1L])))
  expect_true(any(grepl("intercept +0\\.261", out)))

  reg <- capture.output(print(usdt_population(dprime = 0.15, intercept = 0.25,
                                              sd_dprime = 0.3,
                                              cor_dprime = 0.3)))
  expect_true(any(grepl("latent regression", reg[1L])))
  expect_true(any(grepl("follows from the intercept", reg)))

  scenarios <- capture.output(print(update(base_population(),
                                           cor_dprime = c(0.1, 0.5))))
  expect_true(any(grepl("2 scenarios", scenarios[1L])))
  expect_true(any(grepl("varies", scenarios)))
})

test_that("from = takes the estimates of a frequentist fit", {
  set.seed(31)
  trials <- make_trials(n_subj = 40L, n_trials = 80L)
  d <- usdt_data_long(trials, task_col = "task",
                      task_levels      = c(direct = "D", indirect = "I"),
                      subject_col      = "subj",
                      condition_col    = "cond",
                      condition_levels = c(signal = 1, noise = 0),
                      response_col     = "response",
                      response_levels  = c(signal = 1, noise = 0))
  m <- hsdt(d)
  e <- m$pars$est
  pop <- usdt_population(from = m)
  expect_equal(c(pop$mu_D, pop$mu_I), unname(e[c("gamma_D", "gamma_I")]))
  expect_equal(c(pop$sigma_D, pop$sigma_I), sqrt(unname(e[c("s2_D", "s2_I")])))
  expect_equal(pop$rho, m$pars$rho)

  # The criteria return to the classical scale, with the fitted SDs.
  beta <- lme4::fixef(m$fit)
  expect_equal(pop$criterion_D, -beta[["c_D"]])
  expect_equal(pop$sd_criterion_D, m$pars$re$sd[m$pars$re$term == "c_D"])

  # The call replaces single values, and an intercept the indirect mean.
  expect_equal(usdt_population(from = m, cor_dprime = 0.2)$rho, 0.2)
  reg <- usdt_population(from = m, intercept = 0)
  expect_identical(attr(reg, "parameterization"), "regression")
  expect_equal(reg$mu_D, pop$mu_D)
  expect_equal(reg$intercept, 0)
})

test_that("from = rejects Bayesian fits, other objects and boundary estimates", {
  expect_error(usdt_population(from = list()), "must come from hsdt")
  bayes <- structure(list(estimation = "bayesian"), class = "hsdt")
  expect_error(usdt_population(from = bayes), "frequentist")

  boundary <- structure(list(
    estimation = "frequentist",
    pars = list(est = c(gamma_D = 0.2, gamma_I = 0.1, s2_D = 0, s2_I = 0.04,
                        s_DI = 0),
                rho = 0, re = data.frame(term = character(0), sd = numeric(0)),
                criterion_cor = NULL),
    design = list(criteria = character(0)),
    data = list(meta = list(coding = "deviation",
                            labels = c(direct = "D", indirect = "I"))),
    fit = NULL), class = "hsdt")
  local_mocked_bindings(fixef = function(object, ...) c(d_D = 0.2, d_I = 0.1),
                        .package = "lme4")
  expect_error(usdt_population(from = boundary), "`sd_dprime`.*boundary")
  expect_equal(usdt_population(from = boundary, sd_dprime = 0.2)$sigma_D, 0.2)
})
