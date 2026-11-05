# test-plan.R
# This script tests usdt_plan() and usdt_mdes(): the design space and costs,
# the searches, the chosen designs against usdt_power(), the smallest
# effects, and the use of a plan by usdt_power() and usdt_mdes().
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# This function gives a population with moderate reliability.
plan_population <- function(...) {
  usdt_population(dprime    = c(awareness = 0.15, influence = 0.30),
                  sd_dprime = c(awareness = 0.35, influence = 0.30),
                  cor_dprime = 0.3, criterion = 0.1, ...)
}

# This function gives the approximate power of one hypothesis at one design.
power_at <- function(pop, N, TA, TI, hypothesis, ...) {
  usdt_power(pop, subjects = N, trials = c(awareness = TA, influence = TI),
             hypothesis = hypothesis, dichotomize = "influence", ...)$power
}

test_that("the design space reads fixed values, ranges and whole trial counts", {
  s <- .plan_space(c(20, 300), list(awareness = 64, influence = c(100, 801)),
                   0.5)
  expect_identical(s$free, c(subjects = TRUE, awareness = FALSE,
                             influence = TRUE))
  expect_identical(unique(s$grid$awareness), 64L)
  expect_true(all(s$grid$influence %% 2L == 0L))
  expect_identical(range(s$grid$influence), c(100L, 800L))

  linked <- .plan_space(50, c(30, 90), 0.5)
  expect_true(linked$linked)
  expect_identical(linked$grid$awareness, linked$grid$influence)

  # Large ranges are searched on at most 400 values, ends included.
  wide <- .plan_space(50, list(awareness = c(10, 5000), influence = 100), 0.5)
  expect_lte(length(unique(wide$grid$awareness)), 400L)
  expect_identical(range(wide$grid$awareness), c(10L, 5000L))

  # A share of signal trials of a quarter needs multiples of four.
  quarter <- .plan_space(50, list(awareness = c(10, 30), influence = 40), 0.25)
  expect_true(all(quarter$grid$awareness %% 4L == 0L))
  expect_error(.plan_space(50, list(awareness = 30, influence = 40), 0.25),
               "whole numbers")
  expect_error(.plan_space(2, 100, 0.5), "at least 3")
  expect_error(.plan_space(50, c(10, 20, 30), 0.5), "range")
})

test_that("costs are linear, a function, or the single free value", {
  linear <- .plan_cost(c(subject = 10, influence = 0.5),
                       c(subjects = TRUE, awareness = TRUE, influence = FALSE))
  expect_equal(linear(4, 20, 100), 4 * (10 + 0.5 * 100))
  custom <- .plan_cost(function(subjects, awareness, influence)
                         subjects + awareness,
                       c(subjects = TRUE, awareness = TRUE, influence = FALSE))
  expect_equal(custom(c(1, 2), c(10, 20), c(0, 0)), c(11, 22))
  only <- .plan_cost(NULL, c(subjects = FALSE, awareness = TRUE,
                             influence = FALSE))
  expect_equal(only(50, 64, 200), 64)
  expect_error(.plan_cost(c(subjects = 1),
                          c(subjects = TRUE, awareness = FALSE,
                            influence = FALSE)), "named vector")
  expect_error(.plan_cost(function(subjects, awareness, influence) 1,
                          c(subjects = TRUE))(1:2, 1:2, 1:2), "one number")
})

test_that("the searches find the smallest N that reaches the target", {
  # Three designs whose power grows with sqrt(N) at different rates; the last
  # does not reach 0.8 below 2000 subjects.
  effects <- c(0.1, 0.2, 0.05)
  f <- function(N) stats::pnorm(sqrt(N) * effects - 1.96)
  brute <- vapply(effects, function(e) {
    N <- 3:2000
    N[which(stats::pnorm(sqrt(N) * e - 1.96) >= 0.8)[1L]]
  }, 0)
  expect_equal(.min_subjects(f, 3, 2000, 3L, 0.8), brute)
  expect_true(is.na(brute[3L]))

  g <- function(N) as.numeric(N >= 137)
  for (guess in c(10, 136, 137, 138, 400)) {
    expect_equal(.refine_subjects(g, guess, 3, 500, 0.5), 137, info = guess)
  }
  expect_true(is.na(.refine_subjects(g, 50, 3, 100, 0.5)))
  expect_equal(.refine_subjects(function(N) 1, 50, 3, 100, 0.5), 3)
})

test_that("the subjects for given trials reach the power, and one fewer does not", {
  pop <- plan_population()
  p <- usdt_plan(pop, hypothesis = "diff = 0", subjects = c(10, 1000),
                 trials = c(awareness = 64, influence = 200),
                 dichotomize = "influence")
  expect_s3_class(p, "usdt_plan")
  N <- p$design$subjects
  expect_gte(power_at(pop, N, 64, 200, "diff = 0"), 0.8)
  expect_lt(power_at(pop, N - 1, 64, 200, "diff = 0"), 0.8)
  expect_equal(p$design$power, power_at(pop, N, 64, 200, "diff = 0"))
  expect_true(is.na(p$design$limit))
})

test_that("with costs, the trials of H1 follow the optimal-design formula", {
  pop <- plan_population()
  p <- usdt_plan(pop, hypothesis = "diff = 0", power = 0.9,
                 subjects = c(10, 1000),
                 trials = list(awareness = c(8, 400), influence = c(20, 1000)),
                 cost = c(subject = 600, awareness = 4, influence = 1.5),
                 dichotomize = "influence")

  # T* = sqrt((c0 / ct) (a / tau^2)), with a from the error variance.
  tau2 <- 0.35^2 + 0.30^2 - 2 * 0.3 * 0.35 * 0.30
  scenario <- as.list(as.data.frame(pop)[1L, ])
  v <- .error_variance(scenario,
                       list(signal = c(awareness = 50, influence = 50),
                            noise = c(awareness = 50, influence = 50)),
                       c(awareness = FALSE, influence = TRUE))
  optimal <- sqrt(c(600 / 4, 600 / 1.5) * 100 * v / tau2)
  expect_equal(c(p$design$trials_awareness, p$design$trials_influence),
               unname(optimal), tolerance = 0.2)
  expect_gte(power_at(pop, p$design$subjects, p$design$trials_awareness,
                      p$design$trials_influence, "diff = 0"), 0.9)

  # The frontier gives every number of trials its smallest N.
  expect_true(all(c("trials_awareness", "trials_influence", "subjects",
                    "cost") %in% names(p$frontier)))
  expect_lte(p$design$cost, min(p$frontier$cost, na.rm = TRUE) * 1.2)
})

test_that("a wider range never gives a dearer design than the one it contains", {
  # The closed form ranks designs with few awareness trials too high for the
  # intercept, so a search led by it alone missed cheaper designs.
  skip_on_cran()
  pop <- plan_population()
  run <- function(influence) {
    usdt_plan(pop, hypothesis = "intercept > 0", power = 0.9,
              subjects = c(20, 300),
              trials = list(awareness = c(32, 160), influence = influence),
              cost = c(subject = 600, awareness = 4, influence = 1.5),
              dichotomize = "influence")$design$cost
  }
  expect_lte(run(c(100, 400)), run(c(100, 300)) * 1.01)
})

test_that("a budget gives the most powerful design that costs no more", {
  pop <- plan_population()
  args <- list(pop, hypothesis = "rho = 0", subjects = c(20, 400),
               trials = list(awareness = c(32, 256), influence = c(100, 800)),
               cost = c(subject = 600, awareness = 4, influence = 1.5),
               dichotomize = "influence")
  small <- do.call(usdt_plan, c(args, budget = 20 * 3600))
  large <- do.call(usdt_plan, c(args, budget = 40 * 3600))
  expect_lte(small$design$cost, 20 * 3600)
  expect_lte(large$design$cost, 40 * 3600)
  expect_gte(large$design$power, small$design$power)
  expect_equal(small$design$power,
               power_at(pop, small$design$subjects,
                        small$design$trials_awareness,
                        small$design$trials_influence, "rho = 0"))
  expect_error(usdt_plan(pop, subjects = c(20, 400), trials = 100,
                         budget = 1000), "`cost`")
})

test_that("a precision target gives the narrowest affordable interval", {
  pop <- plan_population()
  p <- usdt_plan(pop, precision = c(diff = 0.2, rho = 0.6),
                 trials = c(awareness = 64, influence = 200),
                 dichotomize = "influence")
  expect_true(all(p$design$width <= c(0.2, 0.6)))
  diff_row <- p$design[startsWith(p$design$hypothesis, "diff"), ]
  se <- usdt_power(pop, subjects = c(diff_row$subjects - 1, diff_row$subjects),
                   trials = c(awareness = 64, influence = 200),
                   hypothesis = "diff = 0", dichotomize = "influence")$se
  expect_gt(2 * stats::qnorm(0.975) * se[1L], 0.2)
  expect_lte(2 * stats::qnorm(0.975) * se[2L], 0.2)
})

test_that("all hypotheses at once need at least the subjects of each", {
  pop <- plan_population()
  each <- usdt_plan(pop, trials = c(awareness = 64, influence = 200),
                    dichotomize = "influence")
  all <- usdt_plan(pop, trials = c(awareness = 64, influence = 200),
                   combine = "all", dichotomize = "influence")
  expect_gte(all$design$subjects, max(each$design$subjects))
  joint <- usdt_power(pop, subjects = all$design$subjects,
                      trials = c(awareness = 64, influence = 200),
                      combine = "all", dichotomize = "influence")
  expect_gte(joint$power[joint$hypothesis == "all"], 0.8)
})

test_that("unreachable targets and designs at a limit are reported", {
  pop <- plan_population()
  none <- usdt_plan(pop, hypothesis = "rho = 0", subjects = c(10, 40),
                    trials = c(awareness = 64, influence = 200),
                    dichotomize = "influence")
  expect_true(is.na(none$design$subjects))
  expect_output(print(none), "No design in the ranges")

  edge <- usdt_plan(pop, hypothesis = "diff = 0", subjects = c(10, 1000),
                    trials = list(awareness = c(8, 16), influence = 200),
                    cost = c(subject = 600, awareness = 4, influence = 1.5),
                    dichotomize = "influence")
  expect_match(edge$design$limit, "trials\\$awareness")
  expect_output(print(edge), "at the limit")
})

test_that("usdt_plan() checks what it is given", {
  pop <- plan_population()
  expect_error(usdt_plan(pop), "`trials` is required")
  expect_error(usdt_plan(pop, subjects = c(10, 100),
                         trials = list(awareness = c(10, 100),
                                       influence = 100)),
               "`cost`")
  expect_error(usdt_plan(pop, subjects = 50, trials = 100), "nothing is free")
  expect_error(usdt_plan(pop, trials = 100, precision = c(beta = 0.1)),
               "precision")
  expect_error(usdt_plan(pop, trials = 100, power = 1), "power")
  expect_error(usdt_plan(update(pop, sd_signal = 1.2), trials = 100),
               "equal-variance")
})

test_that("a moved quantity keeps fixed what the parameterization fixes", {
  pop <- as.list(as.data.frame(plan_population())[1L, ])
  moved <- .scenario_with(pop, "rho", 0.6, "correlation")
  expect_equal(moved$mu_I, pop$mu_I)
  expect_false(isTRUE(all.equal(moved$intercept, pop$intercept)))

  reg <- as.list(as.data.frame(update(plan_population(), intercept = 0.2))[1L, ])
  moved <- .scenario_with(reg, "rho", 0.6, "regression")
  expect_equal(moved$intercept, 0.2)
  expect_equal(moved$rho, 0.6)

  expect_equal(.scenario_with(pop, "diff", 0.4, "correlation")$diff, 0.4)
  expect_equal(.scenario_with(pop, "intercept", 0.1, "correlation")$intercept,
               0.1)
  expect_equal(.scenario_with(pop, "slope", 0.5, "correlation")$slope, 0.5)
})

test_that("the smallest effect is where the power reaches the target", {
  pop <- plan_population()
  m <- usdt_mdes(pop, subjects = 80,
                 trials = c(awareness = 64, influence = 200),
                 hypothesis = c("diff = 0", "intercept > 0"),
                 dichotomize = "influence")
  expect_s3_class(m, "usdt_mdes")
  expect_identical(m$varies, c("diff", "intercept"))
  scenario <- update(pop, dprime = c(influence = 0.15 + m$mdes[1L]))
  expect_equal(power_at(scenario, 80, 64, 200, "diff = 0"), 0.8,
               tolerance = 0.02)
  expect_length(attr(m, "curve"), 2L)

  # The narrowest equivalence margin.
  same <- update(pop, dprime = c(awareness = 0.15, influence = 0.15))
  tost <- usdt_mdes(same, subjects = 80,
                    trials = c(awareness = 64, influence = 200),
                    hypothesis = "diff in [-0.2, 0.2]",
                    dichotomize = "influence")
  expect_identical(tost$varies, "margin")
  margin <- tost$mdes
  expect_equal(power_at(same, 80, 64, 200,
                        sprintf("diff in [%f, %f]", -margin, margin)),
               0.8, tolerance = 0.02)

  # Another value of the population, and a target out of reach.
  sd_d <- usdt_mdes(pop, subjects = 80,
                    trials = c(awareness = 64, influence = 200),
                    hypothesis = "intercept > 0", vary = "sigma_A",
                    range = c(0.05, 0.6), dichotomize = "influence")
  expect_identical(sd_d$varies, "sigma_A")
  expect_true(is.finite(sd_d$mdes))
  far <- usdt_mdes(pop, subjects = 10,
                   trials = c(awareness = 16, influence = 40),
                   hypothesis = "rho = 0", power = 0.99,
                   dichotomize = "influence")
  expect_true(is.na(far$mdes))
  expect_output(print(m), "Smallest")
  expect_error(usdt_mdes(pop, subjects = 80, trials = 100, vary = "sigma_A"),
               "range")
})

test_that("usdt_power() and usdt_mdes() take the designs of a plan", {
  pop <- plan_population()
  p <- usdt_plan(pop, hypothesis = c("diff = 0", "intercept = 0"),
                 trials = c(awareness = 64, influence = 200),
                 dichotomize = "influence")
  v <- usdt_power(p)
  expect_s3_class(v, "usdt_power")
  expect_setequal(unique(v$subjects), p$design$subjects)
  expect_error(usdt_power(p, subjects = 10), "leave them out")

  m <- usdt_mdes(p)
  expect_equal(nrow(m), 2L)
  expect_equal(m$mdes, c(pop$diff, pop$intercept), tolerance = 0.1)
})
