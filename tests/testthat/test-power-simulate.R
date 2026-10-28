# test-power-simulate.R
# This script tests the simulated power of usdt_power(): the simulated
# studies, the decisions on each one, the random number streams, the failures,
# the precision summary and the parallel run.
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# This function gives the scenario and the design most tests simulate.
simulation_case <- function(...) {
  pop <- usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
                         sd_dprime = c(direct = 0.35, indirect = 0.30),
                         cor_dprime = 0.3, criterion = 0.1,
                         sd_criterion = 0.3, ...)
  list(population = pop,
       scenario = as.list(as.data.frame(pop)[1L, ]),
       cells = list(signal = c(direct = 32, indirect = 60),
                    noise = c(direct = 32, indirect = 60)),
       split = c(direct = FALSE, indirect = TRUE),
       labels = c(direct = "Direct", indirect = "Indirect"))
}

test_that("a simulated study has the layout of usdt_data_tasks()", {
  case <- simulation_case()
  set.seed(4)
  d <- .simulate_usdt_data(case$scenario, 30L, case$cells, case$split,
                           case$labels)
  expect_s3_class(d, "usdt_data")
  expect_equal(nrow(d$agg), 4L * 30L)
  expect_named(d$agg, c("subj", "task", "sig", "y", "n", "cond", "c_D", "c_I",
                        "d_D", "d_I"))
  expect_equal(d$agg$n[d$agg$task == "D"], rep(32L, 60L))

  # The median split gives every subject half signal responses, so the
  # criterion of that task is fixed and the other is estimated.
  split <- d$agg[d$agg$task == "I", ]
  expect_equal(as.vector(tapply(split$y, split$subj, sum)), rep(60, 30L))
  expect_identical(d$meta$criterion_zero, c(direct = FALSE, indirect = TRUE))
  expect_output(print(d), "simulated study")
  expect_identical(hsdt(d)$design$criteria, "c_D")
})

test_that("simulated studies recover the population", {
  case <- simulation_case()
  set.seed(8)
  d <- .simulate_usdt_data(case$scenario, 400L,
                           list(signal = c(direct = 60, indirect = 150),
                                noise = c(direct = 60, indirect = 150)),
                           case$split, case$labels)
  # With 400 subjects the means and SDs have standard errors near 0.02 and
  # the correlation near 0.1, so each must lie within three of them.
  e <- hsdt(d)$pars$est
  expect_true(all(abs(e[c("gamma_D", "gamma_I")] - c(0.15, 0.30)) < 0.07))
  expect_true(all(abs(sqrt(e[c("s2_D", "s2_I")]) - c(0.35, 0.30)) < 0.07))
  expect_lt(abs(e[["s_DI"]] / sqrt(e[["s2_D"]] * e[["s2_I"]]) - 0.3), 0.3)
})

test_that("each simulated study is decided as usdt_tests() and usdt_equivalence() decide it", {
  case <- simulation_case()
  hypotheses <- c("diff = 0", "rho = 0", "intercept > 0.1", "diff in [-0.3, 0.3]")
  tests <- lapply(hypotheses, .parse_hypothesis)
  streams <- .power_streams(3L, seed = 11)
  for (s in streams) {
    out <- .power_replicate(s, case$scenario, 40L, case$cells, case$split,
                            case$labels, tests, alpha = 0.05)

    # The same stream gives the same study, analysed here directly.
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
    assign(".Random.seed", s, envir = globalenv())
    d <- .simulate_usdt_data(case$scenario, 40L, case$cells, case$split,
                             case$labels)
    if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, envir = globalenv())
    m <- suppressWarnings(hsdt(d))
    t <- usdt_tests(m)
    e <- suppressWarnings(usdt_equivalence(m, hypotheses[3:4]))

    expect_equal(out[["est_diff"]], t$estimate[1L])
    expect_equal(out[["reject_1"]], as.numeric(t$p.value[1L] < 0.05))
    expect_equal(out[["reject_2"]], as.numeric(isTRUE(t$p.value[2L] < 0.05)))
    expect_equal(out[["reject_3"]], as.numeric(identical(e$decision[1L],
                                                         "reject H0")))
    expect_equal(out[["reject_4"]], as.numeric(identical(e$decision[2L],
                                                         "reject H0")))
  }
})

test_that("a simulation reports power with its Monte Carlo error and precision", {
  case <- simulation_case()
  p <- usdt_power(case$population, subjects = 30, trials = c(direct = 64, indirect = 120),
                  dichotomize = "indirect", method = "simulation", nsim = 6,
                  seed = 2, combine = "all")
  expect_s3_class(p, "usdt_power")
  expect_true(all(c("mc_se", "power.low", "power.high", "failures") %in%
                    names(p)))
  expect_equal(p$mc_se, sqrt(p$power * (1 - p$power) / 6))
  expect_true(all(p$power.low <= p$power & p$power <= p$power.high))
  expect_output(print(p), "simulated studies")

  s <- summary(p)
  expect_s3_class(s, "summary.usdt_power")
  expect_setequal(s$quantity, c("diff", "rho", "intercept", "slope"))
  expect_equal(s$true[s$quantity == "rho"], 0.3)
  expect_output(print(s), "Coverage")

  # The same seed gives the same studies, without touching the session's
  # random numbers.
  set.seed(1)
  before <- stats::runif(1L)
  set.seed(1)
  again <- usdt_power(case$population, subjects = 30, trials = c(direct = 64, indirect = 120),
                      dichotomize = "indirect", method = "simulation",
                      nsim = 6, seed = 2, combine = "all")
  expect_identical(stats::runif(1L), before)
  expect_identical(again$power, p$power)
})

test_that("failed fits count as failures to reject", {
  case <- simulation_case()
  local_mocked_bindings(.power_fit = function(data, level) stop("no fit"))
  p <- usdt_power(case$population, subjects = 30, trials = c(direct = 64, indirect = 120),
                  dichotomize = "indirect", method = "simulation", nsim = 4,
                  seed = 1)
  expect_true(all(p$power == 0))
  expect_true(all(p$estimable == 0))
  expect_true(all(p$failures == 4L))
  expect_output(print(p), "did not converge")
})

test_that("a population with unequal variances can only be simulated", {
  case <- simulation_case(sd_signal = c(direct = 1, indirect = 1.3))
  expect_error(usdt_power(case$population, subjects = 30, trials = 100),
               "method = \"simulation\"")
  local_mocked_bindings(.power_fit = function(data, level) stop("no fit"))
  p <- usdt_power(case$population, subjects = 30, trials = 100,
                  dichotomize = "indirect", method = "simulation", nsim = 2,
                  seed = 1)
  expect_true(all(is.na(p$se)))
  expect_true(all(is.na(p$reliability_direct)))
})

test_that("the summary needs a simulation, and the arguments are checked", {
  case <- simulation_case()
  expect_error(summary(usdt_power(case$population, subjects = 30,
                                  trials = 100)), "method = \"simulation\"")
  run <- function(...) {
    usdt_power(case$population, subjects = 30, trials = 100,
               method = "simulation", ...)
  }
  expect_error(run(nsim = 0), "nsim")
  expect_error(run(nsim = 2.5), "nsim")
  expect_error(run(ncores = 0), "ncores")
  expect_error(run(seed = -1), "seed")
  expect_error(run(progress = NA), "progress")
})

test_that("a parallel run gives the studies of a sequential one", {
  skip_on_cran()
  case <- simulation_case()
  run <- function(ncores) {
    usdt_power(case$population, subjects = 30, trials = c(direct = 64, indirect = 120),
               dichotomize = "indirect", method = "simulation", nsim = 4,
               seed = 5, ncores = ncores)
  }
  expect_identical(run(2L)$power, run(1L)$power)
})
