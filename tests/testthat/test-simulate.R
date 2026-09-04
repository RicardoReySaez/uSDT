# test-simulate.R
# This script tests the input checks of the simulator.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

test_that("the simulator rejects undefined model parameters", {
  expect_error(usdt_simulate(sd_D = 0), "invalid value")
  expect_error(usdt_simulate(sd_I = 0), "invalid value")
  expect_error(usdt_simulate(rho = 1), "invalid value")
  expect_error(usdt_simulate(sd_crit_I = -1), "invalid value")
  expect_error(usdt_simulate(n_subj = 2.5), "invalid value")
  expect_error(usdt_simulate(n_trials = 3), "must be even")
  expect_error(usdt_simulate(gamma_D = Inf), "invalid value")
  expect_error(usdt_simulate(rt = 1), "TRUE.*FALSE")
})
