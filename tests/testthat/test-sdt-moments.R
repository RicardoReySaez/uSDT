# test-sdt-moments.R
# This script tests subject-level signal detection measures.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

test_that("the Hautus flag records an applied correction", {
  df <- data.frame(
    subj = rep(1:2, each = 4L),
    cond = rep(c(1L, 1L, 0L, 0L), 2L),
    resp = c(1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L))

  none <- sdt_moments(df, subject_col = "subj", condition_col = "cond",
                      condition_levels = c(signal = 1, noise = 0),
                      response_col = "resp",
                      response_levels = c(signal = 1, noise = 0),
                      correction = "none")
  hautus <- sdt_moments(df, subject_col = "subj", condition_col = "cond",
                        condition_levels = c(signal = 1, noise = 0),
                        response_col = "resp",
                        response_levels = c(signal = 1, noise = 0),
                        correction = "hautus")

  expect_false(any(none$corrected))
  expect_identical(hautus$corrected, c(TRUE, FALSE))
})
