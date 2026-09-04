# test-meyen-split.R
# This script tests the median split used for response times.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

test_that("meyen_split reproduces the published idiom", {
  set.seed(1)
  rt   <- stats::rlnorm(400, 6, .3)
  subj <- rep(1:4, each = 100)

  # This expression reproduces the split used in the literature.
  ref <- stats::ave(rt, subj, FUN = function(r) as.integer(stats::median(r) > r))
  expect_identical(meyen_split(rt, by = subj), as.integer(ref))
})

test_that("an even, balanced split puts the criterion at exactly zero", {
  set.seed(2)
  n    <- 200L
  rt   <- stats::rlnorm(n, 6, .3)
  sig  <- rep(c(0L, 1L), each = n / 2L)
  resp <- meyen_split(rt, by = rep("s1", n))

  # Half of the trials become signal responses.
  expect_equal(mean(resp), 0.5)

  # Balanced conditions make the deviation-coded criterion equal to zero.
  hr  <- mean(resp[sig == 1L])
  far <- mean(resp[sig == 0L])
  expect_equal((stats::qnorm(hr) + stats::qnorm(far)) / 2, 0, tolerance = 1e-12)
})

test_that("the split direction can be reversed", {
  set.seed(3)
  rt <- stats::rlnorm(100, 6, .3)
  g  <- rep("s1", 100)
  expect_identical(meyen_split(rt, g, signal = "faster"),
                   1L - meyen_split(rt, g, signal = "slower"))
})

test_that("ties are handled as documented", {
  # An odd count makes the median equal to an observed value.
  x <- c(1, 2, 3, 4, 5)
  g <- rep("s1", 5)

  # The noise option sends the tied trial to the noise side.
  expect_identical(meyen_split(x, g, ties = "noise"), c(1L, 1L, 0L, 0L, 0L))

  # The random option resolves ties without favoring either side.
  set.seed(4)
  p <- replicate(2000, mean(meyen_split(x, g, ties = "random")))
  expect_setequal(unique(p), c(0.4, 0.6))
  expect_equal(mean(p), 0.5, tolerance = 0.03)
})

test_that("invalid input is rejected with a useful message", {
  expect_error(meyen_split(c("a", "b"), c(1, 1)), "must be numeric")
  expect_error(meyen_split(1:10, 1:3), "length")
  expect_error(meyen_split(rep(NA_real_, 5), rep(1, 5)), "entirely missing")
})

test_that("missing measures propagate", {
  x <- c(1, 2, NA, 4, 5, 6)
  expect_true(is.na(meyen_split(x, rep("s1", 6))[3L]))
})
