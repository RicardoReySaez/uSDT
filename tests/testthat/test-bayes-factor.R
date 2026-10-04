test_that("hypotheses are read into a quantity, an operator and values", {
  h <- .bf_parse("intercept >= 0.1")
  expect_identical(c(h$quantity, h$op), c("intercept", ">"))
  expect_identical(h$value, 0.1)
  expect_identical(.bf_parse("diff in (0.1, -0.1)")$value, c(-0.1, 0.1))

  slope <- .bf_parse("slope < 0")
  expect_identical(c(slope$quantity, slope$summary), c("rho", "slope"))
  expect_identical(.bf_parse("slope in [-1, 1]")$quantity, "slope")

  expect_error(.bf_parse("beta = 0"), "quantities")
  expect_error(.bf_parse("rho = 1"), "inside")
  expect_error(.bf_parse("diff >> 0"), "cannot read")
})

test_that("Bayes factors match a case with normal prior and posterior", {
  skip_if_not_installed("logspline")
  set.seed(1)
  draws <- stats::rnorm(2e4, 0.3, 0.1)
  prior <- list(log_density = function(x) stats::dnorm(x, 0, 1, log = TRUE),
                cdf = function(x) stats::pnorm(x, 0, 1), bounds = NULL)

  point <- .bf_one(.bf_parse("diff = 0"), draws, draws, prior, 0.95)$row
  expected <- stats::dnorm(0, 0, 1, log = TRUE) -
    stats::dnorm(0, 0.3, 0.1, log = TRUE)
  expect_equal(point$log_BF10, expected, tolerance = 0.05)

  directional <- .bf_one(.bf_parse("diff > 0"), draws, draws, prior, 0.95)$row
  expect_equal(directional$post.prob, stats::pnorm(3), tolerance = 0.005)
  expect_equal(directional$log_BF10, stats::qlogis(stats::pnorm(3)),
               tolerance = 0.05)
})

test_that("the induced priors keep the symmetry of the correlation prior", {
  set.seed(2)
  priors <- .bf_priors(usdt_priors(), c("intercept", "slope"))
  expect_equal(priors$intercept$cdf(0), 0.5, tolerance = 0.01)
  expect_identical(priors$slope$cdf(0), priors$rho$cdf(0))
})
