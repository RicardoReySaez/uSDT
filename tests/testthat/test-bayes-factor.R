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

  # The tested value lies one posterior SD from the mean. Far in a tail, the
  # logspline density changes with the draws and with the platform: at three
  # SDs its log BF10 is off by up to 0.2 (point) and 0.45 (directional).
  point <- .bf_one(.bf_parse("diff = 0.2"), draws, draws, prior, 0.95)$row
  expected <- stats::dnorm(0.2, 0, 1, log = TRUE) -
    stats::dnorm(0.2, 0.3, 0.1, log = TRUE)
  expect_equal(point$log_BF10, expected, tolerance = 0.05)

  directional <- .bf_one(.bf_parse("diff > 0.2"), draws, draws, prior,
                         0.95)$row
  expect_equal(directional$post.prob, stats::pnorm(1), tolerance = 0.02)
  odds <- stats::qlogis(stats::pnorm(1)) - stats::qlogis(stats::pnorm(-0.2))
  expect_equal(directional$log_BF10, odds, tolerance = 0.05)
})

test_that("the induced priors keep the symmetry of the correlation prior", {
  set.seed(2)
  priors <- .bf_priors(usdt_priors(), c("intercept", "slope"))
  expect_equal(priors$intercept$cdf(0), 0.5, tolerance = 0.01)
  expect_identical(priors$slope$cdf(0), priors$rho$cdf(0))
})

test_that("each hypothesis prints as a report; a row subset plots its own", {
  skip_if_not_installed("logspline")
  set.seed(3)
  draws <- stats::rnorm(2e4, 0.3, 0.1)
  prior <- list(log_density = function(x) stats::dnorm(x, 0, 1, log = TRUE),
                cdf = function(x) stats::pnorm(x, 0, 1), bounds = NULL)
  results <- lapply(c("diff = 0", "diff > 0"), function(text) {
    .bf_one(.bf_parse(text), draws, draws, prior, 0.95)
  })
  b <- structure(do.call(rbind, lapply(results, `[[`, "row")),
                 curves = lapply(results, `[[`, "curve"), level = 0.95,
                 class = c("usdt_bf", "data.frame"))
  rownames(b) <- NULL

  out <- utils::capture.output(print(b))
  expect_length(grep("Savage-Dickey density ratio test", out, fixed = TRUE), 1L)
  expect_length(grep("Posterior probability (H1)", out, fixed = TRUE), 1L)
  expect_length(grep("95% CrI", out, fixed = TRUE), 2L)
  expect_output(print(b[, c("hypothesis", "log_BF10")]), "d' > 0",
                fixed = TRUE)

  # Panel titles are plotmath, so each one must parse.
  panel <- levels(plot(b[2, ])$data$panel)
  expect_identical(panel,
                   "bold(Delta*d*\"'\" > \"0\" ~~ \"(directional)\")")
  expect_no_error(parse(text = panel))
  expect_error(plot(b[, 1:3]), "lost")
})

test_that("every kind of hypothesis gives a panel title that parses", {
  for (text in c("diff = 0", "rho < -0.2", "slope > 0", "intercept in [-1, 1]",
                 "diff out [-0.1, 0.1]")) {
    expect_no_error(parse(text = .bf_math(.bf_parse(text), "test")))
  }
})
