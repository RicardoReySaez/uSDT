test_that("hypotheses are read into a quantity, an operator and values", {
  parse <- .bf_parse
  h <- parse("intercept >= 0.1")
  expect_identical(c(h$quantity, h$op), c("intercept", ">"))
  expect_identical(h$value, 0.1)
  expect_identical(parse("diff in (0.1, -0.1)")$value, c(-0.1, 0.1))

  # A comparison of the slope with zero is tested on the correlation, which
  # has its sign and a prior of its own.
  slope <- parse("slope < 0")
  expect_identical(c(slope$quantity, slope$summary), c("rho", "slope"))
  expect_identical(parse("slope in [-1, 1]")$quantity, "slope")
  expect_identical(parse("slope > 0.2")$quantity, "slope")
  rho <- parse("rho = 0")
  expect_identical(c(rho$quantity, rho$summary), c("rho", "rho"))
  expect_false(rho$swapped)

  expect_error(parse("beta = 0"), "quantities")
  expect_error(parse("rho = 1"), "inside")
  expect_error(parse("diff >> 0"), "cannot read")
})

test_that("Bayes factors match a case with normal prior and posterior", {
  set.seed(1)
  draws <- stats::rnorm(2e4, 0.3, 0.1)
  prior <- list(log_density = function(x) stats::dnorm(x, 0, 1, log = TRUE),
                cdf = function(x) stats::pnorm(x, 0, 1), bounds = NULL,
                induced = FALSE)

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
  priors <- .quantity_priors(usdt_priors(), c("intercept", "slope"), "correlation")
  expect_equal(priors$intercept$cdf(0), 0.5, tolerance = 0.01)
  expect_identical(priors$slope$cdf(0), priors$rho$cdf(0))
})

test_that("each hypothesis prints as a report; a row subset plots its own", {
  set.seed(3)
  draws <- stats::rnorm(2e4, 0.3, 0.1)
  prior <- list(log_density = function(x) stats::dnorm(x, 0, 1, log = TRUE),
                cdf = function(x) stats::pnorm(x, 0, 1), bounds = NULL,
                induced = FALSE)
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

test_that("the priors of the tested quantities match simulated prior draws", {
  set.seed(4)
  priors <- usdt_priors(dprime = list(direct = "student_t(5, 0.5, 1)",
                                      indirect = "skew_normal(0, 0.5, 3)"),
                        sd_dprime = list(direct = "gamma(3, 6)",
                                         indirect = "student_t(4, 0, 0.5)"),
                        cor_dprime = "lkj(2)",
                        intercept = "student_t(5, 0.2, 0.4)")
  quantities <- c("diff", "rho", "intercept", "slope")
  for (parameterization in c("correlation", "regression")) {
    exact <- .quantity_priors(priors, quantities, parameterization)
    Q <- .usdt_quantities(.bivariate_primitives(
      .prior_draws(priors, 1e5, parameterization)))
    # The correlation has its own prior under both, and the intercept under
    # the latent regression only.
    own <- c("rho", if (parameterization == "regression") "intercept")
    expect_identical(names(Filter(function(p) !p$induced, exact)), own,
                     label = parameterization)
    for (q in quantities) {
      label <- paste(parameterization, q)
      at <- stats::quantile(Q[, q], c(0.1, 0.3, 0.5, 0.7, 0.9), names = FALSE)
      expect_equal(exact[[q]]$cdf(at), c(0.1, 0.3, 0.5, 0.7, 0.9),
                   tolerance = 0.01, label = label)
      # The density is the derivative of the distribution function, which
      # checks the Jacobian of a transformed prior. It is read at the 70%
      # quantile, away from the peak at zero that the slope prior has.
      h <- 1e-3 * (at[4L] - at[2L])
      slope_of_cdf <- diff(exact[[q]]$cdf(at[4L] + c(-h, h))) / (2 * h)
      expect_equal(exp(exact[[q]]$log_density(at[4L])), slope_of_cdf,
                   tolerance = 1e-3, label = label)
    }
  }
})

test_that("very heavy-tailed draws still give a posterior density", {
  # Cauchy draws like these defeat both logspline algorithms on some
  # platforms; whichever fit converges, a density must come back, and the
  # asinh fit, the last resort, must be close to the exact one.
  set.seed(17)
  x <- 0.08 + 0.1 * stats::rt(5000, 1)
  exact <- log(stats::dt(-0.8, 1) / 0.1)
  expect_true(is.finite(.bf_posterior(x)$log_density(0)))
  reserve <- .bf_posterior_asinh(x)
  expect_equal(reserve$log_density(0), exact, tolerance = 0.15)
  expect_equal(reserve$cdf(0), stats::pt(-0.8, 1), tolerance = 0.02)
})

test_that("brms priors become the priors of the tested quantities", {
  variables <- c("b_direct", "b_cond_D", "b_cond_I", "sd_subj__direct",
                 "sd_subj__cond_D", "sd_subj__cond_I",
                 "cor_subj__direct__cond_D", "cor_subj__direct__cond_I",
                 "cor_subj__cond_D__cond_I")
  table <- data.frame(
    prior = c("", "normal(0, 1)", "normal(0.5, 2)", "", "student_t(3, 0, 2.5)",
              "", "cauchy(0, 1)", "", "lkj(2)", ""),
    class = c("b", "b", "b", "b", "sd", "sd", "sd", "sd", "cor", "cor"),
    coef = c("", "", "cond_I", "cond_D", "", "", "cond_I", "cond_D", "", ""),
    group = c("", "", "", "", "", "subj", "subj", "subj", "", "subj"),
    stringsAsFactors = FALSE)[-1L, ]
  p <- .brms_priors(table, variables, "cond_D", "cond_I")

  expect_identical(.prior_rows(p, "dprime")$args, list(c(0, 1), c(0.5, 2)))
  sd <- .prior_rows(p, "sd_dprime")
  expect_identical(sd$family, c("student_t", "cauchy"))
  expect_identical(sd$args, list(c(3, 0, 2.5), c(0, 1)))
  # A 3 x 3 LKJ(2) block: one correlation has the marginal of a 2 x 2
  # LKJ(2 - 1 + 3/2), whether the prior sits on the correlations or, as in a
  # fitted model, on their Cholesky factor.
  expect_identical(.prior_rows(p, "cor_dprime", "both")$args[[1L]], 2.5)
  fitted <- table
  fitted$class[fitted$class == "cor"] <- "L"
  fitted$prior[fitted$prior == "lkj(2)"] <- "lkj_corr_cholesky(2)"
  expect_identical(.prior_rows(.brms_priors(fitted, variables, "cond_D",
                                            "cond_I"),
                               "cor_dprime", "both")$args[[1L]], 2.5)

  flat <- table
  flat$prior[flat$class == "b"] <- ""
  expect_error(.brms_priors(flat, variables, "cond_D", "cond_I"),
               "flat default")

  # An SD prior uSDT cannot simulate stops only the quantities that need it.
  circular <- table
  circular$prior[circular$class == "sd" & circular$coef == ""] <-
    c("von_mises(0, 1)", "")
  p <- .brms_priors(circular, variables, "cond_D", "cond_I")
  expect_no_error(.quantity_priors(p, c("diff", "rho"), "correlation"))
  expect_error(.quantity_priors(p, "slope", "correlation"), "cannot simulate")
})
