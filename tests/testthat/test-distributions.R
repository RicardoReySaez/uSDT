test_that("each prior family integrates to one, and its CDF and draws agree", {
  examples <- list(
    normal = c(0.3, 1.2), std_normal = numeric(0), student_t = c(4, 0.2, 0.7),
    cauchy = c(-0.1, 0.5), double_exponential = c(0.2, 0.8),
    logistic = c(0.1, 0.6), gumbel = c(0.3, 0.9),
    skew_normal = c(0.1, 1.3, 2.5), exp_mod_normal = c(0.2, 0.7, 1.5),
    skew_double_exponential = c(0.1, 0.9, 0.3), lognormal = c(-0.5, 0.6),
    exponential = 2.5, gamma = c(2, 3), inv_gamma = c(3, 1.5),
    weibull = c(1.7, 0.8), frechet = c(2.5, 0.7), rayleigh = 0.6,
    chi_square = 3, inv_chi_square = 4, scaled_inv_chi_square = c(5, 0.4),
    loglogistic = c(0.8, 3), pareto = c(0.2, 2.5),
    pareto_type_2 = c(-0.3, 1.2, 3), uniform = c(-0.5, 1.5), beta = c(2, 3),
    lkj = 2.5)
  expect_setequal(names(examples), names(.distributions))

  set.seed(1)
  for (f in names(examples)) {
    family <- .distributions[[f]]
    a <- examples[[f]]
    support <- family$support(a)
    density <- function(x) exp(family$logd(x, a))
    draws <- family$r(2e4, a)
    at <- stats::quantile(draws, c(0.2, 0.5, 0.8), names = FALSE)
    expect_equal(stats::integrate(density, support[1L], support[2L])$value, 1,
                 tolerance = 1e-5, label = f)
    expect_equal(family$p(at, a), vapply(at, function(v) {
      stats::integrate(density, support[1L], v)$value
    }, 0), tolerance = 1e-5, label = f)
    # runif() has 32-bit resolution, so 20,000 draws can repeat a value; one
    # tie does not affect the test.
    ks <- suppressWarnings(stats::ks.test(draws, function(q) family$p(q, a)))
    expect_gt(ks$p.value, 1e-3, label = f)
  }
})

test_that("a standard deviation prior on the real line is truncated at zero", {
  row <- .prior_rows(usdt_priors(sd_dprime = "normal(-0.2, 0.5)"),
                     "sd_dprime")[1L, ]
  prior <- .prior_fun(row)
  set.seed(2)
  draws <- prior$r(2e4)
  expect_gte(min(draws), 0)
  expect_identical(prior$cdf(0), 0)
  expect_equal(stats::integrate(function(x) exp(prior$log_density(x)), 0,
                                Inf)$value, 1, tolerance = 1e-6)
  expect_gt(stats::ks.test(draws, prior$cdf)$p.value, 1e-3)
})
