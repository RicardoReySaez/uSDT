test_that("priors take one value or one per task from the families that fit", {
  p <- usdt_priors(dprime = list(awareness = "normal(0, 1)",
                                 influence = "skew_normal(0.2, 0.5, 2)"),
                   sd_dprime = "gamma(2, 4)")
  influence <- .prior_rows(p, "dprime", "influence")
  expect_identical(influence$family, "skew_normal")
  expect_identical(influence$args[[1L]], c(0.2, 0.5, 2))
  expect_identical(p$task[p$parameter == "cor_dprime"], "both")

  expect_error(usdt_priors(sd_dprime = "beta(2, 2)"), "any positive value")
  expect_error(usdt_priors(dprime = "lognormal(0, 1)"), "any real value")
  expect_error(usdt_priors(cor_dprime = "normal(0, 1)"), "lkj\\(eta\\)")
  expect_error(usdt_priors(cor_dprime = "lkj(0)"), "eta > 0")
  expect_error(usdt_priors(dprime = "von_mises(0, 1)"), "`student_t`")
})

test_that("median-split prior predictions balance hits and false alarms", {
  set.seed(31)
  for (latent in c("correlation", "regression")) {
    for (unequal in c(FALSE, TRUE)) {
      predictions <- .prior_predictive(usdt_priors(), 2000,
        list(free_c = c(0L, 0L), unequal_variances = unequal, latent = latent))
      expect_equal(predictions$hit + predictions$fa, matrix(1, 2000, 2))
      expect_true(all(is.finite(predictions$hit)))
    }
  }
})

test_that("prior plots show the stated and induced hypothesis distributions", {
  set.seed(32)
  priors <- usdt_priors(intercept = "normal(0, 0.3)")
  p <- plot(priors)
  curves <- p$data
  difference <- curves[curves$panel == "H1: difference in mean d'", ]
  correlation <- curves[curves$panel == "H2: correlation between d'", ]
  expect_lt(max(abs(difference$density -
                      dnorm(difference$x, sd = sqrt(2)))), 0.005)
  expect_equal(correlation$density, rep(0.5, nrow(correlation)))
  rates <- curves[!is.na(curves$rate), ]
  for (r in split(rates, interaction(rates$task, rates$rate))) {
    mass <- sum(diff(r$x) * (head(r$density, -1) + tail(r$density, -1)) / 2)
    expect_equal(mass, 1, tolerance = 0.02)
  }

  regression <- plot(priors, latent = "regression", caption = FALSE)
  intercept <- regression$data[
    regression$data$panel == "H3: intercept of the latent regression", ]
  expect_equal(intercept$density, dnorm(intercept$x, sd = 0.3))
  expect_null(ggplot2::get_labs(regression)$caption)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(print(p))
  expect_no_error(print(regression))
})
