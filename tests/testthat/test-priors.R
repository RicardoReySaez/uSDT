test_that("priors take one value or one per task from the families that fit", {
  p <- usdt_priors(dprime = list(direct = "normal(0, 1)",
                                 indirect = "skew_normal(0.2, 0.5, 2)"),
                   sd_dprime = "gamma(2, 4)")
  indirect <- .prior_rows(p, "dprime", "indirect")
  expect_identical(indirect$family, "skew_normal")
  expect_identical(indirect$args[[1L]], c(0.2, 0.5, 2))
  expect_identical(p$task[p$parameter == "cor_dprime"], "both")

  expect_error(usdt_priors(sd_dprime = "beta(2, 2)"), "any positive value")
  expect_error(usdt_priors(dprime = "lognormal(0, 1)"), "any real value")
  expect_error(usdt_priors(cor_dprime = "normal(0, 1)"), "lkj\\(eta\\)")
  expect_error(usdt_priors(cor_dprime = "lkj(0)"), "eta > 0")
  expect_error(usdt_priors(dprime = "von_mises(0, 1)"), "`student_t`")
})
