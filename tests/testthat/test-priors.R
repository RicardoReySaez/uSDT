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

test_that("the Stan hyperparameters follow the estimated criteria and variances", {
  p <- usdt_priors(criterion = list(direct = "normal(0.1, 1)",
                                    indirect = "normal(0.2, 1)"),
                   sd_dprime = "student_t(4, 0.3, 0.2)")
  design <- function(free_c, UV) {
    list(free_c = free_c, unequal_variances = UV,
         parameterization = "correlation")
  }
  one <- .stan_priors(p, design(c(1L, 0L), FALSE))
  expect_identical(as.vector(one$mu_c_loc), 0.1)
  expect_identical(as.vector(one$sd_d_loc), c(0.3, 0.3))
  expect_length(one$sd_s_scale, 0L)

  both <- .stan_priors(p, design(c(1L, 1L), TRUE))
  expect_identical(as.vector(both$mu_c_loc), c(0.1, 0.2))
  expect_identical(as.vector(both$sd_s_scale), c(0.5, 0.5))
  expect_identical(dim(both$mu_d_loc), 2L)
})
