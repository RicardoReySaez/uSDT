test_that("priors accept one value or one per task and reject other families", {
  p <- usdt_priors(dprime = list(direct = "normal(0, 1)",
                                 indirect = "normal(0.2, 0.5)"))
  indirect <- p[p$parameter == "dprime" & p$task == "indirect", ]
  expect_identical(c(indirect$location, indirect$scale), c(0.2, 0.5))
  expect_identical(p$task[p$parameter == "cor_dprime"], "both")

  expect_error(usdt_priors(sd_dprime = "normal(0, 1)"), "student_t prior")
  expect_error(usdt_priors(cor_dprime = "scaled_beta(0, 1)"), "positive")
})

test_that("the Stan hyperparameters follow the estimated criteria and variances", {
  p <- usdt_priors(criterion = list(direct = "normal(0.1, 1)",
                                    indirect = "normal(0.2, 1)"),
                   sd_dprime = "student_t(4, 0.3, 0.2)")
  one <- .stan_priors(p, free_c = c(1L, 0L), UV = 0L)
  expect_identical(as.vector(one$mu_c_loc), 0.1)
  expect_identical(as.vector(one$sd_d_loc), c(0.3, 0.3))
  expect_length(one$sd_s_scale, 0L)

  both <- .stan_priors(p, free_c = c(1L, 1L), UV = 1L)
  expect_identical(as.vector(both$mu_c_loc), c(0.1, 0.2))
  expect_identical(as.vector(both$sd_s_scale), c(0.5, 0.5))
  expect_identical(dim(both$mu_d_loc), 2L)
})
