# test-plot-power.R
# This script tests the plots of populations, power, smallest effects and
# plans: their layers and data, and what they refuse to draw.
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# This function gives the population the plots draw.
plot_population <- function() {
  usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
                  sd_dprime = c(direct = 0.35, indirect = 0.30),
                  cor_dprime = 0.3, criterion = 0.1)
}

test_that("a population plot draws the true and the measured d'", {
  pop <- plot_population()
  figure <- plot(pop, trials = c(direct = 64, indirect = 200),
                 dichotomize = "indirect")
  expect_s3_class(figure, "ggplot")
  expect_setequal(unique(figure$data$kind), c("true", "measured"))
  expect_no_error(ggplot2::ggplot_build(figure))

  # The 90% ellipse of the true d' holds 90% of a bivariate normal.
  true <- figure$data[figure$data$kind == "true" &
                        figure$data$coverage == "0.9", ]
  S <- matrix(c(0.35^2, 0.3 * 0.35 * 0.30, 0.3 * 0.35 * 0.30, 0.30^2), 2L)
  d <- cbind(true$x - 0.15, true$y - 0.30)
  expect_equal(rowSums((d %*% solve(S)) * d), rep(stats::qchisq(0.9, 2L),
                                                  nrow(d)), tolerance = 1e-8)

  several <- plot(update(pop, cor_dprime = c(0.1, 0.5)))
  expect_identical(levels(several$data$scenario),
                   c("cor(d') 0.10", "cor(d') 0.50"))
  expect_error(plot(update(pop, cor_dprime = seq(-0.6, 0.7, 0.1))),
               "at most 12")
  expect_error(plot(pop, trials = c(64, 128)), "one design")
})

test_that("a power plot puts subjects or trials on the axis", {
  pop <- plot_population()
  by_n <- plot(usdt_power(pop, subjects = c(40, 80, 120),
                          trials = c(direct = 64, indirect = 200),
                          dichotomize = "indirect"))
  expect_s3_class(by_n, "ggplot")
  expect_identical(rlang::as_label(by_n$mapping$x), "subjects")
  expect_no_error(ggplot2::ggplot_build(by_n))

  by_t <- plot(usdt_power(pop, subjects = 80,
                          trials = list(direct = c(32, 64, 128), indirect = 200),
                          dichotomize = "indirect"))
  expect_identical(rlang::as_label(by_t$mapping$x),
                   "trials_direct")

  one <- plot(usdt_power(pop, subjects = 80, trials = 100, combine = "all"))
  expect_identical(rlang::as_label(one$mapping$x), "power")
  expect_true("All at once" %in% as.character(one$data$label))
  expect_no_error(ggplot2::ggplot_build(one))
  expect_error(plot(usdt_power(pop, subjects = 80, trials = 100), 1),
               "not used")
})

test_that("a simulated power plot adds its intervals", {
  pop <- plot_population()
  local_mocked_bindings(.power_fit = function(data, level) stop("no fit"))
  sims <- usdt_power(pop, subjects = c(30, 40), trials = 100,
                     method = "simulation", nsim = 2, seed = 1)
  figure <- plot(sims)
  geoms <- vapply(figure$layers, function(l) class(l$geom)[1L], "")
  expect_true("GeomErrorbar" %in% geoms)
  expect_no_error(ggplot2::ggplot_build(figure))
})

test_that("an MDES plot draws the curves and marks the smallest effects", {
  m <- usdt_mdes(plot_population(), subjects = 80,
                 trials = c(direct = 64, indirect = 200),
                 hypothesis = c("diff = 0", "diff in [-0.4, 0.4]"),
                 dichotomize = "indirect")
  figure <- plot(m)
  expect_s3_class(figure, "ggplot")
  expect_identical(levels(figure$data$panel),
                   c("Delta d' = 0", "Delta d' in [-0.4, 0.4]"))
  expect_equal(nrow(figure$data), 2L * 25L)
  expect_no_error(ggplot2::ggplot_build(figure))
})

test_that("a plan plot draws costs over the trials, or power along one value", {
  pop <- plot_population()
  two <- usdt_plan(pop, hypothesis = "diff = 0", subjects = c(10, 1000),
                   trials = list(direct = c(16, 96), indirect = c(50, 300)),
                   cost = c(subject = 600, direct = 4, indirect = 1.5),
                   dichotomize = "indirect")
  map <- plot(two)
  expect_true(all(map$data$relative >= 1))
  expect_equal(min(map$data$relative), 1)
  expect_gt(nrow(two$checked), 0L)
  expect_no_error(ggplot2::ggplot_build(map))

  one <- usdt_plan(pop, hypothesis = "diff = 0", subjects = c(10, 1000),
                   trials = list(direct = c(16, 96), indirect = 200),
                   cost = c(subject = 600, direct = 4, indirect = 1.5),
                   dichotomize = "indirect")
  line <- plot(one)
  expect_identical(rlang::as_label(line$mapping$x),
                   "trials_direct")
  expect_no_error(ggplot2::ggplot_build(line))

  curve <- plot(usdt_plan(pop, hypothesis = "diff = 0",
                          trials = c(direct = 64, indirect = 200),
                          dichotomize = "indirect"))
  expect_true(all(diff(curve$data$y) >= -0.01))
  expect_no_error(ggplot2::ggplot_build(curve))

  width <- plot(usdt_plan(pop, precision = c(diff = 0.2),
                          trials = c(direct = 64, indirect = 200),
                          dichotomize = "indirect"))
  expect_true(all(diff(width$data$y) <= 0))

  budget <- usdt_plan(pop, hypothesis = "diff = 0", subjects = c(10, 200),
                      trials = list(direct = c(16, 96), indirect = 200),
                      cost = c(subject = 600, direct = 4, indirect = 1.5),
                      budget = 50000, dichotomize = "indirect")
  expect_error(plot(budget), "budget")
})
