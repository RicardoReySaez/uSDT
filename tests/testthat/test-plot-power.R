# test-plot-power.R
# This script tests the plots of populations, power, smallest effects and
# plans: their layers and data, and what they refuse to draw.
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# This function gives the population the plots draw.
plot_population <- function() {
  usdt_population(dprime    = c(awareness = 0.15, influence = 0.30),
                  sd_dprime = c(awareness = 0.35, influence = 0.30),
                  cor_dprime = 0.3, criterion = 0.1)
}

test_that("a population plot draws the true and the measured d'", {
  pop <- plot_population()
  figure <- plot(pop, trials = c(awareness = 64, influence = 200),
                 dichotomize = "influence")
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
                          trials = c(awareness = 64, influence = 200),
                          dichotomize = "influence"))
  expect_s3_class(by_n, "ggplot")
  expect_identical(rlang::as_label(by_n$mapping$x), "subjects")
  expect_no_error(ggplot2::ggplot_build(by_n))

  by_t <- plot(usdt_power(pop, subjects = 80,
                          trials = list(awareness = c(32, 64, 128),
                                        influence = 200),
                          dichotomize = "influence"))
  expect_identical(rlang::as_label(by_t$mapping$x),
                   "trials_awareness")

  one <- plot(usdt_power(pop, subjects = 80, trials = 100, combine = "all"))
  expect_identical(rlang::as_label(one$mapping$x), "power")
  expect_true("All at once" %in% as.character(one$data$label))
  expect_no_error(ggplot2::ggplot_build(one))
  expect_error(plot(usdt_power(pop, subjects = 80, trials = 100), size = 1),
               "not used")
  # Both endings of the method leave the caption out on request.
  expect_null(plot(usdt_power(pop, subjects = 80, trials = 100),
                   caption = FALSE)$labels$caption)
  expect_null(plot(usdt_power(pop, subjects = c(40, 80), trials = 100),
                   caption = FALSE)$labels$caption)
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
                 trials = c(awareness = 64, influence = 200),
                 hypothesis = c("diff = 0", "diff in [-0.4, 0.4]"),
                 dichotomize = "influence")
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
                   trials = list(awareness = c(16, 96), influence = c(50, 300)),
                   cost = c(subject = 600, awareness = 4, influence = 1.5),
                   dichotomize = "influence")
  map <- plot(two)
  expect_true(all(map$data$relative >= 1))
  expect_equal(min(map$data$relative), 1)
  expect_gt(nrow(two$checked), 0L)
  expect_no_error(ggplot2::ggplot_build(map))

  one <- usdt_plan(pop, hypothesis = "diff = 0", subjects = c(10, 1000),
                   trials = list(awareness = c(16, 96), influence = 200),
                   cost = c(subject = 600, awareness = 4, influence = 1.5),
                   dichotomize = "influence")
  line <- plot(one)
  expect_identical(rlang::as_label(line$mapping$x),
                   "trials_awareness")
  expect_no_error(ggplot2::ggplot_build(line))

  curve <- plot(usdt_plan(pop, hypothesis = "diff = 0",
                          trials = c(awareness = 64, influence = 200),
                          dichotomize = "influence"))
  expect_true(all(diff(curve$data$y) >= -0.01))
  expect_no_error(ggplot2::ggplot_build(curve))

  width <- plot(usdt_plan(pop, precision = c(diff = 0.2),
                          trials = c(awareness = 64, influence = 200),
                          dichotomize = "influence"))
  expect_true(all(diff(width$data$y) <= 0))

  budget <- usdt_plan(pop, hypothesis = "diff = 0", subjects = c(10, 200),
                      trials = list(awareness = c(16, 96), influence = 200),
                      cost = c(subject = 600, awareness = 4, influence = 1.5),
                      budget = 50000, dichotomize = "influence")
  expect_error(plot(budget), "budget")
})
