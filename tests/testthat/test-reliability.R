# test-reliability.R
# This script tests reliability estimates for both tasks.
# Author: Ricardo Rey-Sáez
# Last modified: 04-09-2026

# Both criteria are estimated here.
fitted_example <- function() {
  set.seed(21)
  trials <- make_trials(n_subj = 24L, n_trials = 60L)
  data <- usdt_data_long(
    trials,
    task_col = "task",
    task_levels = c(direct = "D", indirect = "I"),
    subject_col = "subj",
    condition_col = "cond",
    condition_levels = c(signal = 1, noise = 0),
    response_col = "response",
    response_levels = c(signal = 1, noise = 0)
  )
  list(data = data, fit = suppressWarnings(hsdt(data)))
}

# A Meyen split removes the indirect criterion.
split_example <- function() {
  set.seed(21)
  trials <- make_trials(n_subj = 24L, n_trials = 60L, rt = TRUE)
  data <- suppressMessages(usdt_data_long(
    trials,
    task_col = "task",
    task_levels = c(direct = "D", indirect = "I"),
    subject_col = "subj",
    condition_col = "cond",
    condition_levels = c(signal = 1, noise = 0),
    response_col = c(direct = "response", indirect = "rt"),
    response_levels = list(direct   = c(signal = 1, noise = 0),
                           indirect = c(signal = "faster", noise = "slower")),
    dichotomize = "indirect"
  ))
  list(data = data, fit = suppressWarnings(hsdt(data)))
}

# This function names the design columns of one task.
task_columns <- function(fit, task) {
  c(intersect(paste0("c_", task), fit$design$criteria), paste0("d_", task))
}

test_that("reliability follows the documented formulas", {
  example <- fitted_example()
  result <- usdt_reliability(example$fit)

  expect_s3_class(result, "usdt_reliability")
  expect_equal(nrow(result$tasks), 2L)
  expect_equal(nrow(result$subjects), 48L)

  for (task in result$tasks$task) {
    row <- result$tasks[result$tasks$task == task, ]
    subject <- result$subjects[result$subjects$task == task, ]

    # Group reliability uses the mean subject variance.
    expect_equal(row$reliability, row$tau2 / (row$tau2 + row$mean_variance))
    expect_equal(row$mean_variance, mean(subject$variance))
    expect_equal(subject$reliability, row$tau2 / (row$tau2 + subject$variance))
  }

  # A variance ratio cannot leave the unit interval.
  expect_true(all(result$subjects$reliability > 0 &
                    result$subjects$reliability < 1))
})

test_that("reliability requires a fitted model", {
  example <- fitted_example()
  expect_error(usdt_reliability(example$data), "hsdt")
})

test_that("reliability uses every usable bootstrap refit", {
  example <- fitted_example()
  point <- usdt_reliability(example$fit)
  statistic <- .boot_statistic(example$fit$fit, "d_D", "d_I")
  fitted <- statistic(example$fit$fit)
  samples <- matrix(rep(fitted, times = .boot_min),
                    nrow = .boot_min, byrow = TRUE,
                    dimnames = list(NULL, names(fitted)))
  effects <- .boot_effects(fitted, samples, attr(statistic, "layout"))
  example$fit$boot <- list(
    variance = effects$variance,
    population = effects$population,
    subjects = effects$subjects,
    ok = rep(TRUE, .boot_min), summary_available = TRUE,
    complete = TRUE, usable = .boot_min, level = 0.95, type = "perc"
  )

  result <- usdt_reliability(example$fit)
  output <- capture.output(summary(result))

  expect_equal(result$tasks$boot_median, point$tasks$reliability)
  expect_equal(result$tasks$boot_se, c(0, 0), tolerance = 1e-12)
  expect_equal(result$subjects$boot_median, point$subjects$reliability)
  expect_identical(dim(result$boot$samples$tasks), c(.boot_min, 2L))
  expect_true(any(grepl("Group boot median", output, fixed = TRUE)))
  expect_false(any(grepl("Components", output, fixed = TRUE)))
})

test_that("the profiled information reproduces the closed form", {
  for (example in list(fitted_example(), split_example())) {
    fit <- example$fit
    agg <- example$data$agg
    fixed <- lme4::fixef(fit$fit)
    random <- lme4::ranef(fit$fit, condVar = FALSE)[["subj"]]
    subject <- rownames(random)[1L]

    for (task in c("D", "I")) {
      slope <- paste0("d_", task)
      columns <- task_columns(fit, task)
      cells <- agg[as.character(agg$subj) == subject & agg$task == task, ]
      coef <- fixed[columns] +
        unlist(random[subject, columns, drop = FALSE], use.names = FALSE)
      eta <- drop(as.matrix(cells[, columns, drop = FALSE]) %*% coef)

      # Rebuild the two cell weights by hand.
      p <- stats::pnorm(eta)
      w <- cells$n * stats::dnorm(eta)^2 / (p * (1 - p))

      # The calculation profiles an estimated criterion.
      expected <- if (paste0("c_", task) %in% fit$design$criteria) {
        sum(1 / w)
      } else {
        4 / sum(w)
      }
      expect_equal(.sensitivity_variance(cells, columns, slope, eta), expected)
    }
  }
})

test_that("the profiled information is the error variance of a per-subject GLM", {

  # A separate probit model provides the reference variance.
  for (example in list(fitted_example(), split_example())) {
    fit <- example$fit
    agg <- example$data$agg

    for (task in c("D", "I")) {
      slope <- paste0("d_", task)
      columns <- task_columns(fit, task)
      formula <- stats::as.formula(paste(
        "cbind(y, n - y) ~",
        if (paste0("c_", task) %in% fit$design$criteria) slope
        else paste("0 +", slope)))

      rows <- agg[agg$task == task, ]
      ours <- theirs <- numeric(0)
      for (s in unique(as.character(rows$subj))) {
        cells <- rows[as.character(rows$subj) == s, , drop = FALSE]

        # A cell at floor or ceiling separates the per-subject fit.
        if (any(cells$y == 0 | cells$y == cells$n)) next

        # A tight fit aligns the final model weights.
        one <- stats::glm(formula, family = stats::binomial("probit"),
                          data = cells,
                          control = stats::glm.control(epsilon = 1e-12,
                                                       maxit = 200L))
        eta <- unname(stats::predict(one, type = "link"))
        ours <- c(ours, .sensitivity_variance(cells, columns, slope, eta))
        theirs <- c(theirs, unname(stats::vcov(one)[slope, slope]))
      }

      # The remaining tolerance reflects the final model iteration.
      expect_gt(length(ours), 0L)
      expect_equal(ours, theirs, tolerance = 1e-6)
    }
  }
})

test_that("the profiled information is var_gg at the observed rates", {

  # Both routes evaluate the same variance at the same rates.
  example <- fitted_example()
  moments <- sdt_moments(example$data, correction = "hautus", variances = TRUE)
  moments <- moments[moments$task == "Direct", ]

  rows <- example$data$agg[example$data$agg$task == "D", ]
  ours <- vapply(as.character(moments$subj), function(s) {
    cells <- rows[as.character(rows$subj) == s, , drop = FALSE]
    cells <- cells[order(cells$sig), ]
    row <- moments[as.character(moments$subj) == s, ]
    .sensitivity_variance(cells, c("c_D", "d_D"), "d_D",
                          stats::qnorm(c(row$far, row$hr)))
  }, 0)

  expect_equal(unname(ours), moments$var_gg)
})

test_that("Miller variance uses the observed number of trials", {
  result <- .sdt_table(
    subj = "1", hit = 2, miss = 0, fa = 0, cr = 2,
    coding = "deviation", correction = "hautus", variances = TRUE
  )
  expected <- .sdt_variances(
    nr = 2, ns = 2, pi_fa = result$far, pi_a = result$hr
  )

  expect_equal(result$var_miller, expected$var_miller)
})
