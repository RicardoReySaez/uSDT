# test-usdt-freq.R
# This script tests the frequentist model and its hypotheses.
# Author: Ricardo Rey-Sáez
# Last modified: 07-09-2026

# This function fits a stable model for the comparisons below.
fit_reference <- function(seed = 4L, ...) {
  set.seed(seed)
  df <- usdt_simulate(n_subj = 60L, n_trials = 120L, gamma_D = 0.8,
                      gamma_I = 0.35, sd_D = 0.5, sd_I = 0.3, rho = 0.5)
  d  <- usdt_data_long(df, task_col = "task",
                       task_levels      = c(direct = "D", indirect = "I"),
                       subject_col      = "subj",
                       condition_col    = "cond",
                       condition_levels = c(signal = 1, noise = 0),
                       response_col     = "response",
                       response_levels  = c(signal = 1, noise = 0))
  hsdt(d, ...)
}

test_that("the joint covariance uses the fitted Hessian", {
  m <- fit_reference()
  H <- m$fit@optinfo$derivs$Hessian
  expected <- solve((H + t(H)) / 4)
  dimnames(expected) <- dimnames(m$pars$V_full)

  expect_identical(attr(m$pars$V_full, "source"), "lme4")
  expect_equal(m$pars$V_full, expected, tolerance = 1e-12,
               ignore_attr = TRUE)
  fixed <- names(lme4::fixef(m$fit))
  expect_equal(m$pars$V_full[fixed, fixed],
               as.matrix(stats::vcov(m$fit)), tolerance = 1e-12)
})

test_that("the numerical Hessian remains an interior fallback", {
  m <- fit_reference()
  reference <- as.matrix(stats::vcov(m$fit))
  m$fit@optinfo$derivs$Hessian <- NULL
  V <- .full_vcov(m$fit, m$devfun)
  fixed <- names(lme4::fixef(m$fit))

  expect_identical(attr(V, "source"), "numerical")
  expect_equal(V[fixed, fixed], reference, tolerance = 1e-5)
})

test_that("subject errors include the fixed and random covariance", {
  predict_method <- utils::getS3method("predict", "merMod")
  skip_if_not("se.fit" %in% names(formals(predict_method)))

  m <- fit_reference()
  ag <- m$data$agg
  form <- cbind(y, n - y) ~ 0 + d_D + d_I +
    (0 + d_D + d_I | subj)
  fit <- lme4::glmer(
    form, data = ag, family = stats::binomial("probit"), nAGQ = 1L,
    control = lme4::glmerControl(optimizer = "bobyqa")
  )
  devfun <- lme4::glmer(
    form, data = ag, family = stats::binomial("probit"), nAGQ = 1L,
    control = lme4::glmerControl(optimizer = "bobyqa"), devFunOnly = TRUE
  )

  phi <- c(lme4::getME(fit, "theta"), lme4::fixef(fit))
  V <- matrix(0, length(phi), length(phi),
              dimnames = list(names(phi), names(phi)))
  fixed <- names(lme4::fixef(fit))
  V[fixed, fixed] <- as.matrix(stats::vcov(fit))
  ours <- .conditional_se(fit, V, devfun)

  subjects <- levels(lme4::getME(fit, "flist")$subj)
  newdata <- data.frame(subj = factor(subjects, levels = subjects),
                        d_D = 1, d_I = 0)
  reference <- suppressWarnings(
    stats::predict(fit, newdata = newdata, se.fit = TRUE)$se.fit
  )
  direct <- ours[ours$term == "d_D", ]
  direct <- direct$se[match(subjects, as.character(direct$grp))]

  expect_equal(unname(direct), unname(reference), tolerance = 1e-3)
})

test_that("subject errors follow the full two-block mixed system", {
  m <- fit_reference()
  fit <- m$fit

  phi <- c(lme4::getME(fit, "theta"), lme4::fixef(fit))
  V <- matrix(0, length(phi), length(phi),
              dimnames = list(names(phi), names(phi)))
  fixed <- names(lme4::fixef(fit))
  V[fixed, fixed] <- as.matrix(stats::vcov(fit))
  ours <- .conditional_se(fit, V, m$devfun)

  # The joint system gives an independent fixed-theta reference.
  Z <- t(as.matrix(lme4::getME(fit, "Zt")))
  X <- as.matrix(lme4::getME(fit, "X"))
  Lambda <- t(as.matrix(lme4::getME(fit, "Lambdat")))
  weight <- fit@resp$sqrtWrkWt()^2
  ZL <- Z %*% Lambda
  A <- crossprod(ZL, weight * ZL) + diag(ncol(ZL))
  B <- crossprod(ZL, weight * X)
  D <- crossprod(X, weight * X)
  joint <- solve(rbind(cbind(A, B), cbind(t(B), D)))

  subjects <- levels(lme4::getME(fit, "flist")$subj)
  n_subjects <- length(subjects)
  sensitivity_start <- 2L * n_subjects
  reference <- unlist(lapply(seq_along(subjects), function(j) {
    v_direct <- c(Lambda[sensitivity_start + 2L * j - 1L, ],
                  as.numeric(fixed == "d_D"))
    v_indirect <- c(Lambda[sensitivity_start + 2L * j, ],
                    as.numeric(fixed == "d_I"))
    c(d_D = sqrt(drop(crossprod(v_direct, joint %*% v_direct))),
      d_I = sqrt(drop(crossprod(v_indirect, joint %*% v_indirect))))
  }))
  names(reference) <- NULL

  actual <- unlist(lapply(seq_along(subjects), function(j) {
    rows <- ours$grp == subjects[j] & ours$term %in% c("d_D", "d_I")
    values <- ours[rows, ]
    values$se[match(c("d_D", "d_I"), values$term)]
  }))

  expect_equal(actual, reference, tolerance = 1e-3)
})

test_that("the variance components match VarCorr", {
  m  <- fit_reference()
  vc <- lme4::VarCorr(m$fit)
  b  <- vc[[which(vapply(vc, function(z) all(c("d_D", "d_I") %in% colnames(z)),
                         TRUE))]]
  expect_equal(unname(m$pars$est[["s2_D"]]), b["d_D", "d_D"], tolerance = 1e-9)
  expect_equal(unname(m$pars$est[["s2_I"]]), b["d_I", "d_I"], tolerance = 1e-9)
  expect_equal(unname(m$pars$est[["s_DI"]]), b["d_D", "d_I"], tolerance = 1e-9)
  expect_equal(m$pars$rho, attr(b, "correlation")["d_D", "d_I"], tolerance = 1e-7)

  criteria <- vc[[which(vapply(vc, function(z) all(c("c_D", "c_I") %in%
                                                    colnames(z)), TRUE))]]
  expect_equal(m$pars$criterion_cor$estimate,
               attr(criteria, "correlation")["c_D", "c_I"], tolerance = 1e-7)
  expect_identical(m$design$n_fixed, 4L)
  expect_equal(m$design$n_var, 6)
})

test_that("the printed model uses the classical criterion sign", {
  m <- fit_reference()
  output <- capture.output(summary(m))
  line <- output[grepl("^  criterion +Direct", output)]
  estimate <- sprintf("%.4f", -unname(lme4::fixef(m$fit)[["c_D"]]))
  z <- -unname(lme4::fixef(m$fit)[["c_D"]]) /
    sqrt(diag(as.matrix(stats::vcov(m$fit))))[["c_D"]]

  expect_match(line, estimate, fixed = TRUE)
  expect_match(line, sprintf("%.2f", z), fixed = TRUE)
  expect_output(summary(m), "cor\\(c\\)")
  expect_output(summary(m), paste0(.usdt_chars()$delta, "d' \\(D - I\\)"))
})

test_that("the model plots use the expected subject sensitivities", {
  m <- fit_reference()
  figure <- plot(m, type = "shrinkage")
  values <- figure$data
  moments <- sdt_moments(m$data)
  random <- lme4::ranef(m$fit, condVar = FALSE)$subj
  fixed <- lme4::fixef(m$fit)

  direct <- moments[moments$task == m$data$meta$labels[["direct"]], ]
  indirect <- moments[moments$task == m$data$meta$labels[["indirect"]], ]
  direct_row <- match(values$subject, as.character(direct$subj))
  indirect_row <- match(values$subject, as.character(indirect$subj))
  random_row <- match(values$subject, rownames(random))

  expect_s3_class(figure, "ggplot")
  expect_equal(values$observed_direct, direct$dprime[direct_row])
  expect_equal(values$observed_indirect, indirect$dprime[indirect_row])
  expect_equal(values$model_direct,
               fixed[["d_D"]] + random$d_D[random_row])
  expect_equal(values$model_indirect,
               fixed[["d_I"]] + random$d_I[random_row])
  expect_no_error(ggplot2::ggplot_build(figure))

  caterpillar <- plot(m, type = "caterpillar")
  intervals <- caterpillar$data
  observed_intervals <- intervals[
    intervals$method == "Observed" & intervals$task == "Direct", ]
  model_intervals <- intervals[
    intervals$method == "Model-estimated" & intervals$task == "Direct", ]
  moments <- sdt_moments(m$data, variances = TRUE)
  moments <- moments[moments$task == "Direct", ]
  moment_row <- match(observed_intervals$subject, as.character(moments$subj))
  conditional <- as.data.frame(lme4::ranef(m$fit, condVar = TRUE))
  conditional <- conditional[conditional$term == "d_D", ]
  conditional_row <- match(model_intervals$subject,
                           as.character(conditional$grp))
  critical <- stats::qnorm(1 - (1 - m$level) / 2)

  expect_s3_class(caterpillar, "ggplot")
  expect_equal(observed_intervals$estimate, moments$dprime[moment_row])
  expect_equal(observed_intervals$se, moments$se_gg[moment_row])

  # `observed_se` selects the other standard error without touching the rest.
  miller <- plot(m, type = "caterpillar", observed_se = "miller")$data
  miller <- miller[miller$method == "Observed" & miller$task == "Direct", ]
  expect_equal(miller$se, moments$se_miller[moment_row])
  expect_equal(miller$estimate, observed_intervals$estimate)
  expect_error(plot(m, type = "shrinkage", observed_se = "gg"),
               "only used by")
  expect_equal(model_intervals$estimate,
               fixed[["d_D"]] + conditional$condval[conditional_row])

  # The model interval includes all fitted uncertainty.
  condsd <- conditional$condsd[conditional_row]
  expect_true(all(model_intervals$se > condsd))

  joint <- .conditional_se(m$fit, m$pars$V_full, m$devfun)
  joint <- joint[joint$term == "d_D", ]
  expect_equal(model_intervals$se,
               joint$se[match(model_intervals$subject, as.character(joint$grp))])

  # A fixed-only covariance reproduces the conditional-theta calculation.
  fixed_vcov <- matrix(0, nrow(m$pars$V_full), ncol(m$pars$V_full),
                       dimnames = dimnames(m$pars$V_full))
  fixed_names <- names(fixed)
  fixed_vcov[fixed_names, fixed_names] <- as.matrix(stats::vcov(m$fit))
  fixed_only <- .conditional_se(m$fit, fixed_vcov, m$devfun)
  fixed_only <- fixed_only[fixed_only$term == "d_D", ]
  fixed_only <- fixed_only$se[match(model_intervals$subject,
                                    as.character(fixed_only$grp))]
  expect_true(all(fixed_only < model_intervals$se))
  expect_true(all(fixed_only > condsd))
  expect_equal(intervals$conf.low,
               intervals$estimate - critical * intervals$se)
  expect_equal(intervals$conf.high,
               intervals$estimate + critical * intervals$se)

  # Every panel reports the descriptive share of its displayed intervals.
  zero_summary <- .caterpillar_zero_summary(intervals)
  expected <- vapply(levels(intervals$panel), function(panel) {
    100 * mean(intervals$status[intervals$panel == panel] == "Includes zero")
  }, numeric(1L))
  expect_identical(nrow(zero_summary), 4L)
  expect_equal(zero_summary$percentage, unname(expected))
  expect_true(all(grepl("Intervals including 0:", zero_summary$label,
                         fixed = TRUE)))
  expect_no_error(ggplot2::ggplot_build(caterpillar))
  expect_error(plot(m, type = "unknown"), "caterpillar")
})

test_that("the regression plot separates observed and model relationships", {
  m <- fit_reference()
  figure <- plot(m, type = "regression")
  values <- .regression_data(m)
  naive <- stats::lm(observed_indirect ~ observed_direct, data = values$points)
  latent <- m$tests[m$tests$hypothesis == "H3", ]

  expect_s3_class(figure, "ggplot")
  expect_identical(levels(figure$data$panel),
                   c("Observed estimates", "Model-estimated sensitivities"))
  expect_equal(values$naive_intercept, unname(stats::coef(naive)[1L]))
  expect_equal(values$naive_slope, unname(stats::coef(naive)[2L]))
  expect_equal(unname(values$naive_p),
               unname(summary(naive)$coefficients[, "Pr(>|t|)"]))
  expect_equal(values$intercept,
               latent$estimate[latent$term == "intercept"])
  expect_equal(values$slope, latent$estimate[latent$term == "slope"])
  expect_true(all(vapply(c("italic(b)[0]", "italic(b)[1]", "italic(p)"), function(x) {
    all(grepl(x, values$annotations$label, fixed = TRUE))
  }, logical(1L))))
  expect_no_error(lapply(values$annotations$label, function(x) parse(text = x)))
  expect_no_error(ggplot2::ggplot_build(figure))
})

test_that("ROC plots transform the fitted sensitivities and criteria", {
  m <- fit_reference()
  fixed <- lme4::fixef(m$fit)
  population <- plot(m, type = "roc")
  values <- population$data
  direct <- values[values$task == "Direct", ]
  midpoint <- which(direct$false_alarm == 0.5)

  expect_s3_class(population, "ggplot")
  expect_equal(direct$hit[midpoint], stats::pnorm(fixed[["d_D"]]))
  expect_equal(unique(direct$auc), stats::pnorm(fixed[["d_D"]] / sqrt(2)))
  expect_identical(unique(values$interval), "Wald")
  expect_true(all(values$conf.low <= values$hit))
  expect_true(all(values$conf.high >= values$hit))

  operating <- .roc_data(m, band = FALSE)$points
  direct_point <- operating[operating$task == "Direct", ]
  expect_equal(direct_point$false_alarm,
               stats::pnorm(fixed[["c_D"]] - fixed[["d_D"]] / 2))
  expect_equal(direct_point$hit,
               stats::pnorm(fixed[["c_D"]] + fixed[["d_D"]] / 2))

  id <- as.character(m$data$agg$subj[1L])
  subject <- plot(m, type = "roc", subject_id = id)
  random <- lme4::ranef(m$fit, condVar = FALSE)$subj
  expected <- fixed[["d_D"]] + random[id, "d_D"]
  subject_direct <- subject$data[subject$data$task == "Direct", ]
  midpoint <- which(subject_direct$false_alarm == 0.5)

  expect_equal(subject_direct$hit[midpoint], stats::pnorm(expected))
  expect_identical(unique(subject$data$interval), "conditional")
  expect_no_error(ggplot2::ggplot_build(population))
  expect_no_error(ggplot2::ggplot_build(subject))
  expect_error(plot(m, type = "roc", subject_id = "not-a-subject"),
               "was not found")
})

test_that("the bootstrap separates population and subject parameters", {
  m <- fit_reference(seed = 8L)
  statistic <- .boot_statistic(m$fit, "d_D", "d_I")
  fitted <- statistic(m$fit)
  samples <- rbind(fitted, fitted)
  effects <- .boot_effects(fitted, samples, attr(statistic, "layout"))
  fixed <- lme4::fixef(m$fit)
  random <- lme4::ranef(m$fit, condVar = FALSE)$subj

  expect_identical(names(effects$population$t0),
                   c("c_D", "c_I", "d_D", "d_I"))
  expect_identical(names(effects$variance$t0), c("s2_D", "s2_I"))
  expect_identical(dim(effects$variance$t), c(2L, 2L))
  expect_identical(dim(effects$population$t), c(2L, 4L))
  expect_identical(dim(effects$subjects$t),
                   c(2L, nrow(random), 4L))
  expect_equal(effects$population$t0[["d_D"]], fixed[["d_D"]])
  expect_equal(unname(effects$subjects$t0[, "d_D"]),
               unname(fixed[["d_D"]] + random$d_D))
})

test_that("the block is found whichever order the random terms are written in", {
  m  <- fit_reference()
  ag <- m$data$agg
  f2 <- lme4::glmer(cbind(y, n - y) ~ 0 + c_D + c_I + d_D + d_I +
                      (0 + d_D + d_I | subj) + (0 + c_D + c_I | subj),
                    data = ag, family = stats::binomial("probit"),
                    control = lme4::glmerControl(optimizer = "bobyqa"))
  f3 <- lme4::glmer(cbind(y, n - y) ~ 0 + c_D + c_I + d_D + d_I +
                      (0 + c_D + c_I | subj) + (0 + d_D + d_I | subj),
                    data = ag, family = stats::binomial("probit"),
                    control = lme4::glmerControl(optimizer = "bobyqa"))
  p2 <- .usdt_pars(f2, "d_D", "d_I")
  p3 <- .usdt_pars(f3, "d_D", "d_I")
  expect_equal(unname(p2$est), unname(p3$est), tolerance = 1e-5)
})

test_that("the hypotheses are internally consistent", {
  m <- fit_reference()
  e <- m$pars$est
  t <- m$tests

  # The slope has two equivalent expressions.
  slope <- t$estimate[t$term == "slope"]
  expect_equal(slope, e[["s_DI"]] / e[["s2_D"]], tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(slope, m$pars$rho * sqrt(e[["s2_I"]] / e[["s2_D"]]),
               tolerance = 1e-8, ignore_attr = TRUE)

  # The intercept adjusts the indirect mean by the slope.
  expect_equal(t$estimate[t$term == "intercept"],
               e[["gamma_I"]] - slope * e[["gamma_D"]], tolerance = 1e-10,
               ignore_attr = TRUE)

  # The difference compares the two fixed effects.
  expect_equal(t$estimate[t$hypothesis == "H1"],
               e[["gamma_D"]] - e[["gamma_I"]], tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("the correlation and the slope share one null hypothesis", {
  m <- fit_reference()
  t <- m$tests
  expect_equal(t$p.value[t$hypothesis == "H2"],
               t$p.value[t$term == "slope"], tolerance = 1e-12)
  expect_equal(t$statistic[t$hypothesis == "H2"],
               t$statistic[t$term == "slope"], tolerance = 1e-12)
})

test_that("each hypothesis uses its fixed interval method", {
  m <- fit_reference()
  expect_identical(m$tests$ci_method[m$tests$hypothesis == "H1"], "Wald")
  expect_identical(m$tests$ci_method[m$tests$hypothesis == "H2"], "Fisher-z")
  expect_identical(m$tests$ci_method[m$tests$term == "intercept"], "delta")
  expect_identical(m$tests$ci_method[m$tests$term == "slope"], "Wald")

  # Fisher-z keeps the correlation inside its valid range.
  expect_lte(m$tests$conf.high[m$tests$hypothesis == "H2"], 1)
  expect_gte(m$tests$conf.low[m$tests$hypothesis == "H2"], -1)
})

test_that("the layer-2 functions agree with the fitted table", {
  m <- fit_reference()

  # The functions accept the fitted model.
  expect_equal(latent_regression(m$fit)$estimate,
               m$tests$estimate[m$tests$hypothesis == "H3"], tolerance = 1e-8)
  expect_equal(latent_cor(m$fit)$estimate, m$pars$rho, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(sensitivity_diff(m$fit)$estimate,
               m$tests$estimate[m$tests$hypothesis == "H1"], tolerance = 1e-8)

  # The functions also accept the full uSDT result.
  expect_equal(latent_regression(m)$estimate, latent_regression(m$fit)$estimate,
               tolerance = 1e-8)
  expect_identical(nrow(usdt_tests(m)), 4L)
})

test_that("a boundary fit returns stable hypothesis rows", {
  # Very little direct variation makes the slope hard to identify.
  set.seed(12)
  df <- usdt_simulate(n_subj = 40L, n_trials = 100L, gamma_D = 0.6,
                      gamma_I = 0.2, sd_D = 0.02, sd_I = 0.05, rho = 0.3)
  d  <- usdt_data_long(df, task_col = "task",
                       task_levels      = c(direct = "D", indirect = "I"),
                       subject_col      = "subj",
                       condition_col    = "cond",
                       condition_levels = c(signal = 1, noise = 0),
                       response_col     = "response",
                       response_levels  = c(signal = 1, noise = 0))
  m <- suppressWarnings(hsdt(d))
  s <- m$tests[m$tests$term == "slope", ]

  expect_identical(nrow(m$tests), 4L)
  expect_identical(s$status, "not estimable")
  expect_true(is.na(s$conf.low))
  expect_true(nzchar(s$reason))
  expect_output(summary(m), "Not estimable")
  expect_error(plot(m, type = "caterpillar"), "unavailable")
})

test_that("extra glmer arguments reach the fitted model", {
  m0 <- fit_reference()
  weights <- rep(2, nrow(m0$data$agg))
  m <- fit_reference(weights = weights)
  expect_equal(stats::model.weights(stats::model.frame(m$fit)), weights)

  agq1 <- fit_reference(nAGQ = 1)
  expect_identical(agq1$fit@call$nAGQ, 1L)
  expect_error(fit_reference(nAGQ = 0), "must be 1")
})

test_that("the model structure and interval methods cannot be changed", {
  expect_error(fit_reference(re = "diagonal"), "not a uSDT option")
  expect_error(fit_reference(ci = "wald"), "not a uSDT option")
})

test_that("failure with every optimizer stops the analysis", {
  expect_error(fit_reference(not_a_glmer_argument = TRUE),
               "did not converge with any available optimizer")
})

test_that("hsdt rejects anything that is not a usdt_data object", {
  expect_error(hsdt(data.frame(x = 1)), "usdt_data_long")
})
