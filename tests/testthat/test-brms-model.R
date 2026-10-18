designs <- expand.grid(parameterization = c("correlation", "regression"),
                       unequal_variances = c(FALSE, TRUE),
                       indirect_free = c(TRUE, FALSE),
                       stringsAsFactors = FALSE)
design_of <- function(i) {
  free_c <- c(1L, as.integer(designs$indirect_free[i]))
  list(criteria = c("c_D", "c_I")[free_c == 1L], free_c = free_c,
       unequal_variances = designs$unequal_variances[i],
       parameterization = designs$parameterization[i])
}

test_that("the brms formula gives every row its SDT predictor", {
  agg <- expand.grid(subj = c("s1", "s2"), sig = c(TRUE, FALSE),
                     task = factor(c("D", "I")), stringsAsFactors = FALSE)
  agg$y <- 3L
  agg$n <- 10L
  data <- .brms_data(agg)
  v <- list(dD = 0.8, dI = 0.3, cD = 0.2, cI = -0.1, sD = 1.4, sI = 0.7)
  direct <- data$D == 1

  for (i in seq_len(nrow(designs))) {
    design <- design_of(i)
    eta <- eval(.brms_formula(design)$formula[[3L]], c(as.list(data), v))

    # A fixed criterion makes HR + FAR = 1, which is 0 under equal variances.
    d <- ifelse(direct, v$dD, v$dI)
    s <- if (design$unequal_variances) ifelse(direct, v$sD, v$sI) else 1
    criterion <- ifelse(direct, v$cD, if (design$free_c[2L] == 1L) v$cI else
      d * (1 - s) / (2 * (1 + s)))
    expect_equal(eta, ifelse(data$sig == 1, (d / 2 - criterion) / s,
                             -d / 2 - criterion), label = paste(designs[i, ]))
  }
})

test_that("each prior a design uses becomes one brms prior of its own", {
  set.seed(6)
  d <- suppressWarnings(usdt_data_long(
    make_trials(n_subj = 8, n_trials = 20), task_col = "task",
    task_levels = c(direct = "D", indirect = "I"), subject_col = "subj",
    condition_col = "cond", condition_levels = c(signal = 1, noise = 0),
    response_col = "response", response_levels = c(signal = 1, noise = 0)))
  priors <- usdt_priors(dprime = list(direct = "normal(0, 1)",
                                      indirect = "skew_normal(0, 1, 2)"),
                        sd_dprime = "gamma(2, 4)",
                        sd_signal = "student_t(4, 1, 0.5)",
                        intercept = "std_normal()")

  for (i in seq_len(nrow(designs))) {
    design <- design_of(i)
    regression <- design$parameterization == "regression"
    model <- .brms_model(d, design, priors)
    table <- brms::validate_prior(model$prior, model$formula,
                                  data = model$data,
                                  family = stats::binomial("probit"))
    used <- .used_priors(priors, design)
    expect_identical(any(used$parameter == "intercept"), regression)
    stated <- used[used$parameter != "intercept", ]
    user <- table[table$source == "user" & nzchar(table$class), ]
    label <- paste(designs[i, ])
    expect_identical(nrow(user), nrow(stated), label = label)
    # brms keeps an LKJ prior on the Cholesky factor of the correlations.
    expect_setequal(user$prior, sub("^lkj\\(", "lkj_corr_cholesky(",
                                    stated$prior))
    # Every mean has a prior but the indirect one under the latent
    # regression, whose prior is that of the intercept.
    flat <- table$class == "b" & !nzchar(table$coef) & nzchar(table$nlpar) &
      table$prior %in% c("", "(flat)")
    expect_identical(table$nlpar[flat], if (regression) "dI" else character(0),
                     label = label)

    code <- suppressMessages(brms::stancode(
      model$formula, data = model$data, family = stats::binomial("probit"),
      prior = model$prior))
    lpdf <- paste0(sub("^lkj$", "lkj_corr_cholesky", unique(stated$family)),
                   "_lpdf")
    for (f in lpdf) expect_match(code, f, fixed = TRUE, label = label)

    # The intercept prior reads the direct and indirect d' of the first group
    # of brms, in that order.
    expect_match(code, "r_1_dD_1 = r_1[, 1];", fixed = TRUE, label = label)
    expect_match(code, "r_1_dI_2 = r_1[, 2];", fixed = TRUE, label = label)
    expect_identical(grepl(paste0(
      "target += std_normal_lpdf(b_dI[1] - L_1[2, 1] * sd_1[2] / sd_1[1] * ",
      "b_dD[1]);"), code, fixed = TRUE), regression, label = label)
  }
})

test_that("the intercept prior passes the arguments of its family", {
  row <- .prior_rows(usdt_priors(intercept = "student_t(5, 0.2, 0.4)"),
                     "intercept", "both")
  expect_identical(.intercept_prior(row)$prior, paste0(
    "target += student_t_lpdf(b_dI[1] - L_1[2, 1] * sd_1[2] / sd_1[1] * ",
    "b_dD[1] | 5, 0.2, 0.4)"))
})
