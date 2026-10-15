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
  skip_if_not_installed("brms")
  agg <- expand.grid(subj = c("s1", "s2"), sig = c(TRUE, FALSE),
                     task = factor(c("D", "I")), stringsAsFactors = FALSE)
  agg$y <- 3L
  agg$n <- 10L
  data <- .brms_data(agg)
  v <- list(dD = 0.8, dI = 0.3, cD = 0.2, cI = -0.1, sD = 1.4, sI = 0.7,
            b0 = 0.1, b1 = 0.5, e = 0.05)
  direct <- data$D == 1

  for (i in seq_len(nrow(designs))) {
    design <- design_of(i)
    eta <- eval(.brms_formula(design)$formula[[3L]], c(as.list(data), v))

    # A fixed criterion makes HR + FAR = 1, which is 0 under equal variances.
    d <- ifelse(direct, v$dD, if (design$parameterization == "regression")
      v$b0 + v$b1 * v$dD + v$e else v$dI)
    s <- if (design$unequal_variances) ifelse(direct, v$sD, v$sI) else 1
    criterion <- ifelse(direct, v$cD, if (design$free_c[2L] == 1L) v$cI else
      d * (1 - s) / (2 * (1 + s)))
    expect_equal(eta, ifelse(data$sig == 1, (d / 2 - criterion) / s,
                             -d / 2 - criterion), label = paste(designs[i, ]))
  }
})

test_that("each prior a design uses becomes one brms prior of its own", {
  skip_if_not_installed("brms")
  set.seed(6)
  d <- suppressWarnings(usdt_data_long(
    make_trials(n_subj = 8, n_trials = 20), task_col = "task",
    task_levels = c(direct = "D", indirect = "I"), subject_col = "subj",
    condition_col = "cond", condition_levels = c(signal = 1, noise = 0),
    response_col = "response", response_levels = c(signal = 1, noise = 0)))
  priors <- usdt_priors(dprime = list(direct = "normal(0, 1)",
                                      indirect = "skew_normal(0, 1, 2)"),
                        sd_dprime = "gamma(2, 4)",
                        sd_signal = "student_t(4, 1, 0.5)")

  for (i in seq_len(nrow(designs))) {
    design <- design_of(i)
    model <- .brms_model(d, design, priors)
    table <- brms::validate_prior(model$prior, model$formula,
                                  data = model$data,
                                  family = stats::binomial("probit"))
    used <- .used_priors(priors, design)
    user <- table[table$source == "user", ]
    label <- paste(designs[i, ])
    expect_identical(nrow(user), nrow(used), label = label)
    # brms keeps an LKJ prior on the Cholesky factor of the correlations.
    expect_setequal(user$prior, sub("^lkj\\(", "lkj_corr_cholesky(",
                                    used$prior))
    expect_false(any(table$class == "b" & table$prior == "(flat)" &
                       nzchar(table$nlpar)), label = label)

    code <- suppressMessages(brms::stancode(
      model$formula, data = model$data, family = stats::binomial("probit"),
      prior = model$prior))
    lpdf <- paste0(sub("^lkj$", "lkj_corr_cholesky", unique(used$family)),
                   "_lpdf")
    for (f in lpdf) expect_match(code, f, fixed = TRUE, label = label)
  }
})
