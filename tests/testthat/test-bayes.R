bayes_data <- function(df) {
  usdt_data_long(df, task_col = "task",
                 task_levels      = c(direct = "D", indirect = "I"),
                 subject_col      = "subj",
                 condition_col    = "cond",
                 condition_levels = c(signal = 1, noise = 0),
                 response_col     = "response",
                 response_levels  = c(signal = 1, noise = 0))
}

test_that("the Stan data keep every count, with zero trials where a task is missing", {
  set.seed(1)
  df <- make_trials(n_subj = 6, n_trials = 20)
  df <- df[!(df$subj == 6 & df$task == "I"), ]
  d <- suppressWarnings(bayes_data(df))
  stan <- .stan_data(d, c(1L, 1L), FALSE, usdt_priors())$data

  expect_identical(dim(stan$hit), c(6L, 2L))
  expect_identical(sum(stan$hit) + sum(stan$fa), sum(d$agg$y))
  expect_identical(sum(stan$N_hit) + sum(stan$N_fa), sum(d$agg$n))
  expect_identical(c(stan$N_hit[6, 2], stan$N_fa[6, 2]), c(0L, 0L))
})

test_that("the Bayesian fit reproduces the frequentist group-level estimates", {
  # Compiling and sampling Stan takes minutes, so this check runs on request.
  skip_if_not(identical(Sys.getenv("USDT_TEST_STAN"), "true"),
              "set USDT_TEST_STAN=true to fit the Stan model")
  skip_if_not_installed("rstan")
  old <- options(uSDT.cache_dir = file.path(tempdir(), "usdt-stan"))
  on.exit(options(old), add = TRUE)

  set.seed(2)
  d <- bayes_data(make_trials(n_subj = 40, n_trials = 100))
  m <- hsdt(d, estimation = "bayesian", chains = 2, iter = 1500,
            warmup = 500, seed = 1)
  f <- hsdt(d)

  expect_identical(m$tests$term, f$tests$term)
  expect_equal(m$tests$estimate[1L], f$tests$estimate[1L], tolerance = 0.05)
  expect_true(all(m$tests$rhat < 1.05))
  expect_identical(usdt_tests(m), m$tests)

  # The latent line meets the H3 table at zero, and every plot draws.
  origin <- .regression_data(m)$origin
  expect_equal(origin$fit, m$tests$estimate[3L])
  for (type in c("regression", "shrinkage", "caterpillar", "roc")) {
    expect_s3_class(plot(m, type = type), "ggplot")
  }

  # Reliability lies inside its credible interval, close to the frequentist one.
  r <- usdt_reliability(m)
  expect_true(all(r$tasks$conf.low <= r$tasks$reliability &
                    r$tasks$reliability <= r$tasks$conf.high))
  expect_equal(r$tasks$reliability, usdt_reliability(f)$tasks$reliability,
               tolerance = 0.2)

  # The default Bayes factors test the three point nulls and draw.
  skip_if_not_installed("logspline")
  b <- usdt_bf(m, plot = FALSE)
  expect_identical(b$hypothesis,
                   c(paste0(.usdt_chars()$delta, "d' = 0"), "rho = 0",
                     "intercept = 0"))
  expect_equal(b$BF10, 1 / b$BF01)
  expect_s3_class(plot(b), "ggplot")
})

test_that("brms variables are found whatever the order of the correlation", {
  variables <- c("b_direct", "b_taskD:cond", "b_taskI:cond", "sd_subj__direct",
                 "sd_subj__taskD:cond", "sd_subj__taskI:cond",
                 "cor_subj__taskI:cond__taskD:cond", "lp__")
  expect_identical(
    .brms_variables(variables, "taskD:cond", "taskI:cond"),
    c("b_taskD:cond", "b_taskI:cond", "sd_subj__taskD:cond",
      "sd_subj__taskI:cond", "cor_subj__taskI:cond__taskD:cond"))
  expect_error(.brms_variables(variables, "d_D", "d_I"), "not found")
  expect_error(.brms_variables(variables, "taskD:cond", "direct"),
               "no group-level term")
})

test_that("aggregated counts keep every trial and response", {
  counts <- usdt_aggregate(vadillo_awareness, response = "judged.old",
                           by = c("subj", "condition", "set.size"))
  expect_identical(sum(counts$n), nrow(vadillo_awareness))
  expect_identical(sum(counts$y), as.integer(sum(vadillo_awareness$judged.old)))
  expect_error(usdt_aggregate(vadillo_awareness, "condition", "subj"), "binary")
})

test_that("a brms model is tested from its posterior", {
  skip_if_not(identical(Sys.getenv("USDT_TEST_STAN"), "true"),
              "set USDT_TEST_STAN=true to fit the brms model")
  skip_if_not_installed("brms")

  set.seed(3)
  df <- make_trials(n_subj = 30, n_trials = 60)
  df$cond_D <- (df$cond - 0.5) * (df$task == "D")
  df$cond_I <- (df$cond - 0.5) * (df$task == "I")
  df$direct <- as.integer(df$task == "D")
  df$indirect <- as.integer(df$task == "I")
  counts <- usdt_aggregate(df, "response",
                           c("subj", "direct", "indirect", "cond_D", "cond_I"))
  fit <- suppressWarnings(brms::brm(
    y | trials(n) ~ 0 + direct + indirect + cond_D + cond_I +
      (0 + direct + indirect | subj) + (0 + cond_D + cond_I | subj),
    data = counts, family = binomial("probit"),
    prior = brms::prior(normal(0, 1), class = b), chains = 2, iter = 1000,
    seed = 1, refresh = 0, silent = 2))

  tests <- usdt_tests(fit, direct = "cond_D", indirect = "cond_I")
  expect_identical(tests$hypothesis, c("H1", "H2", "H3", "H3"))
  expect_identical(latent_cor(fit, "cond_D", "cond_I")$p.value,
                   tests$p.value[2L])
})
