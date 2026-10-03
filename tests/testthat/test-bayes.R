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
})
