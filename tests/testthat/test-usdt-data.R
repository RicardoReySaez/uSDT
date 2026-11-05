# test-usdt-data.R
# This script tests the preparation of model data.
# Author: Ricardo Rey-Sáez
# Last modified: 07-09-2026

# This function creates a small data set for both input routes.
make_binary <- function(seed = 1L, n_subj = 30L, n_trials = 60L) {
  set.seed(seed)
  make_trials(n_subj = n_subj, n_trials = n_trials)
}

# This function keeps the common input options in one place.
long_args <- function(df, ...) {
  args <- list(data = df, task_col = "task",
               task_levels      = c(awareness = "A", influence = "I"),
               subject_col      = "subj",
               condition_col    = "cond",
               condition_levels = c(signal = 1, noise = 0),
               response_col     = "response",
               response_levels  = c(signal = 1, noise = 0))
  do.call(usdt_data_long, utils::modifyList(args, list(...)))
}

test_that("both entry points give the same object", {
  df <- make_binary()
  d1 <- long_args(df)
  d2 <- usdt_data_tasks(awareness = df[df$task == "A", ],
                        influence = df[df$task == "I", ],
                        subject_col      = "subj",
                        condition_col    = "cond",
                        condition_levels = c(signal = 1, noise = 0),
                        response_col     = "response",
                        response_levels  = c(signal = 1, noise = 0))
  keep <- c("subj", "task", "sig", "y", "n", "cond", "c_A", "c_I", "d_A", "d_I")
  expect_identical(d1$agg[, keep], d2$agg[, keep])
})

test_that("counts, proportions and SDT tables reproduce trial-level input", {
  df  <- make_binary()
  ref <- long_args(df)

  # This input contains cell counts.
  cnt <- ref$agg
  cnt$tk <- as.character(cnt$task)
  d_cnt <- usdt_data_long(cnt, task_col = "tk",
                          task_levels      = c(awareness = "A",
                                               influence = "I"),
                          subject_col      = "subj",
                          condition_col    = "sig",
                          condition_levels = c(signal = TRUE, noise = FALSE),
                          successes_col    = "y", trials_col = "n")
  expect_identical(d_cnt$agg$y, ref$agg$y)
  expect_identical(d_cnt$agg$n, ref$agg$n)

  # This input contains response proportions.
  prp <- cnt; prp$p <- prp$y / prp$n
  d_prp <- usdt_data_long(prp, task_col = "tk",
                          task_levels      = c(awareness = "A",
                                               influence = "I"),
                          subject_col      = "subj",
                          condition_col    = "sig",
                          condition_levels = c(signal = TRUE, noise = FALSE),
                          successes_col    = "p", trials_col = "n",
                          successes_type   = "proportions")
  expect_identical(d_prp$agg$y, ref$agg$y)

  # This input contains the four signal detection counts.
  w <- stats::reshape(ref$agg[, c("subj", "task", "sig", "y", "n")],
                      idvar = c("subj", "task"), timevar = "sig",
                      direction = "wide")
  names(w) <- c("subj", "task", "fa", "n_noise", "hit", "n_sig")
  w$miss <- w$n_sig - w$hit
  w$cr   <- w$n_noise - w$fa
  w$task <- as.character(w$task)
  d_sdt <- usdt_data_long(w, task_col = "task",
                          task_levels = c(awareness = "A", influence = "I"),
                          subject_col = "subj",
                          sdt_cols    = c(hit = "hit", miss = "miss",
                                          fa = "fa", cr = "cr"))
  expect_identical(d_sdt$agg$y, ref$agg$y)
  expect_identical(d_sdt$agg$n, ref$agg$n)
})

test_that("the two codings differ only where they should", {
  df  <- make_binary()
  dev <- long_args(df, coding = "deviation")
  trt <- long_args(df, coding = "treatment")

  expect_setequal(unique(dev$agg$cond), c(-0.5, 0.5))
  expect_setequal(unique(trt$agg$cond), c(0, 1))

  # The coding does not change the cell counts.
  expect_identical(dev$agg$y, trt$agg$y)
})

test_that("a Meyen-split task has its criterion fixed under deviation coding only", {
  set.seed(5)
  df <- make_trials(n_subj = 40L, n_trials = 100L, rt = TRUE)
  # The median split replaces this temporary response.
  df$response[df$task == "I"] <- 0L

  mk <- function(coding) {
    usdt_data_tasks(
      awareness = df[df$task == "A", ],
      influence = df[df$task == "I", ],
      subject_col      = "subj",
      condition_col    = "cond",
      condition_levels = c(signal = 1, noise = 0),
      response_col     = c(awareness = "response", influence = "rt"),
      response_levels  = list(awareness = c(signal = 1, noise = 0),
                              influence = c(signal = "faster",
                                            noise = "slower")),
      dichotomize = "influence", coding = coding)
  }

  dev <- mk("deviation")
  trt <- suppressMessages(mk("treatment"))

  # Deviation coding makes the criterion equal to zero in this design.
  expect_true(dev$meta$criterion_zero[["influence"]])
  expect_lt(dev$meta$criterion[["influence"]]$mean_abs, 0.02)

  # Treatment coding makes the criterion equal to d'/2, so it is estimated.
  expect_false(trt$meta$criterion_zero[["influence"]])

  # The awareness task keeps its criterion with both coding choices.
  expect_false(dev$meta$criterion_zero[["awareness"]])
})

test_that("the formula follows the criterion decision", {
  set.seed(6)
  df <- make_trials(n_subj = 30L, n_trials = 60L)
  d  <- long_args(df)

  # Both criteria are estimated when no task was dichotomized.
  f_all <- .usdt_formula(d)
  expect_setequal(f_all$criteria, c("c_A", "c_I"))
  expect_identical(f_all$n_fixed, 4L)

  # Fixing the influence criterion at zero removes one term.
  d$meta$criterion_zero[["influence"]] <- TRUE
  f_one <- .usdt_formula(d)
  expect_identical(f_one$criteria, "c_A")
  expect_identical(f_one$dropped, "c_I")
  expect_identical(f_one$n_fixed, 3L)
  expect_identical(f_one$n_var, 4)
})

test_that("input problems are caught with actionable messages", {
  df <- make_binary()

  # A missing column produces an error.
  expect_error(long_args(df, subject_col = "missing_subject_column"), "not found")

  # A condition with more than two levels produces an error.
  bad <- df; bad$cond[1:10] <- 2L
  expect_error(long_args(bad, condition_levels = NULL), "exactly 2 distinct")

  # Incorrect role names produce an error.
  expect_error(long_args(df, condition_levels = c(a = 1, b = 0)),
               "named `signal` and `noise`")

  # An unknown task label produces an error.
  expect_error(long_args(df, task_levels = c(awareness = "A", influence = "Z")),
               "no rows")
  expect_error(long_args(df, task_levels = c(awareness = "A", influence = "A")),
               "different values")

  # Missing subject identifiers produce an error.
  bad_id <- df
  bad_id$subj[1L] <- NA
  expect_error(long_args(bad_id), "missing subject")

  # Display labels must identify two different tasks.
  expect_error(long_args(df, labels = c(awareness = "Task",
                                        influence = "Task")),
               "two different names")

  # Aggregated data cannot be dichotomized.
  df$one <- 1L
  expect_error(
    usdt_data_long(df, task_col = "task",
                   task_levels      = c(awareness = "A", influence = "I"),
                   subject_col      = "subj",
                   condition_col    = "cond",
                   condition_levels = c(signal = 1, noise = 0),
                   successes_col    = "response", trials_col = "one",
                   successes_type   = "counts",
                   dichotomize      = "influence"),
    "needs one row per trial")
})

test_that("aggregated values are validated before conversion", {
  ref <- long_args(make_binary())
  cnt <- ref$agg
  cnt$tk <- as.character(cnt$task)
  prepare <- function(x, type = "counts", successes = "y") {
    usdt_data_long(x, task_col = "tk",
                   task_levels = c(awareness = "A", influence = "I"),
                   subject_col = "subj",
                   condition_col = "sig",
                   condition_levels = c(signal = TRUE, noise = FALSE),
                   successes_col = successes, trials_col = "n",
                   successes_type = type)
  }

  fractional <- cnt
  fractional$y[1L] <- 0.5
  expect_error(prepare(fractional), "whole numbers")

  negative <- cnt
  negative$y[1L] <- -1
  expect_error(prepare(negative), "negative")

  too_many <- cnt
  too_many$y[1L] <- too_many$n[1L] + 1L
  expect_error(prepare(too_many), "exceeds")

  bad_trials <- cnt
  bad_trials$n[1L] <- 0
  expect_error(prepare(bad_trials), "positive values")

  ambiguous <- cnt
  ambiguous$y <- as.integer(ambiguous$y > 0)
  expect_error(prepare(ambiguous, type = "auto"), "format is ambiguous")

  proportions <- cnt
  proportions$p <- as.numeric(proportions$sig)
  out <- prepare(proportions, type = "proportions", successes = "p")
  expect_identical(out$agg$y, as.integer(proportions$p * proportions$n))

  wide <- stats::reshape(cnt[, c("subj", "tk", "sig", "y", "n")],
                         idvar = c("subj", "tk"), timevar = "sig",
                         direction = "wide")
  names(wide) <- c("subj", "task", "fa", "n_noise", "hit", "n_signal")
  wide$miss <- wide$n_signal - wide$hit
  wide$cr <- wide$n_noise - wide$fa
  wide$hit[1L] <- 0.5
  expect_error(
    usdt_data_long(wide, task_col = "task",
                   task_levels = c(awareness = "A", influence = "I"),
                   subject_col = "subj",
                   sdt_cols = c(hit = "hit", miss = "miss",
                                fa = "fa", cr = "cr")),
    "whole numbers")
})

test_that("a factor cannot be used as a response time", {
  set.seed(7)
  df <- make_trials(n_subj = 20L, n_trials = 60L, rt = TRUE)
  awareness <- df[df$task == "A", ]
  influence <- df[df$task == "I", ]
  influence$rt <- factor(influence$rt)
  expect_error(
    usdt_data_tasks(
      awareness = awareness, influence = influence, subject_col = "subj",
      condition_col = "cond", condition_levels = c(signal = 1, noise = 0),
      response_col = c(awareness = "response", influence = "rt"),
      response_levels = list(awareness = c(signal = 1, noise = 0),
                             influence = c(signal = "faster",
                                           noise = "slower")),
      dichotomize = "influence"),
    "must be numeric")
})

test_that("subjects missing from one task are reported", {
  df   <- make_binary()
  drop <- df$task == "I" & df$subj %in% 1:3
  expect_warning(long_args(df[!drop, ]), "only one task")
})

test_that("missing signal and noise roles are an error, never a guess", {
  df <- make_binary()
  expect_error(long_args(df, condition_levels = NULL),
               "`condition_levels` is required")
  expect_error(long_args(df, response_levels = NULL),
               "`response_levels` is required")

  # A median-split task names the side of the median, and must still name it.
  set.seed(1)
  rt <- make_trials(n_subj = 30L, n_trials = 60L, rt = TRUE)
  expect_error(
    long_args(rt,
              response_col    = c(awareness = "response", influence = "rt"),
              response_levels = NULL,
              dichotomize     = "influence"),
    "`response_levels` is required")
})

# Per-task arguments

# This function relabels one task so the two designs differ.
make_split_tasks <- function(seed = 3L) {
  df <- make_binary(seed)
  awareness <- df[df$task == "A", ]
  influence <- df[df$task == "I", ]
  names(awareness)[names(awareness) == "cond"] <- "seen"
  names(influence)[names(influence) == "cond"] <- "cue"
  awareness$seen <- ifelse(awareness$seen == 1, "old", "new")
  influence$cue  <- ifelse(influence$cue == 1, "cued", "uncued")
  list(awareness = awareness, influence = influence)
}

test_that("each task may use its own condition column and levels", {
  p <- make_split_tasks()
  d <- usdt_data_tasks(
    awareness = p$awareness, influence = p$influence, subject_col = "subj",
    condition_col    = list(awareness = "seen", influence = "cue"),
    condition_levels = list(awareness = c(signal = "old",  noise = "new"),
                            influence = c(signal = "cued", noise = "uncued")),
    response_col     = "response",
    response_levels  = c(signal = 1, noise = 0))

  expect_identical(d$meta$tasks$awareness$condition_col, "seen")
  expect_identical(d$meta$tasks$influence$condition_col, "cue")
  expect_identical(unname(d$meta$tasks$influence$condition_levels[["signal"]]),
                   "cued")

  # The relabelling changes no count, so the aggregate must not move.
  ref <- long_args(make_binary(3L))
  expect_equal(d$agg[, c("subj", "task", "sig", "y", "n")],
               ref$agg[, c("subj", "task", "sig", "y", "n")])
})

test_that("each task may use its own subject column and input format", {
  p <- make_split_tasks()
  tab <- do.call(rbind, lapply(split(p$influence, p$influence$subj), function(z)
    data.frame(id = z$subj[1L],
               H  = sum(z$cue == "cued"   & z$response == 1),
               M  = sum(z$cue == "cued"   & z$response == 0),
               FA = sum(z$cue == "uncued" & z$response == 1),
               CR = sum(z$cue == "uncued" & z$response == 0))))
  d <- usdt_data_tasks(
    awareness = p$awareness, influence = tab,
    subject_col      = list(awareness = "subj", influence = "id"),
    condition_col    = list(awareness = "seen", influence = NULL),
    condition_levels = list(awareness = c(signal = "old", noise = "new"),
                            influence = NULL),
    response_col     = list(awareness = "response", influence = NULL),
    response_levels  = list(awareness = c(signal = 1, noise = 0),
                            influence = NULL),
    sdt_cols         = list(awareness = NULL,
                            influence = c(hit = "H", miss = "M",
                                          fa = "FA", cr = "CR")))

  expect_identical(d$meta$tasks$influence$subject, "id")
  expect_identical(d$meta$tasks$influence$granularity, "SDT table")
  expect_identical(d$meta$tasks$awareness$granularity, "trial level")

  ref <- long_args(make_binary(3L))
  expect_equal(d$agg$y, ref$agg$y)
  expect_equal(d$agg$n, ref$agg$n)
})

test_that("dichotomize accepts task names and one logical per task", {
  df <- make_binary(5L, n_trials = 40L)
  df$rt <- stats::rnorm(nrow(df))
  args <- list(
    response_col    = list(awareness = "response", influence = "rt"),
    response_levels = list(awareness = c(signal = 1, noise = 0),
                           influence = c(signal = "faster", noise = "slower")))

  by_name <- do.call(long_args, c(list(df, dichotomize = "influence"), args))
  by_flag <- do.call(long_args,
                     c(list(df, dichotomize = list(awareness = FALSE,
                                                   influence = TRUE)), args))
  expect_equal(by_name$agg, by_flag$agg)
  expect_false(by_flag$meta$tasks$awareness$dichotomized)
  expect_true(by_flag$meta$tasks$influence$dichotomized)

  both <- do.call(long_args, c(list(df, dichotomize = TRUE),
                               list(response_col = "rt",
                                    response_levels = c(signal = "faster",
                                                        noise = "slower"))))
  expect_true(both$meta$tasks$awareness$dichotomized)
  expect_true(both$meta$tasks$influence$dichotomized)
})

test_that("ties and successes_type may differ between tasks", {
  df <- make_binary(6L, n_trials = 40L)
  df$rt <- stats::rnorm(nrow(df))
  d <- long_args(df,
                 response_col    = "rt",
                 response_levels = c(signal = "faster", noise = "slower"),
                 dichotomize     = "both",
                 ties            = list(awareness = "noise",
                                        influence = "random"))
  expect_identical(unname(d$meta$ties), c("noise", "random"))
  expect_error(long_args(df, ties = list(awareness = "noise",
                                         influence = "half")),
               "`ties` must be one of")
})

test_that("naming only one task is an error rather than a shared value", {
  df <- make_binary(7L)
  expect_error(long_args(df, condition_col = list(awareness = "cond",
                                                  influnce = "cond")),
               "names one task but not the other")
  expect_error(long_args(df, dichotomize = "influencia"),
               "does not accept")
  expect_error(long_args(df, dichotomize = c("none", "awareness")),
               "combines")
})

test_that("the printed mapping keeps one row per task", {
  p <- make_split_tasks()
  d <- usdt_data_tasks(
    awareness = p$awareness, influence = p$influence, subject_col = "subj",
    condition_col    = list(awareness = "seen", influence = "cue"),
    condition_levels = list(awareness = c(signal = "old",  noise = "new"),
                            influence = c(signal = "cued", noise = "uncued")),
    response_col     = "response",
    response_levels  = c(signal = 1, noise = 0))
  out <- utils::capture.output(print(d))

  expect_false(any(grepl("\bboth\b", out)))
  expect_true(any(grepl("^  subject .* Awareness .* subj", out)))
  expect_true(any(grepl("^ +Influence +subj", out)))
  expect_true(any(grepl("^  condition +Awareness +seen +old +new", out)))
  expect_true(any(grepl("^ +Influence +cue +cued +uncued", out)))
})
