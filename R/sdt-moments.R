# sdt-moments.R
# This script calculates signal detection measures for each subject.
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Public functions

#' Signal detection measures for each subject
#'
#' Computes the hit rate, the false-alarm rate, d' and the criterion of every
#' subject, without fitting a model. It can also add the sampling variance of
#' d'. These descriptive values are useful to inspect the data before fitting,
#' and to compare with the model estimates afterwards.
#'
#' @param data A `usdt_data` object, or a plain data frame with one row per
#'   trial. With a `usdt_data` object the function processes both tasks and
#'   adds a `task` column.
#' @param subject_col,condition_col,response_col Names of the columns that hold
#'   the subject, the condition and the response. They are needed only for a
#'   plain data frame.
#' @param condition_levels,response_levels Which value plays each role, as
#'   `c(signal = ..., noise = ...)`. The function guesses them and reports its
#'   choice when they are missing.
#' @param coding Which criterion to report. `"deviation"` gives the classical
#'   criterion `c = -(z(HR) + z(FAR)) / 2`, measured from the midpoint between
#'   the two distributions. `"treatment"` gives `lambda = -z(FAR)`, measured
#'   from the noise distribution. The two are related by `lambda = c + d'/2`.
#'   A `usdt_data` object supplies its own coding.
#' @param correction What to do with rates of exactly 0 or 1, which make `d'`
#'   infinite. `"hautus"` adds 0.5 to the four counts of the affected subject.
#'   `"none"` leaves the infinite values in place.
#' @param variances If `TRUE`, adds the sampling variance of `d'` from
#'   Gourevitch and Galanter (1967) and from Miller (1996), together with the
#'   standard error that follows from the second one.
#'
#' @return A data frame with one row per subject, or one row per subject and
#'   task when `data` is a `usdt_data` object. It holds the four response
#'   counts (`hit`, `miss`, `fa`, `cr`), the two rates (`hr`, `far`) and their
#'   probit values (`zhr`, `zfar`), then `dprime`, `criterion`, and `corrected`
#'   to mark the subjects that received the edge correction. With
#'   `variances = TRUE` it also holds `var_gg`, `var_miller`, `e_miller` and
#'   `se_dprime`.
#'
#' @details
#' The edge correction applies only to the subjects that need it. Applying it
#' to the whole sample would change the estimates of every other subject as
#' well, and those estimates are already usable.
#'
#' The Miller variance treats the observed hit and false-alarm rates as
#' binomial probabilities. Samples that reach a rate of 0 or 1 receive the
#' values `0.5 / n` and `(n - 0.5) / n`. The number of trials stays the
#' original one throughout.
#'
#' `var_gg` comes from the expected information of the two probit cells, with
#' the criterion treated as a nuisance parameter. It equals the standard error
#' that a probit regression would give for `d'` if it were fitted to that
#' subject alone. [usdt_reliability()] uses the same formula, but evaluates it
#' at the rates the hierarchical model predicts instead of the observed ones.
#'
#' @references
#' Gourevitch, V., & Galanter, E. (1967). A significance test for one parameter
#' isosensitivity functions. *Psychometrika*.
#'
#' Hautus, M. J. (1995). Corrections for extreme proportions and their biasing
#' effects on estimated values of d'. *Behavior Research Methods*.
#'
#' Miller, J. (1996). The sampling distribution of d'. *Perception &
#' Psychophysics*.
#'
#' Suero, M., Privado, J., & Botella, J. (2017). Methods to estimate the
#' variance of some indices of the signal detection theory: A simulation study.
#' *Psicologica*.
#'
#' @seealso [hsdt()]
#'
#' @examples
#' set.seed(1)
#' df <- usdt_simulate(n_subj = 20, n_trials = 80)
#' head(sdt_moments(df[df$task == "D", ],
#'                  subject_col   = "subj",
#'                  condition_col = "cond",
#'                  condition_levels = c(signal = 1, noise = 0),
#'                  response_col  = "response",
#'                  response_levels  = c(signal = 1, noise = 0)))
#'
#' @export
sdt_moments <- function(data,
                        subject_col      = NULL,
                        condition_col    = NULL, condition_levels = NULL,
                        response_col     = NULL, response_levels  = NULL,
                        coding           = c("deviation", "treatment"),
                        correction       = c("hautus", "none"),
                        variances        = FALSE) {

  # The function checks the selected options.
  coding     <- match.arg(coding)
  correction <- match.arg(correction)

  # Prepared data already contain the required columns and coding.
  if (inherits(data, "usdt_data")) {
    out <- lapply(c("D", "I"), function(k) {
      cbind(task = data$meta$labels[[if (k == "D") "direct" else "indirect"]],
            .task_moments(data$agg, k, data$meta$coding,
                          correction = correction, variances = variances))
    })
    res <- do.call(rbind, out)
    rownames(res) <- NULL
    return(res)
  }

  # Plain data require the names and meanings of their columns.
  .check_df(data, "data")
  for (a in c("subject_col", "condition_col", "response_col")) {
    if (is.null(get(a))) .usdt_stop("`", a, "` is required when `data` is a ",
                                    "plain data frame.")
  }
  .check_cols(data, c(subject_col, condition_col, response_col),
              "column arguments", "`data`")
  .check_ids(data[[subject_col]], subject_col, "`data`")

  condition_levels <- .check_levels(condition_levels, "condition_levels")
  if (is.null(condition_levels)) {
    condition_levels <- .guess_levels(data[[condition_col]], "condition_levels",
                                      paste0("`", condition_col, "`"))
  }
  response_levels <- .check_levels(response_levels, "response_levels")
  if (is.null(response_levels)) {
    response_levels <- .guess_levels(data[[response_col]], "response_levels",
                                     paste0("`", response_col, "`"))
  }

  # The function maps the roles and counts the four response types.
  sig  <- .role_match(data[[condition_col]], condition_levels, condition_col, "`data`")
  resp <- .role_match(data[[response_col]],  response_levels,  response_col,  "`data`")
  subj <- as.character(data[[subject_col]])
  keep <- !is.na(sig) & !is.na(resp)
  if (any(!keep)) {
    .usdt_warn(sum(!keep), " trials dropped for missing condition or response.")
  }

  cnt <- .count_cells(subj[keep], sig[keep], resp[keep])
  empty <- cnt$hit + cnt$miss == 0 | cnt$fa + cnt$cr == 0
  if (any(empty)) {
    .usdt_stop("every subject needs signal and noise trials.")
  }
  .sdt_table(cnt$subj, cnt$hit, cnt$miss, cnt$fa, cnt$cr,
             coding = coding, correction = correction,
             variances = variances)
}

# Internal functions

# This function calculates one task's measures for each subject.
.task_moments <- function(agg, task, coding, correction = "hautus",
                          variances = FALSE) {
  a <- agg[agg$task == task, ]
  a <- a[order(a$subj, a$sig), ]
  .sdt_table(subj = a$subj[a$sig],
             hit  = a$y[a$sig],  miss = a$n[a$sig]  - a$y[a$sig],
             fa   = a$y[!a$sig], cr   = a$n[!a$sig] - a$y[!a$sig],
             coding = coding, correction = correction, variances = variances)
}

# This function counts hits, misses, false alarms and correct rejections.
.count_cells <- function(subj, sig, resp) {
  ind <- cbind(hit  = as.numeric(sig  & resp), miss = as.numeric(sig  & !resp),
               fa   = as.numeric(!sig & resp), cr   = as.numeric(!sig & !resp))
  m <- rowsum(ind, subj, reorder = TRUE)
  data.frame(subj = rownames(m), hit = m[, "hit"], miss = m[, "miss"],
             fa = m[, "fa"], cr = m[, "cr"], stringsAsFactors = FALSE,
             row.names = NULL)
}

# This function calculates rates, d' and criterion from the four counts.
.sdt_table <- function(subj, hit, miss, fa, cr, coding, correction, variances) {

  # The cell counts give the original response rates.
  ns  <- hit + miss
  nr  <- fa + cr
  hr  <- hit / ns
  far <- fa  / nr

  # The edge correction changes only subjects with an extreme rate.
  edge <- (hr %in% c(0, 1)) | (far %in% c(0, 1))
  if (correction == "hautus" && any(edge)) {
    hit[edge] <- hit[edge] + 0.5; miss[edge] <- miss[edge] + 0.5
    fa[edge]  <- fa[edge]  + 0.5; cr[edge]   <- cr[edge]   + 0.5
    hr  <- hit / (hit + miss)
    far <- fa  / (fa  + cr)
  }

  # The transformed rates give d' and the criterion.
  zhr  <- stats::qnorm(hr)
  zfar <- stats::qnorm(far)
  out  <- data.frame(
    subj = subj, hit = hit, miss = miss, fa = fa, cr = cr,
    hr = hr, far = far, zhr = zhr, zfar = zfar,
    dprime    = zhr - zfar,
    criterion = if (coding == "deviation") -(zhr + zfar) / 2 else -zfar,
    corrected = edge & correction == "hautus",
    stringsAsFactors = FALSE)

  # The function adds the sampling variance when the user requests it.
  if (variances) {
    v <- .sdt_variances(nr = nr, ns = ns, pi_fa = far, pi_a = hr)
    out <- cbind(out, v, se_dprime = sqrt(v$var_miller))
  }
  rownames(out) <- NULL
  out
}

# This function estimates the sampling variance of d' for each subject.
.sdt_variances <- function(nr, ns, pi_fa, pi_a) {

  # The rates stay away from values that produce infinite scores.
  pi_fa <- pmax(pmin(pi_fa, 1 - 1e-5), 1e-5)
  pi_a  <- pmax(pmin(pi_a,  1 - 1e-5), 1e-5)

  # This formula follows Gourevitch and Galanter (1967).
  var_gg <- pi_fa * (1 - pi_fa) / (nr * stats::dnorm(stats::qnorm(pi_fa))^2) +
            pi_a  * (1 - pi_a)  / (ns * stats::dnorm(stats::qnorm(pi_a))^2)

  # Subjects with the same trial count share one score grid.
  grid <- local({
    cache <- list()
    function(n) {
      key <- as.character(n)
      if (is.null(cache[[key]])) {
        fre <- if (n > 1) c(0.5, seq_len(n - 1L), n - 0.5) else c(0.5, n - 0.5)
        cache[[key]] <<- stats::qnorm(fre / n)
      }
      cache[[key]]
    }
  })

  # This calculation follows the distribution described by Miller (1996).
  k <- length(nr)
  v_esp <- var_m <- numeric(k)
  for (i in seq_len(k)) {
    nri <- max(nr[i], 1); nsi <- max(ns[i], 1)
    zf  <- grid(nri); za <- grid(nsi)
    pf  <- stats::dbinom(0:nri, nri, pi_fa[i]); pf <- pf / sum(pf)
    pa  <- stats::dbinom(0:nsi, nsi, pi_a[i]);  pa <- pa / sum(pa)
    ef  <- sum(zf * pf); ea <- sum(za * pa)
    v_esp[i] <- ea - ef
    var_m[i] <- (sum(za^2 * pa) - ea^2) + (sum(zf^2 * pf) - ef^2)
  }
  data.frame(var_gg = var_gg, e_miller = v_esp, var_miller = var_m)
}
