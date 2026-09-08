# sdt-moments.R
# Calculate subject-level Signal Detection Theory measures
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Public functions

#' Signal detection measures for each subject
#'
#' Computes empirical hit rates, false-alarm rates, sensitivity (\eqn{d'}), and
#' response criteria for each participant without fitting a model. It can also
#' calculate sampling variances, standard errors, and expected values for
#' \eqn{d'}.
#'
#' @param data A `usdt_data` object or a standard trial-level data frame.
#'   When given a `usdt_data` object, the function processes both tasks and
#'   includes a `task` column in the output.
#' @param subject_col,condition_col,response_col Column names for subject,
#'   condition, and response variables. Only required when `data` is a plain
#'   data frame.
#' @param condition_levels,response_levels Named vectors mapping condition and
#'   response labels, like `c(signal = "old", noise = "new")`. Required for a
#'   plain data frame. A `usdt_data` object supplies its own roles and needs
#'   neither.
#' @param coding Criterion definition to report: `"deviation"` measures the
#'   criterion from the midpoint between the signal and noise distributions,
#'   whereas `"treatment"` measures it from the noise distribution.
#'   A `usdt_data` object supplies its own coding.
#' @param correction Handling of extreme rates (0 or 1) that make \eqn{d'}
#'   infinite. `"hautus"` adds 0.5 to all four cell counts for affected
#'   participants. `"none"` leaves infinite values in place.
#' @param variances Logical. If `TRUE`, computes the sampling variance of
#'   \eqn{d'} from Gourevitch and Galanter (1967) and Miller (1996), each with
#'   its own standard error, as well as the expected value of \eqn{d'} under
#'   Miller's distribution.
#'
#' @return A data frame with one row per subject (or per subject and task for
#'   `usdt_data` inputs) containing:
#' * `hit`, `miss`, `fa`, `cr`: Raw response counts.
#' * `hr`, `far`: Observed hit and false-alarm rates.
#' * `zhr`, `zfar`: Probit-transformed rates.
#' * `dprime`, `criterion`: Descriptive SDT estimates.
#' * `corrected`: Logical flag indicating whether the participant received an
#'   edge correction.
#' * `var_gg`, `se_gg`: Asymptotic variance and standard error from Gourevitch
#'   and Galanter (1967), present when `variances = TRUE`.
#' * `var_miller`, `se_miller`, `expected_dprime`: Moments from Miller (1996),
#'   present when `variances = TRUE`.
#'
#' @details
#' Edge corrections apply only to participants with extreme rates (0 or 1)
#' rather than the whole sample, leaving well-defined rates unchanged.
#'
#' When requested, the sampling variance of \eqn{d'} is estimated using the
#' asymptotic approximation of Gourevitch and Galanter (1967) and the
#' binomial-distribution approach of Miller (1996). See Suero et al. (2017)
#' for a comparison between the two approaches.
#'
#' @references
#' Gourevitch, V., & Galanter, E. (1967). A significance test for one parameter
#' isosensitivity functions. \emph{Psychometrika}, 32(1), 25--33.
#' \doi{10.1007/BF02289402}
#'
#' Hautus, M. J. (1995). Corrections for extreme proportions and their biasing
#' effects on estimated values of \eqn{d'}. \emph{Behavior Research Methods,
#' Instruments, & Computers}, 27(1), 46--51. \doi{10.3758/BF03203619}
#'
#' Miller, J. (1996). The sampling distribution of \eqn{d'}. \emph{Perception &
#' Psychophysics}, 58(1), 65--72. \doi{10.3758/BF03205476}
#'
#' Suero, M., Privado, J., & Botella, J. (2017). Methods to estimate the
#' variance of some indices of the signal detection theory: A simulation study.
#' \emph{Psicologica}, 38(1), 77--109.
#'
#' @seealso [hsdt()]
#'
#' @examples
#' # 1. From a prepared usdt_data object (both tasks at once)
#' d <- usdt_data_tasks(
#'   direct   = vadillo_awareness,
#'   indirect = vadillo_cuing,
#'   subject_col      = "subj",
#'   condition_col    = "condition",
#'   condition_levels = c(signal = "old", noise = "new"),
#'   response_col     = list(direct = "judged.old", indirect = "rt"),
#'   response_levels  = list(direct   = c(signal = 1, noise = 0),
#'                           indirect = c(signal = "faster", noise = "slower")),
#'   dichotomize      = list(direct = FALSE, indirect = TRUE)
#' )
#'
#' head(sdt_moments(d))
#'
#' # 2. From raw trials with sampling variances and standard errors
#' head(sdt_moments(vadillo_awareness,
#'                  subject_col      = "subj",
#'                  condition_col    = "condition",
#'                  condition_levels = c(signal = "old", noise = "new"),
#'                  response_col     = "judged.old",
#'                  response_levels  = c(signal = 1, noise = 0),
#'                  variances        = TRUE))
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
    .require_levels(data[[condition_col]], "condition_levels",
                    paste0("`", condition_col, "`"))
  }
  response_levels <- .check_levels(response_levels, "response_levels")
  if (is.null(response_levels)) {
    .require_levels(data[[response_col]], "response_levels",
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
    out <- cbind(out, .sdt_variances(nr = nr, ns = ns, pi_fa = far, pi_a = hr))
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
  # Each variance travels next to its own standard error, so that the choice
  # between the two is made by name at the point of use.
  data.frame(var_gg          = var_gg,
             se_gg           = sqrt(var_gg),
             var_miller      = var_m,
             se_miller       = sqrt(var_m),
             expected_dprime = v_esp)
}
