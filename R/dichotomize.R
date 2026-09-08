# dichotomize.R
# Dichotomize continuous measures (e.g., response times) into binary responses
# Author: Ricardo Rey-Sáez
# Last modified: 08-09-2026

# Public functions

#' Dichotomize response times into binary choices
#'
#' Splits response times (or other continuous measures) at each subject's
#' overall median, following the preprocessing approach of Meyen et al. (2022).
#' The median is calculated across all trials for each participant without
#' distinguishing between stimulus conditions or other covariates. This produces
#' a binary outcome that allows response times to be mapped onto a Signal
#' Detection Theory sensitivity metric (\eqn{d'}).
#'
#' @param x Numeric vector of continuous values, typically response times.
#' @param by Vector identifying the subject for each observation in `x`.
#'   Medians are computed independently for each participant.
#' @param signal Character string specifying which side of the median will be
#'   treated as the "signal" response under an SDT framework. Use `"faster"`
#'   when the target condition speeds up responses (e.g., facilitatory priming,
#'   spatial cueing) or `"slower"` when it slows responses down (e.g.,
#'   interference, Stroop-like effects).
#' @param ties How to handle trials that match the subject's median exactly.
#'   `"noise"` assigns them to the noise category (0). `"random"` breaks ties at
#'   random, keeping cell proportions as balanced as possible.
#'
#' @return An integer vector of `0`s (noise response) and `1`s (signal response)
#'   matching the length of `x`. Missing values (`NA`) are preserved.
#'
#' @details
#' With an odd number of trials (\eqn{n}), a dataset cannot be split into two
#' equal halves because the median falls exactly on an observed trial.
#' Setting `ties = "noise"` assigns this middle trial to noise, producing a
#' signal proportion of \eqn{(n - 1) / (2\cdot n)} and slightly shifting the response
#' criterion. In practice, this difference (\eqn{1 / (2\cdot n)}) is negligible, but
#' setting `ties = "random"` resolves ties probabilistically to avoid any
#' systematic directional bias.
#'
#' @references
#' Meyen, S., Zerweck, I. A., Amado, C., von Luxburg, U., & Franz, V. H.
#' (2022). Advancing research on unconscious priming: When can scientists claim
#' an indirect task advantage? \emph{Journal of Experimental Psychology: General},
#' 151(1), 65--81. \doi{10.1037/xge0001065}
#'
#' @examples
#' rt   <- c(320, 410, 295, 500, 380, 450)
#' subj <- rep(c("s1", "s2"), each = 3)
#' meyen_split(rt, by = subj)
#'
#' @export
meyen_split <- function(x, by,
                        signal = c("faster", "slower"),
                        ties   = c("noise", "random")) {

  # The function checks the requested options and the input data.
  signal <- match.arg(signal)
  ties   <- match.arg(ties)
  if (!is.numeric(x)) {
    .usdt_stop("`x` must be numeric; `meyen_split()` needs the raw continuous ",
               "measure (e.g. response times), not an already-binary response.")
  }
  if (length(by) != length(x)) {
    .usdt_stop("`by` has length ", length(by), " but `x` has length ",
               length(x), "; they must match.")
  }
  if (all(is.na(x))) .usdt_stop("`x` is entirely missing.")

  # Each subject receives one median across both conditions.
  grp <- as.character(by)
  med <- stats::ave(x, grp, FUN = function(v) stats::median(v, na.rm = TRUE))

  # Values on the selected side become signal responses.
  out <- if (signal == "faster") as.integer(x < med) else as.integer(x > med)

  # The random option redistributes values that equal the median.
  if (ties == "random") out <- .break_ties(x, med, grp, out)

  # Missing input values remain missing.
  out[is.na(x)] <- NA_integer_
  out
}

# Internal functions

# This function assigns ties while keeping each split nearly even.
.break_ties <- function(x, med, grp, out) {

  # The function returns early when no value equals the median.
  tied <- which(!is.na(x) & x == med)
  if (!length(tied)) return(out)

  # Each subject's ties are handled together.
  for (g in unique(grp[tied])) {
    idx  <- which(grp == g & !is.na(x))
    tidx <- intersect(tied, idx)
    n    <- length(idx)

    # The target keeps the split as even as possible.
    target <- if (n %% 2L == 0L) n %/% 2L else n %/% 2L + sample.int(2L, 1L) - 1L

    # A random set of tied trials fills the target.
    out[tidx] <- 0L
    need <- max(min(target - sum(out[idx], na.rm = TRUE), length(tidx)), 0L)
    if (need > 0L) out[tidx[sample.int(length(tidx), need)]] <- 1L
  }
  out
}
