# aggregate.R
# Aggregate binary trials into binomial counts
# Author: Ricardo Rey-Sáez
# Last modified: 04-10-2026

# Public functions

#' Count binary responses by group
#'
#' Aggregates trial-level binary responses into the number of signal responses
#' (`y`) and the number of trials (`n`) for every combination of the grouping
#' columns. A binomial model fitted to these counts has the same likelihood,
#' and therefore the same estimates and posterior, as a Bernoulli model fitted
#' to the trials, and it runs much faster. Use it to prepare custom models
#' fitted with brms or lme4 (see [usdt_hypotheses]).
#'
#' @param data A data frame with one row per trial.
#' @param response Name of the binary response column, coded `0`/`1` or
#'   `FALSE`/`TRUE`.
#' @param by Names of the grouping columns. Every predictor of the model, the
#'   subject included, must be among them, so that all the trials of a row share
#'   their predictors.
#'
#' @return A data frame with one row per combination of the `by` columns that
#'   occurs in `data`, those columns, and the counts `y` and `n`.
#'
#' @seealso [usdt_hypotheses]
#'
#' @examples
#' # Direct task of Vadillo et al. (2025): "old" judgements by subject,
#' # condition and set size
#' counts <- usdt_aggregate(vadillo_awareness, response = "judged.old",
#'                          by = c("subj", "condition", "set.size"))
#' head(counts)
#'
#' # The counts keep every trial
#' sum(counts$n) == nrow(vadillo_awareness)
#'
#' @export
usdt_aggregate <- function(data, response, by) {
  .check_df(data, "data")
  .check_string(response, "response")
  if (!is.character(by) || !length(by) || anyNA(by) || anyDuplicated(by)) {
    .usdt_stop("`by` must name one or more distinct columns.")
  }
  if (any(c("y", "n", response) %in% by)) {
    .usdt_stop("`by` cannot include `y`, `n` or the response column, because ",
               "the counts take those places.")
  }
  .check_cols(data, c(response, by), "response` or `by", "`data`")

  r <- data[[response]]
  if (is.logical(r)) r <- as.integer(r)
  if (!is.numeric(r) || anyNA(r) || any(!r %in% c(0, 1))) {
    .usdt_stop("`", response, "` must hold binary responses coded 0/1 or ",
               "FALSE/TRUE, with no missing values.")
  }
  if (anyNA(data[by])) {
    .usdt_stop("the `by` columns contain missing values.")
  }
  .agg_counts(r, data[by])
}
