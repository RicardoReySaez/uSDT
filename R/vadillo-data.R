#' Awareness data from a probabilistic cuing experiment
#'
#' Trial-level direct-awareness data from Experiment 2 of Vadillo et al.
#' (2024). The object reproduces the source CSV file without filtering or
#' recoding and can be used as the `direct` input to [usdt_data_tasks()].
#'
#' @format A data frame with 6,656 rows, 104 participants, and 10 variables:
#' \describe{
#'   \item{`subj`}{Participant identifier.}
#'   \item{`file`}{Original participant file name.}
#'   \item{`trial`}{Trial number.}
#'   \item{`pattId`}{Pattern identifier.}
#'   \item{`judgm`}{Awareness judgement on the original response scale.}
#'   \item{`judged.old`}{Whether the pattern was judged old (1) or new (0).}
#'   \item{`condition`}{Whether the pattern was old or new.}
#'   \item{`offset`}{Spatial-offset condition.}
#'   \item{`color`}{Colour condition.}
#'   \item{`set.size`}{Set-size condition.}
#' }
#' @source Vadillo, M. A., Malejka, S., & Shanks, D. R. (2024). Mapping the
#'   reliability multiverse of contextual cuing. *Journal of Experimental
#'   Psychology: Learning, Memory, and Cognition*. \doi{10.1037/xlm0001410}.
#'   Data retrieved from \url{https://osf.io/jp3gx/}.
#' @examples
#' data(vadillo_awareness)
#' str(vadillo_awareness)
"vadillo_awareness"

#' Cuing data from a probabilistic cuing experiment
#'
#' Trial-level indirect-task data from Experiment 2 of Vadillo et al. (2024).
#' The object reproduces the source CSV file without filtering or recoding and
#' can be used as the `indirect` input to [usdt_data_tasks()].
#'
#' @format A data frame with 39,936 rows, 104 participants, and 13 variables:
#' \describe{
#'   \item{`experiment`}{Experiment label.}
#'   \item{`subj`}{Participant identifier.}
#'   \item{`file`}{Original participant file name.}
#'   \item{`trial`}{Trial number.}
#'   \item{`block`}{Block number.}
#'   \item{`epoch`}{Epoch number.}
#'   \item{`pattId`}{Pattern identifier.}
#'   \item{`acc`}{Response accuracy (1 = correct, 0 = incorrect).}
#'   \item{`rt`}{Response time in milliseconds.}
#'   \item{`condition`}{Whether the pattern was old or new.}
#'   \item{`offset`}{Spatial-offset condition.}
#'   \item{`color`}{Colour condition.}
#'   \item{`set.size`}{Set-size condition.}
#' }
#' @source Vadillo, M. A., Malejka, S., & Shanks, D. R. (2024). Mapping the
#'   reliability multiverse of contextual cuing. *Journal of Experimental
#'   Psychology: Learning, Memory, and Cognition*. \doi{10.1037/xlm0001410}.
#'   Data retrieved from \url{https://osf.io/jp3gx/}.
#' @examples
#' data(vadillo_cuing)
#' str(vadillo_cuing)
"vadillo_cuing"
