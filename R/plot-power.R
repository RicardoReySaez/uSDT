# plot-power.R
# Plots of planned populations, power, smallest effects and plans
# Author: Ricardo Rey-Sáez
# Last modified: 28-10-2026

# Public functions

#' Plot a population, its power, its smallest effects or a plan
#'
#' Draws the objects of the power analysis with the style of the other uSDT
#' plots.
#'
#' * A `usdt_population` shows, for each scenario, the bivariate normal
#'   distribution of the true \eqn{d'} across subjects (ellipses holding 50%
#'   and 90% of them) with the latent regression of indirect on direct
#'   \eqn{d'} and its intercept at \eqn{d'_D = 0}. With `trials`, a dashed
#'   ellipse adds the distribution of the \eqn{d'} that those trials measure,
#'   wider by the measurement error, which is what attenuates the observed
#'   correlation.
#' * A `usdt_power` shows power against the number of subjects, or of trials
#'   when the subjects are fixed, one line per hypothesis. A simulation adds
#'   its 95% Wilson intervals. A single design shows one point per
#'   hypothesis.
#' * A `usdt_mdes` shows the power along the values searched, with the
#'   smallest detectable effect and the population value marked.
#' * A `usdt_plan` with subjects and trials free shows the cost of the cheapest
#'   closed-form design for each number of trials (a line for one free task,
#'   a map for two), with the chosen design marked; with one value free, it
#'   shows the approximate power (or the expected interval width) along it.
#'
#' @param x A `usdt_population`, `usdt_power`, `usdt_mdes` or `usdt_plan`
#'   object.
#' @param trials For a population, the number of trials per subject in each
#'   task, as in [usdt_power()], to add the distribution of the measured
#'   \eqn{d'}.
#' @param signal,dichotomize For a population with `trials`, as in
#'   [usdt_power()].
#' @param ... Ignored.
#'
#' @return A `ggplot` object. Its underlying data frame is stored in `$data`.
#'
#' @seealso [usdt_population()], [usdt_power()], [usdt_mdes()], [usdt_plan()]
#'
#' @examples
#' pop <- usdt_population(dprime    = c(direct = 0.15, indirect = 0.30),
#'                        sd_dprime = c(direct = 0.35, indirect = 0.30),
#'                        cor_dprime = 0.3, criterion = 0.1)
#' plot(pop, trials = c(direct = 64, indirect = 200), dichotomize = "indirect")
#'
#' plot(usdt_power(pop, subjects = seq(20, 300, 40),
#'                 trials = c(direct = 64, indirect = 200),
#'                 dichotomize = "indirect"))
#'
#' plot(usdt_mdes(pop, subjects = 80, trials = c(direct = 64, indirect = 200),
#'                hypothesis = "intercept > 0", dichotomize = "indirect"))
#'
#' plot(usdt_plan(pop, hypothesis = "diff = 0", subjects = c(20, 300),
#'                trials = list(direct = c(32, 96), indirect = 200),
#'                cost = c(subject = 600, direct = 4, indirect = 1.5),
#'                dichotomize = "indirect"))
#'
#' @name plot_power
NULL

#' @rdname plot_power
#' @export
plot.usdt_population <- function(x, trials = NULL, signal = 0.5,
                                 dichotomize = c("none", "indirect", "direct",
                                                 "both"), ...) {
  if (length(list(...))) .usdt_stop("`...` is not used by this plot.")
  if (nrow(x) > 12L) {
    .usdt_stop("the population has ", nrow(x), " scenarios; plot at most 12, ",
               "e.g. a population made with fewer values.")
  }
  labels <- attr(x, "labels")
  scenarios <- .scenario_labels(x)
  split <- .check_dichotomize(dichotomize)
  design <- if (!is.null(trials)) .design_trials(trials, signal)
  if (!is.null(design) && nrow(design$trials) > 1L) {
    .usdt_stop("`trials` must give one design for the plot.")
  }

  ellipses <- list()
  lines <- list()
  notes <- character(0)
  for (s in seq_len(nrow(x))) {
    pop <- as.list(as.data.frame(x)[s, ])
    Sigma <- matrix(c(pop$sigma_D^2, pop$rho * pop$sigma_D * pop$sigma_I,
                      pop$rho * pop$sigma_D * pop$sigma_I, pop$sigma_I^2), 2L)
    mu <- c(pop$mu_D, pop$mu_I)
    true <- .covariance_ellipse(mu, Sigma, c(0.5, 0.9))
    true$kind <- "true"
    observed <- NULL
    if (!is.null(design)) {
      v <- .error_variance(pop, list(signal = unlist(design$signal[1L, ]),
                                     noise = unlist(design$noise[1L, ])), split)
      observed <- .covariance_ellipse(mu, Sigma + diag(v), 0.9)
      observed$kind <- "measured"
      reliability <- diag(Sigma) / (diag(Sigma) + v)
      if (s == 1L) {
        notes <- sprintf("reliability %.2f %s, %.2f %s", reliability[1L],
                         labels[["direct"]], reliability[2L],
                         labels[["indirect"]])
      }
    }
    both <- rbind(true, observed)
    both$scenario <- scenarios[s]
    ellipses[[s]] <- both
    span <- range(both$x, 0)
    lines[[s]] <- data.frame(scenario = scenarios[s],
                             x = seq(span[1L], span[2L], length.out = 50L),
                             intercept = pop$intercept, slope = pop$slope)
    lines[[s]]$y <- pop$intercept + pop$slope * lines[[s]]$x
  }
  ellipses <- do.call(rbind, ellipses)
  ellipses$scenario <- factor(ellipses$scenario, levels = scenarios)
  lines <- do.call(rbind, lines)
  lines$scenario <- factor(lines$scenario, levels = scenarios)
  points <- unique(lines[c("scenario", "intercept")])
  points$x <- 0
  centres <- data.frame(scenario = factor(scenarios, levels = scenarios),
                        x = x$mu_D, y = x$mu_I)

  model_colour <- "#116B60"
  model_fill <- "#3FA88E"
  plot <- ggplot2::ggplot(ellipses) +
    ggplot2::geom_hline(yintercept = 0, colour = "#555B60",
                        linetype = "dotted", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, colour = "#555B60",
                        linetype = "dotted", linewidth = 0.5) +
    ggplot2::geom_polygon(
      data = ellipses[ellipses$kind == "true", ],
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   group = .data[["coverage"]]),
      fill = model_fill, alpha = 0.12, colour = model_colour, linewidth = 0.6)
  if (any(ellipses$kind == "measured")) {
    plot <- plot + ggplot2::geom_path(
      data = ellipses[ellipses$kind == "measured", ],
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = "#687076", linetype = "dashed", linewidth = 0.6)
  }
  caption <- .wrap_caption(
    "Ellipses hold 50% and 90% of the subjects' true d'. The line is the ",
    "latent regression of indirect on direct d' (H3), and the open circle ",
    "its intercept at a direct d' of zero; the filled point is the mean. ",
    if (length(notes)) paste0(
      "The dashed ellipse holds 90% of the d' that the trials measure (",
      notes, "): the measurement error widens it and weakens the observed ",
      "correlation.") else "")
  plot +
    ggplot2::geom_line(
      data = lines, ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = model_colour, linewidth = 1) +
    ggplot2::geom_point(
      data = centres, ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = model_colour, size = 2.4) +
    ggplot2::geom_point(
      data = points, ggplot2::aes(x = .data[["x"]], y = .data[["intercept"]]),
      shape = 21, size = 3.2, stroke = 1.05, colour = model_colour,
      fill = "white") +
    ggplot2::facet_wrap(ggplot2::vars(.data[["scenario"]])) +
    ggplot2::labs(
      x = bquote(.(labels[["direct"]]) ~ italic(d) * minute),
      y = bquote(.(labels[["indirect"]]) ~ italic(d) * minute),
      caption = caption) +
    .power_theme() +
    ggplot2::theme(aspect.ratio = 1,
                   strip.background = if (nrow(x) == 1L)
                     ggplot2::element_blank() else
                       ggplot2::element_rect(fill = "#F5F7F6", colour = NA),
                   strip.text = if (nrow(x) == 1L) ggplot2::element_blank() else
                     ggplot2::element_text(colour = "#252A2D", face = "bold",
                                           size = 11))
}

#' @rdname plot_power
#' @export
plot.usdt_power <- function(x, ...) {
  if (length(list(...))) .usdt_stop("`...` is not used by this plot.")
  labels <- attr(x, "labels")
  d <- as.data.frame(x)
  d$label <- .ascii_hypothesis(d$hypothesis, d$quantity)
  d$label <- factor(d$label, levels = unique(d$label))
  simulated <- identical(attr(x, "method"), "simulation")

  # The axis is the number of subjects when it varies, otherwise the trials
  # of the task that varies; the rest of each design names the panel.
  varies <- vapply(c("subjects", "trials_direct", "trials_indirect"),
                   function(k) length(unique(d[[k]])) > 1L, TRUE)
  axis <- if (varies[["subjects"]]) "subjects" else
    if (varies[["trials_direct"]]) "trials_direct" else
      if (varies[["trials_indirect"]]) "trials_indirect" else NA_character_
  axis_title <- switch(axis %||% "none",
    subjects = "Subjects",
    trials_direct = paste(labels[["direct"]], "trials per subject"),
    trials_indirect = paste(labels[["indirect"]], "trials per subject"),
    "Power")
  others <- setdiff(c("subjects", "trials_direct", "trials_indirect"), axis)
  d$panel <- .design_panel(d, others, labels,
                           length(unique(d$scenario)) > 1L)

  colours <- .hypothesis_colours(levels(d$label))
  caption <- .wrap_caption(
    if (simulated) paste0("Simulated power from ", .fmt_int(attr(x, "nsim")),
                          " studies per design, with 95% Wilson intervals. ")
    else "Approximate power, as in usdt_power(). ",
    "The dashed line marks a power of 0.80",
    if (any(d$in_null)) paste0("; a hypothesis whose population lies in H0 ",
                               "shows its type I error rate") else "",
    if (any(d$estimable < 0.99, na.rm = TRUE))
      ". Fits on the boundary count as failures to reject" else "",
    ".")

  if (is.na(axis)) {
    plot <- ggplot2::ggplot(d, ggplot2::aes(x = .data[["power"]],
                                            y = .data[["label"]],
                                            colour = .data[["label"]])) +
      ggplot2::geom_vline(xintercept = 0.8, colour = "#697076",
                          linetype = "dashed", linewidth = 0.5)
    if (simulated) {
      plot <- plot + ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data[["power.low"]], xmax = .data[["power.high"]]),
        width = 0.2, linewidth = 0.7, orientation = "y")
    }
    return(plot +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_colour_manual(values = colours, guide = "none") +
      ggplot2::scale_x_continuous(limits = c(0, 1),
                                  expand = ggplot2::expansion(mult = 0.02)) +
      ggplot2::scale_y_discrete(limits = rev(levels(d$label))) +
      ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]])) +
      ggplot2::labs(x = "Power", y = NULL, caption = caption) +
      .power_theme() + .single_panel(d$panel))
  }

  plot <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[axis]],
                                          y = .data[["power"]],
                                          colour = .data[["label"]])) +
    ggplot2::geom_hline(yintercept = 0.8, colour = "#697076",
                        linetype = "dashed", linewidth = 0.5)
  if (simulated) {
    plot <- plot +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[["power.low"]], ymax = .data[["power.high"]]),
        width = 0, linewidth = 0.6, alpha = 0.7) +
      ggplot2::geom_point(size = 1.8)
  }
  plot +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]])) +
    ggplot2::labs(x = axis_title, y = "Power", caption = caption) +
    .power_theme() + .single_panel(d$panel)
}

#' @rdname plot_power
#' @export
plot.usdt_mdes <- function(x, ...) {
  if (length(list(...))) .usdt_stop("`...` is not used by this plot.")
  labels <- attr(x, "labels")
  d <- as.data.frame(x)
  several <- nrow(unique(d[c("scenario", "subjects", "trials_direct",
                             "trials_indirect")])) > 1L
  d$panel <- .ascii_hypothesis(d$hypothesis, d$quantity)
  if (several) {
    d$panel <- paste0(d$panel, "\n", .design_panel(
      d, c("subjects", "trials_direct", "trials_indirect"), labels,
      length(unique(d$scenario)) > 1L))
  }
  d$panel <- factor(d$panel, levels = unique(d$panel))
  curves <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    cbind(panel = d$panel[i], attr(x, "curve")[[i]])
  }))
  marks <- d[is.finite(d$mdes), , drop = FALSE]
  target <- attr(x, "power")
  model_colour <- "#116B60"

  ggplot2::ggplot(curves, ggplot2::aes(x = .data[["value"]],
                                       y = .data[["power"]])) +
    ggplot2::geom_hline(yintercept = target, colour = "#697076",
                        linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(data = marks,
                        ggplot2::aes(xintercept = .data[["mdes"]]),
                        colour = "#B6543A", linetype = "dashed",
                        linewidth = 0.6) +
    ggplot2::geom_line(colour = model_colour, linewidth = 1) +
    ggplot2::geom_point(
      data = d[d$varies != "margin", , drop = FALSE],
      ggplot2::aes(x = .data[["population"]],
                   y = .data[["power_population"]]),
      shape = 21, size = 3, stroke = 1, colour = model_colour, fill = "white") +
    ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]]), scales = "free_x") +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::labs(
      x = "Value searched", y = "Power",
      caption = .wrap_caption(
        "Approximate power along the values searched, from H0 outwards. The ",
        "dashed horizontal line is the power to reach (", format(target),
        "), the dashed vertical line the smallest value that reaches it, and ",
        "the open circle the population value. For an equivalence test the ",
        "value is the half-width of the region.")) +
    .power_theme()
}

#' @rdname plot_power
#' @export
plot.usdt_plan <- function(x, ...) {
  if (length(list(...))) .usdt_stop("`...` is not used by this plot.")
  st <- x$settings
  if (st$mode == "budget") {
    .usdt_stop("plot() draws plans with a power or precision target. For a ",
               "budget, plot usdt_power() over the designs around the chosen ",
               "one.")
  }
  labels <- attr(st$population, "labels")
  free <- st$space$free
  trials_free <- c("direct", "indirect")[free[c("direct", "indirect")]]
  chosen <- x$design[!is.na(x$design$subjects), , drop = FALSE]
  if (!nrow(chosen)) .usdt_stop("no design of the plan reaches its target.")
  several <- length(unique(x$design$scenario)) > 1L
  panel_of <- function(d) {
    p <- .ascii_hypothesis(d$hypothesis)
    if (several) p <- paste0("Scenario ", d$scenario, ": ", p)
    p
  }
  chosen$panel <- panel_of(chosen)
  model_colour <- "#116B60"
  axis_title <- function(k) {
    if (st$space$linked) return("Trials per subject in each task")
    paste(labels[[k]], "trials per subject")
  }

  # With subjects and trials free, the cost of the cheapest design for each
  # number of trials: from the closed form everywhere, and from the draws of
  # usdt_power() where the search checked it.
  if (free[["subjects"]] && length(trials_free)) {
    f <- x$frontier[is.finite(x$frontier$cost), , drop = FALSE]
    f$panel <- panel_of(f)
    checked <- x$checked[is.finite(x$checked$cost), , drop = FALSE]
    checked$panel <- panel_of(checked)
    if (length(trials_free) == 1L) {
      k <- paste0("trials_", trials_free)
      return(ggplot2::ggplot(f, ggplot2::aes(x = .data[[k]],
                                             y = .data[["cost"]])) +
        ggplot2::geom_line(colour = "#A8AEB3", linewidth = 0.8) +
        ggplot2::geom_point(data = checked, colour = model_colour, size = 1.6,
                            alpha = 0.8) +
        ggplot2::geom_point(data = chosen, shape = 21, size = 3.4,
                            stroke = 1.1, colour = "#B6543A", fill = "white") +
        ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]]), scales = "free_y") +
        ggplot2::labs(
          x = axis_title(trials_free), y = "Cost",
          caption = .wrap_caption(
            "The grey line is the cost of the cheapest design with each ",
            "number of trials in the closed form; the green points, that of ",
            "the designs the search checked with the draws of usdt_power(); ",
            "the circle, the chosen design. The closed form is optimistic ",
            "where boundary fits are frequent.")) +
        .power_theme())
    }

    # Each panel shows its cost relative to its cheapest design, so plans
    # whose costs differ by far share one scale; costs twice the cheapest
    # and more take the darkest colour. The grid of trials may be thinned
    # unevenly, so each design is a tile.
    f$relative <- f$cost / stats::ave(f$cost, f$panel, FUN = min)
    return(ggplot2::ggplot(f, ggplot2::aes(x = .data[["trials_direct"]],
                                           y = .data[["trials_indirect"]])) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data[["relative"]])) +
      ggplot2::geom_point(data = checked, colour = "#252A2D", size = 0.7,
                          alpha = 0.6) +
      ggplot2::geom_point(data = chosen, shape = 21, size = 3.4, stroke = 1.1,
                          colour = "#B6543A", fill = "white") +
      ggplot2::scale_fill_gradient(low = "#E3F1EC", high = "#0B4A42",
                                   limits = c(1, 2), oob = .squish,
                                   name = "Cost /\ncheapest") +
      ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]])) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(
        x = axis_title("direct"), y = axis_title("indirect"),
        caption = .wrap_caption(
          "Cost of the cheapest design with each number of trials in the ",
          "closed form, relative to the cheapest of the panel; the steps ",
          "come from whole numbers of subjects. Dark points are the designs ",
          "the search checked with the draws of usdt_power(), and the circle ",
          "the chosen one, which can lie away from the closed-form minimum ",
          "where boundary fits are frequent.")) +
      .power_theme() +
      ggplot2::theme(legend.position = "right",
                     legend.key.height = grid::unit(28, "pt"),
                     legend.key.width = grid::unit(12, "pt")))
  }

  # With one value free, the power (or width) along it, from usdt_power().
  axis <- if (free[["subjects"]]) "subjects" else paste0("trials_", trials_free)
  curves <- .plan_curves(x, axis)
  curves$panel <- panel_of(curves)
  chosen$value <- chosen[[axis]]
  chosen$y <- if (st$mode == "precision") chosen$width else chosen$power
  target <- if (st$mode == "precision") NULL else st$power
  ggplot2::ggplot(curves, ggplot2::aes(x = .data[["value"]],
                                       y = .data[["y"]])) +
    (if (!is.null(target)) ggplot2::geom_hline(
      yintercept = target, colour = "#697076", linetype = "dashed",
      linewidth = 0.5)) +
    ggplot2::geom_line(colour = model_colour, linewidth = 1) +
    ggplot2::geom_point(data = chosen, shape = 21, size = 3.2, stroke = 1.05,
                        colour = "#B6543A", fill = "white") +
    ggplot2::facet_wrap(ggplot2::vars(.data[["panel"]]),
                        scales = if (is.null(target)) "free" else "free_x") +
    (if (!is.null(target)) ggplot2::scale_y_continuous(
      limits = c(0, 1), expand = ggplot2::expansion(mult = 0.02))) +
    ggplot2::labs(
      x = if (axis == "subjects") "Subjects" else
        axis_title(sub("trials_", "", axis)),
      y = if (st$mode == "precision") "Expected interval width" else "Power",
      caption = .wrap_caption(
        if (st$mode == "precision")
          "Expected width of the interval of usdt_tests() at the population values"
        else paste0("Approximate power, as in usdt_power(); the dashed line ",
                    "is the power to reach"),
        ". The point is the chosen design.")) +
    .power_theme()
}

# Internal functions

# This function gives the theme of the power plots: that of the faceted plots,
# with the y axis back, since power and cost have a scale.
.power_theme <- function() {
  .panel_theme() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_line(colour = "#292D30", linewidth = 0.65),
      axis.ticks.y = ggplot2::element_line(colour = "#292D30", linewidth = 0.5),
      axis.text.y = ggplot2::element_text(colour = "#3E4347", size = 10),
      axis.title.y = ggplot2::element_text(colour = "#292D30", size = 13.5,
                                           margin = ggplot2::margin(r = 10)))
}

# This function keeps values inside the limits of a scale, as scales::squish()
# does, without importing scales.
.squish <- function(x, range = c(0, 1), only.finite = TRUE) {
  pmin(pmax(x, range[1L]), range[2L])
}

# This function hides the strip of a plot with a single unnamed panel.
.single_panel <- function(panel) {
  if (length(unique(panel)) == 1L && identical(as.character(panel[1L]), "")) {
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_blank())
  }
}

# The colours of the hypotheses, in the order they are given: the model green
# first, then the colours of the other plots.
.hypothesis_colours <- function(labels) {
  palette <- c("#116B60", "#B6543A", "#3A6EA5", "#8A6BB0", "#C28F2C",
               "#5E6B73")
  stats::setNames(rep_len(palette, length(labels)), labels)
}

# This function writes hypothesis labels for a graphics device, which may not
# draw the delta: "Delta d'" stands for the difference.
.ascii_hypothesis <- function(hypothesis, quantity = NULL) {
  delta <- .usdt_chars()$delta
  out <- gsub(paste0(delta, "d'"), "Delta d'", hypothesis, fixed = TRUE)
  out[hypothesis == "all"] <- "All at once"
  out
}

# This function names the panel of each design by what the axis does not
# show: the scenario, the subjects or the trials of each task.
.design_panel <- function(d, columns, labels, scenarios) {
  parts <- lapply(columns, function(k) {
    switch(k,
      subjects = paste(d$subjects, "subjects"),
      trials_direct = paste(d$trials_direct, labels[["direct"]]),
      trials_indirect = paste(d$trials_indirect, labels[["indirect"]]))
  })
  text <- do.call(paste, c(parts, sep = ", "))
  if (scenarios) text <- paste0("Scenario ", d$scenario, ": ", text)
  if (!length(columns) && !scenarios) text <- rep("", nrow(d))
  factor(text, levels = unique(text))
}

# This function names each scenario of a population by the values that differ
# between them.
.scenario_labels <- function(x) {
  if (nrow(x) == 1L) return("")
  inputs <- c(mu_D = "mean d' D", mu_I = "mean d' I", sigma_D = "sd(d') D",
              sigma_I = "sd(d') I", rho = "cor(d')", intercept = "intercept",
              criterion_D = "mean c D", criterion_I = "mean c I",
              sd_criterion_D = "sd(c) D", sd_criterion_I = "sd(c) I",
              cor_criterion = "cor(c)")
  if (!identical(attr(x, "parameterization"), "regression")) {
    inputs <- inputs[names(inputs) != "intercept"]
  }
  varies <- names(inputs)[vapply(names(inputs), function(k) {
    length(unique(x[[k]])) > 1L
  }, TRUE)]
  out <- vapply(seq_len(nrow(x)), function(i) {
    paste(paste(inputs[varies], formatC(unlist(x[i, varies]), format = "f",
                                        digits = 2L)), collapse = ", ")
  }, "")
  make.unique(out, sep = " #")
}

# This function gives the points of the ellipses of a bivariate normal that
# hold the given shares of its mass.
.covariance_ellipse <- function(mu, Sigma, levels, n = 121L) {
  e <- eigen(Sigma, symmetric = TRUE)
  transform <- diag(sqrt(pmax(e$values, 0)), 2L) %*% t(e$vectors)
  angle <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(angle), sin(angle))
  do.call(rbind, lapply(levels, function(level) {
    points <- sqrt(stats::qchisq(level, df = 2L)) * circle %*% transform
    data.frame(x = points[, 1L] + mu[1L], y = points[, 2L] + mu[2L],
               coverage = factor(level, levels = levels))
  }))
}

# This function evaluates a plan with one free value along that value: the
# approximate power of usdt_power() for a power target, or the expected width
# of the interval for a precision target.
.plan_curves <- function(x, axis) {
  st <- x$settings
  sp <- st$space
  d <- x$design[!is.na(x$design$subjects), , drop = FALSE]
  # The subjects of each design run up to three times its chosen number,
  # where the curve has long levelled off.
  if (axis != "subjects") {
    values <- unique(sp$grid[[sub("trials_", "", axis)]])
    values <- values[unique(round(seq(1, length(values),
                                      length.out = min(60L, length(values)))))]
  }
  rows <- lapply(seq_len(nrow(d)), function(i) {
    pop <- .population_rows(st$population, d$scenario[i])
    subjects <- d$subjects[i]
    trials <- list(direct = d$trials_direct[i], indirect = d$trials_indirect[i])
    if (axis == "subjects") {
      top <- min(max(sp$subjects), max(3 * subjects, min(sp$subjects) + 20))
      subjects <- unique(round(seq(min(sp$subjects), top, length.out = 60L)))
    } else if (sp$linked) {
      trials <- values
    } else {
      trials[[sub("trials_", "", axis)]] <- values
    }
    if (st$mode == "precision") {
      q <- sub(" width.*$", "", d$hypothesis[i])
      p <- usdt_power(pop, subjects = subjects, trials = trials,
                      hypothesis = if (q == "rho") "rho in [-0.5, 0.5]" else
                        paste(q, "= 0"),
                      alpha = st$alpha, signal = st$signal,
                      dichotomize = st$split)
      z <- stats::qnorm(1 - st$alpha / 2)
      y <- if (q == "rho") {
        rho <- p$true
        se_z <- p$se / (1 - rho^2)
        tanh(atanh(rho) + z * se_z) - tanh(atanh(rho) - z * se_z)
      } else 2 * z * p$se
    } else {
      hypothesis <- if (d$hypothesis[i] == "all at once") st$hypothesis else
        st$hypothesis[match(d$hypothesis[i], vapply(st$tests, function(h)
          .bf_label(.bf_name(h$quantity), h), ""))]
      p <- usdt_power(pop, subjects = subjects, trials = trials,
                      hypothesis = hypothesis, alpha = st$alpha,
                      signal = st$signal, dichotomize = st$split,
                      combine = if (d$hypothesis[i] == "all at once") "all"
                        else "each")
      if (d$hypothesis[i] == "all at once") p <- p[p$hypothesis == "all", ]
      y <- p$power
    }
    data.frame(scenario = d$scenario[i], hypothesis = d$hypothesis[i],
               value = p[[axis]], y = y, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
