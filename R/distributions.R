# distributions.R
# The Stan distributions that uSDT can simulate and evaluate as priors
# Author: Ricardo Rey-Sáez
# Last modified: 15-10-2026

# Internal functions

# This function describes one family: the names of its arguments in Stan
# order, the condition they must meet (as a check and as text), its support,
# and its log density, distribution function and generator. The last three
# take the arguments as a numeric vector `a` and are vectorised over values
# inside the support.
.family <- function(args, rule, valid, support, logd, p, r) {
  list(args = args, rule = rule, valid = valid, support = support,
       logd = logd, p = p, r = r)
}

# Supports that do not depend on the arguments.
.real_line <- function(a) c(-Inf, Inf)
.positive_line <- function(a) c(0, Inf)

# Nodes and weights of the 64-point Gauss-Legendre rule on (-1, 1), from the
# eigen decomposition of the Jacobi matrix (Golub and Welsch).
.gauss_legendre <- local({
  k <- seq_len(63L)
  J <- matrix(0, 64L, 64L)
  J[cbind(k, k + 1L)] <- J[cbind(k + 1L, k)] <- k / sqrt(4 * k^2 - 1)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1L, ]^2)
})

# This function computes Owen's T(h, a), which gives the skew-normal
# distribution function, written over the angle t = atan(x) so that the
# integrand is smooth and bounded:
# T(h, a) = 1 / (2 pi) * integral from 0 to atan(a) of exp(-h^2 / (2 cos^2 t)).
.owens_t <- function(h, a) {
  theta <- atan(a)
  t <- theta * (.gauss_legendre$x + 1) / 2
  w <- theta * .gauss_legendre$w / (4 * pi)
  as.vector(exp(-outer(h^2 / 2, 1 / cos(t)^2)) %*% w)
}

# The inverse-gamma family, shared by inv_gamma and scaled_inv_chi_square.
.inv_gamma_logd <- function(x, alpha, beta) {
  alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - beta / x
}

# The continuous univariate families of Stan that suit a prior on a mean, a
# standard deviation or a correlation, with Stan's parameterisations. `lkj`
# stands for the marginal of one correlation of a 2 x 2 LKJ matrix, a beta
# distribution stretched to (-1, 1).
.distributions <- list(
  normal = .family(
    c("mu", "sigma"), "sigma > 0", function(a) a[2] > 0, .real_line,
    function(x, a) stats::dnorm(x, a[1], a[2], log = TRUE),
    function(q, a) stats::pnorm(q, a[1], a[2]),
    function(n, a) stats::rnorm(n, a[1], a[2])),
  std_normal = .family(
    character(0), NULL, function(a) TRUE, .real_line,
    function(x, a) stats::dnorm(x, log = TRUE),
    function(q, a) stats::pnorm(q),
    function(n, a) stats::rnorm(n)),
  student_t = .family(
    c("nu", "mu", "sigma"), "nu > 0 and sigma > 0",
    function(a) a[1] > 0 && a[3] > 0, .real_line,
    function(x, a) stats::dt((x - a[2]) / a[3], a[1], log = TRUE) - log(a[3]),
    function(q, a) stats::pt((q - a[2]) / a[3], a[1]),
    function(n, a) a[2] + a[3] * stats::rt(n, a[1])),
  cauchy = .family(
    c("mu", "sigma"), "sigma > 0", function(a) a[2] > 0, .real_line,
    function(x, a) stats::dcauchy(x, a[1], a[2], log = TRUE),
    function(q, a) stats::pcauchy(q, a[1], a[2]),
    function(n, a) stats::rcauchy(n, a[1], a[2])),
  double_exponential = .family(
    c("mu", "sigma"), "sigma > 0", function(a) a[2] > 0, .real_line,
    function(x, a) -log(2 * a[2]) - abs(x - a[1]) / a[2],
    function(q, a) {
      z <- (q - a[1]) / a[2]
      ifelse(z < 0, exp(z) / 2, 1 - exp(-z) / 2)
    },
    function(n, a) a[1] + a[2] * (stats::rexp(n) - stats::rexp(n))),
  logistic = .family(
    c("mu", "sigma"), "sigma > 0", function(a) a[2] > 0, .real_line,
    function(x, a) stats::dlogis(x, a[1], a[2], log = TRUE),
    function(q, a) stats::plogis(q, a[1], a[2]),
    function(n, a) stats::rlogis(n, a[1], a[2])),
  gumbel = .family(
    c("mu", "beta"), "beta > 0", function(a) a[2] > 0, .real_line,
    function(x, a) {
      z <- (x - a[1]) / a[2]
      -log(a[2]) - z - exp(-z)
    },
    function(q, a) exp(-exp(-(q - a[1]) / a[2])),
    function(n, a) a[1] - a[2] * log(stats::rexp(n))),
  skew_normal = .family(
    c("xi", "omega", "alpha"), "omega > 0", function(a) a[2] > 0, .real_line,
    function(x, a) {
      z <- (x - a[1]) / a[2]
      log(2) - log(a[2]) + stats::dnorm(z, log = TRUE) +
        stats::pnorm(a[3] * z, log.p = TRUE)
    },
    function(q, a) {
      z <- (q - a[1]) / a[2]
      stats::pnorm(z) - 2 * .owens_t(z, a[3])
    },
    function(n, a) {
      delta <- a[3] / sqrt(1 + a[3]^2)
      u <- stats::rnorm(n)
      v <- delta * u + sqrt(1 - delta^2) * stats::rnorm(n)
      a[1] + a[2] * ifelse(u >= 0, v, -v)
    }),
  exp_mod_normal = .family(
    c("mu", "sigma", "lambda"), "sigma > 0 and lambda > 0",
    function(a) a[2] > 0 && a[3] > 0, .real_line,
    function(x, a) {
      log(a[3]) + a[3] * (a[1] - x) + (a[3] * a[2])^2 / 2 +
        stats::pnorm((x - a[1]) / a[2] - a[3] * a[2], log.p = TRUE)
    },
    function(q, a) {
      stats::pnorm((q - a[1]) / a[2]) -
        exp(a[3] * (a[1] - q) + (a[3] * a[2])^2 / 2 +
              stats::pnorm((q - a[1]) / a[2] - a[3] * a[2], log.p = TRUE))
    },
    function(n, a) stats::rnorm(n, a[1], a[2]) + stats::rexp(n, a[3])),
  skew_double_exponential = .family(
    c("mu", "sigma", "tau"), "sigma > 0 and 0 < tau < 1",
    function(a) a[2] > 0 && a[3] > 0 && a[3] < 1, .real_line,
    function(x, a) {
      d <- x - a[1]
      log(2 * a[3] * (1 - a[3]) / a[2]) -
        2 / a[2] * ifelse(d < 0, -(1 - a[3]) * d, a[3] * d)
    },
    function(q, a) {
      d <- q - a[1]
      ifelse(d < 0, a[3] * exp(2 * (1 - a[3]) * d / a[2]),
             1 - (1 - a[3]) * exp(-2 * a[3] * d / a[2]))
    },
    function(n, a) {
      u <- stats::runif(n)
      ifelse(u < a[3], a[1] + a[2] / (2 * (1 - a[3])) * log(u / a[3]),
             a[1] - a[2] / (2 * a[3]) * log((1 - u) / (1 - a[3])))
    }),
  lognormal = .family(
    c("mu", "sigma"), "sigma > 0", function(a) a[2] > 0, .positive_line,
    function(x, a) stats::dlnorm(x, a[1], a[2], log = TRUE),
    function(q, a) stats::plnorm(q, a[1], a[2]),
    function(n, a) stats::rlnorm(n, a[1], a[2])),
  exponential = .family(
    "beta", "beta > 0", function(a) a[1] > 0, .positive_line,
    function(x, a) stats::dexp(x, a[1], log = TRUE),
    function(q, a) stats::pexp(q, a[1]),
    function(n, a) stats::rexp(n, a[1])),
  gamma = .family(
    c("alpha", "beta"), "alpha > 0 and beta > 0",
    function(a) a[1] > 0 && a[2] > 0, .positive_line,
    function(x, a) stats::dgamma(x, a[1], rate = a[2], log = TRUE),
    function(q, a) stats::pgamma(q, a[1], rate = a[2]),
    function(n, a) stats::rgamma(n, a[1], rate = a[2])),
  inv_gamma = .family(
    c("alpha", "beta"), "alpha > 0 and beta > 0",
    function(a) a[1] > 0 && a[2] > 0, .positive_line,
    function(x, a) .inv_gamma_logd(x, a[1], a[2]),
    function(q, a) stats::pgamma(a[2] / q, a[1], lower.tail = FALSE),
    function(n, a) 1 / stats::rgamma(n, a[1], rate = a[2])),
  weibull = .family(
    c("alpha", "sigma"), "alpha > 0 and sigma > 0",
    function(a) a[1] > 0 && a[2] > 0, .positive_line,
    function(x, a) stats::dweibull(x, a[1], a[2], log = TRUE),
    function(q, a) stats::pweibull(q, a[1], a[2]),
    function(n, a) stats::rweibull(n, a[1], a[2])),
  frechet = .family(
    c("alpha", "sigma"), "alpha > 0 and sigma > 0",
    function(a) a[1] > 0 && a[2] > 0, .positive_line,
    function(x, a) {
      z <- x / a[2]
      log(a[1]) - log(a[2]) - (1 + a[1]) * log(z) - z^(-a[1])
    },
    function(q, a) exp(-(q / a[2])^(-a[1])),
    function(n, a) a[2] * stats::rexp(n)^(-1 / a[1])),
  rayleigh = .family(
    "sigma", "sigma > 0", function(a) a[1] > 0, .positive_line,
    function(x, a) log(x) - 2 * log(a[1]) - x^2 / (2 * a[1]^2),
    function(q, a) -expm1(-q^2 / (2 * a[1]^2)),
    function(n, a) a[1] * sqrt(2 * stats::rexp(n))),
  chi_square = .family(
    "nu", "nu > 0", function(a) a[1] > 0, .positive_line,
    function(x, a) stats::dchisq(x, a[1], log = TRUE),
    function(q, a) stats::pchisq(q, a[1]),
    function(n, a) stats::rchisq(n, a[1])),
  inv_chi_square = .family(
    "nu", "nu > 0", function(a) a[1] > 0, .positive_line,
    function(x, a) stats::dchisq(1 / x, a[1], log = TRUE) - 2 * log(x),
    function(q, a) stats::pchisq(1 / q, a[1], lower.tail = FALSE),
    function(n, a) 1 / stats::rchisq(n, a[1])),
  scaled_inv_chi_square = .family(
    c("nu", "s"), "nu > 0 and s > 0", function(a) a[1] > 0 && a[2] > 0,
    .positive_line,
    function(x, a) .inv_gamma_logd(x, a[1] / 2, a[1] * a[2]^2 / 2),
    function(q, a) {
      stats::pgamma(a[1] * a[2]^2 / (2 * q), a[1] / 2, lower.tail = FALSE)
    },
    function(n, a) a[1] * a[2]^2 / stats::rchisq(n, a[1])),
  loglogistic = .family(
    c("alpha", "beta"), "alpha > 0 and beta > 0",
    function(a) a[1] > 0 && a[2] > 0, .positive_line,
    function(x, a) {
      z <- x / a[1]
      log(a[2]) - log(a[1]) + (a[2] - 1) * log(z) - 2 * log1p(z^a[2])
    },
    function(q, a) 1 / (1 + (q / a[1])^(-a[2])),
    function(n, a) {
      u <- stats::runif(n)
      a[1] * (u / (1 - u))^(1 / a[2])
    }),
  pareto = .family(
    c("y_min", "alpha"), "y_min > 0 and alpha > 0",
    function(a) a[1] > 0 && a[2] > 0, function(a) c(a[1], Inf),
    function(x, a) log(a[2]) + a[2] * log(a[1]) - (a[2] + 1) * log(x),
    function(q, a) 1 - (a[1] / q)^a[2],
    function(n, a) a[1] * exp(stats::rexp(n) / a[2])),
  pareto_type_2 = .family(
    c("mu", "lambda", "alpha"), "lambda > 0 and alpha > 0",
    function(a) a[2] > 0 && a[3] > 0, function(a) c(a[1], Inf),
    function(x, a) log(a[3]) - log(a[2]) - (a[3] + 1) * log1p((x - a[1]) / a[2]),
    function(q, a) 1 - (1 + (q - a[1]) / a[2])^(-a[3]),
    function(n, a) a[1] + a[2] * expm1(stats::rexp(n) / a[3])),
  uniform = .family(
    c("alpha", "beta"), "alpha < beta", function(a) a[1] < a[2],
    function(a) c(a[1], a[2]),
    function(x, a) stats::dunif(x, a[1], a[2], log = TRUE),
    function(q, a) stats::punif(q, a[1], a[2]),
    function(n, a) stats::runif(n, a[1], a[2])),
  beta = .family(
    c("alpha", "beta"), "alpha > 0 and beta > 0",
    function(a) a[1] > 0 && a[2] > 0, function(a) c(0, 1),
    function(x, a) stats::dbeta(x, a[1], a[2], log = TRUE),
    function(q, a) stats::pbeta(q, a[1], a[2]),
    function(n, a) stats::rbeta(n, a[1], a[2])),
  lkj = .family(
    "eta", "eta > 0", function(a) a[1] > 0, function(a) c(-1, 1),
    function(x, a) stats::dbeta((x + 1) / 2, a[1], a[1], log = TRUE) - log(2),
    function(q, a) stats::pbeta((q + 1) / 2, a[1], a[1]),
    function(n, a) 2 * stats::rbeta(n, a[1], a[1]) - 1)
)

# This function reads a prior written as `family(arguments)`. It returns the
# family name and the arguments as numbers (NA where one is not a number), or
# NULL when the text does not have that form.
.read_prior <- function(x) {
  parts <- regmatches(x, regexec("^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\((.*)\\)\\s*$",
                                 x))[[1L]]
  if (!length(parts)) return(NULL)
  args <- if (grepl("^\\s*$", parts[3L])) numeric(0) else
    suppressWarnings(as.numeric(strsplit(parts[3L], ",", fixed = TRUE)[[1L]]))
  list(family = parts[2L], args = args)
}

# This function tells whether a read prior belongs to a known family and has
# valid arguments.
.known_prior <- function(read) {
  family <- .distributions[[read$family %||% ""]]
  !is.null(family) && length(read$args) == length(family$args) &&
    all(is.finite(read$args)) && isTRUE(family$valid(read$args))
}
