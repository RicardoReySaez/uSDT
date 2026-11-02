# Rebuild the validation of the power approximation shown in the power
# vignette.
#
# For eight designs of two populations, the power of the three point nulls is
# computed three ways: the closed form, the draws of the sample mean and
# covariance that usdt_power() uses by default, and the simulation of whole
# studies analysed with hsdt(). The simulation takes about ten minutes on 14
# cores. Run it from the package root after changing either engine.

pkgload::load_all(".", quiet = TRUE)

ncores <- min(14L, parallel::detectCores() - 2L)
nsim <- 1000L
hypotheses <- c("diff = 0", "rho = 0", "intercept = 0")

# A plausible unconscious-processing population, and the estimates of
# hsdt() on the data of Vadillo et al. (2025), whose awareness task is
# unreliable.
populations <- list(
  A = usdt_population(dprime = c(awareness = 0.15, influence = 0.30),
                      sd_dprime = c(awareness = 0.35, influence = 0.30),
                      cor_dprime = 0.3, criterion = 0.1, sd_criterion = 0.3),
  Vadillo = usdt_population(dprime = c(awareness = 0.235, influence = 0.128),
                            sd_dprime = c(awareness = 0.122, influence = 0.088),
                            cor_dprime = 0.49, criterion = 0.036,
                            sd_criterion = 0.25))
designs <- data.frame(
  population = c(rep("A", 5L), rep("Vadillo", 3L)),
  subjects = c(40, 80, 160, 80, 80, 104, 300, 104),
  trials_awareness = c(64, 64, 64, 32, 128, 64, 64, 256),
  trials_influence = c(200, 200, 200, 100, 400, 384, 384, 384))

split <- c(awareness = FALSE, influence = TRUE)
tests <- lapply(hypotheses, .parse_hypothesis)
base <- .moment_base()

rows <- lapply(seq_len(nrow(designs)), function(i) {
  d <- designs[i, ]
  pop <- populations[[d$population]]
  scenario <- as.list(as.data.frame(pop)[1L, ])
  cells <- list(signal = c(awareness = d$trials_awareness / 2,
                           influence = d$trials_influence / 2),
                noise = c(awareness = d$trials_awareness / 2,
                          influence = d$trials_influence / 2))
  v <- .error_variance(scenario, cells, split)
  approximate <- .power_approximate(scenario, d$subjects, v, tests, 0.05,
                                    base)$rows
  simulated <- usdt_power(pop, subjects = d$subjects,
                          trials = c(awareness = d$trials_awareness,
                                     influence = d$trials_influence),
                          hypothesis = hypotheses, dichotomize = "influence",
                          method = "simulation", nsim = nsim, seed = 2026L + i,
                          ncores = ncores, progress = FALSE)
  message(d$population, " ", d$subjects, " subjects: done")
  data.frame(d, hypothesis = hypotheses,
             closed = approximate$power_closed,
             approximate = approximate$power,
             simulation = simulated$power,
             simulation_low = simulated$power.low,
             simulation_high = simulated$power.high,
             estimable_approximate = approximate$estimable,
             estimable_simulation = simulated$estimable,
             nsim = nsim, row.names = NULL)
})
validation <- do.call(rbind, rows)

dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(validation, file.path("inst", "extdata",
                                       "power-validation.csv"),
                 row.names = FALSE)
