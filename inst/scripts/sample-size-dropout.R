# Regenerates inst/extdata/sample-size-dropout.csv, the precomputed
# result shown in vignette('trial-design'). The full computation takes
# over a minute, which is too long for a vignette built on CRAN, so the
# vignette reads this file and shows the same code with eval = FALSE.
#
# Run from the package root with the package installed:
#   Rscript inst/scripts/sample-size-dropout.R

library(zzrctsim)

sch <- trial_schedule(treatment = 4, interval = 3)
gd <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
b <- c(x_slope = 0.5, x_trt_active = -0.12)
design_fn <- function(n) {
  runin_design(sch,
               factor(rep(c('placebo', 'active'), each = n),
                      levels = c('placebo', 'active')),
               reference = 'placebo')
}
common <- list(
  design_fn = design_fn, dgm = gd, beta = b, intercept = 20,
  estimand = estimand('final contrast',
                      b[['x_trt_active']] * max(sch$time)),
  analyze = list(ancova = function(z) {
    fit_ancova(z, reference = 'placebo')
  }),
  B = 300L, seed = 5L
)
ss_complete <- do.call(sample_size, c(list(
  target = 0.80, n_grid = c(40, 70, 110, 160), confirm_B = 800L
), common))
ss_dropout <- do.call(sample_size, c(list(
  target = 0.80, n_grid = c(60, 110, 160, 220), confirm_B = 800L,
  dropout = list(target = 0.30, psi1 = 0.10)
), common))

out <- data.frame(
  scenario = c('complete data', '30% MAR dropout'),
  n_per_arm = c(ss_complete$n_per_arm, ss_dropout$n_per_arm),
  confirmed_power = c(ss_complete$confirmation$power,
                      ss_dropout$confirmation$power),
  mcse = c(ss_complete$confirmation$mcse, ss_dropout$confirmation$mcse)
)
stopifnot(nrow(out) == 2L, !anyNA(out))
utils::write.csv(out, file.path('inst', 'extdata',
                                'sample-size-dropout.csv'),
                 row.names = FALSE)
