# zzrctsim

Simulation of longitudinal randomized trials, for evaluating design
and analysis choices before a trial is run.

The data-generating mechanism may be declared **either** conditionally,
through random effects, **or** marginally, through a covariance matrix
or a covariance function of time, with exact conversion between the two
parameterizations where one exists and a computable rank certificate
where it does not. Visit schedules span run-in, treatment and
common-close phases, with staggered accrual and a common close-out
date. Missingness is generated from the Diggle and Kenward (1994)
selection model, calibrated to a target dropout proportion, with mask
generation and mask application kept separate so that empirical
patterns can be imported.

Simulations follow Morris, White and Crowther (2019): L'Ecuyer-CMRG
substreams with per-replicate state capture, and performance measures
reported with Monte Carlo standard errors throughout.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rgt47/zzrctsim")
```

Imports are deliberately few: `MASS`, `parallel`, `stats`, `utils`.

## The five steps

A study is always the same five steps, kept separate so that each
corresponds to one element of ADEMP.

```r
library(zzrctsim)

# 1. Design: when are subjects seen?
s <- trial_schedule(treatment = 4, interval = 3)
arm <- factor(rep(c("placebo", "active"), each = 100),
              levels = c("placebo", "active"))
d <- runin_design(s, arm, reference = "placebo")

# The design carries the model columns your `beta` will name. Look at
# them before writing step 3 -- `x_slope` is the common slope and
# `x_trt_active` the treatment effect on the rate of change.
names(d)
#> "id" "arm" "index" "time" "phase" "x_slope" "x_trt_active"

# 2. Data-generating mechanism: a random intercept (variance 9) and a
#    random slope (variance 0.04), in the order the default
#    `z = function(t) cbind(1, t)` gives them, plus residual variance 4.
g <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)

# 3. Generation: draw one trial.
b <- c(x_slope = 0.5, x_trt_active = -0.25)
dat <- generate_outcomes(d, g, beta = b, intercept = 20)

# 4. Analysis: anything returning a fit_result().
fit_ancova(dat, reference = "placebo")

# 5. Replicates and performance, with Monte Carlo standard errors.
theta <- b[["x_trt_active"]] * max(s$time)
res <- run_simulation(
  B = 200,
  generate = function() generate_outcomes(d, g, beta = b,
                                          intercept = 20),
  analyze = list(ancova = function(z)
    fit_ancova(z, reference = "placebo")),
  estimand = estimand("final-visit contrast", theta),
  seed = 42L
)
compute_performance(res)
```

`compute_performance()` reports bias, empirical and model-based
standard errors, coverage, mean squared error and rejection rate, each
with its Monte Carlo standard error, so that an apparent bias can be
told apart from noise from too few replicates.

## Sample size

The power question can be inverted directly. The selected size is
confirmed by an independent run on a different RNG stream, and it is
that confirmed power, with its Monte Carlo standard error, that should
be quoted.

```r
design_fn <- function(n) {
  runin_design(s, factor(rep(c("placebo", "active"), each = n),
                         levels = c("placebo", "active")),
               reference = "placebo")
}
sample_size(
  target = 0.80, n_grid = c(10, 15, 20, 30),
  design_fn = design_fn, dgm = g, beta = b, intercept = 20,
  estimand = estimand("final contrast", theta),
  analyze = list(ancova = function(z)
    fit_ancova(z, reference = "placebo")),
  B = 300L, seed = 5L
)
```

## Documentation

`?zzrctsim` gives the orientation above with links to every stage.

- `vignette("getting-started")` — the five steps end to end
- `vignette("dual-parameterization")` — conditional versus marginal
  DGMs, and when the two are interchangeable
- `vignette("missing-data")` — dropout mechanisms, reusable masks,
  reference-based trajectories
- `vignette("time-to-event")` — single and recurrent events
- `vignette("comparison-lme4")`, `vignette("comparison-simr")` — how
  this relates to tools you may already use

## Scope and limitations

- The response is continuous and Gaussian on the main path.
  Non-Gaussian families are generated through a GLMM construction, but
  only one reference fitter ships.
- `fit_ancova()` reports a single contrast; call it once per contrast
  in a trial with three or more arms.
- No small-sample denominator degrees-of-freedom adjustment
  (Kenward-Roger, Satterthwaite) is applied.
- `accrue(pattern = "poisson")` sets the *expected* accrual span, so a
  realization may overrun `period`.

## Testing

```r
tinytest::test_package("zzrctsim")
```

## References

Diggle, P. and Kenward, M. G. (1994). Informative drop-out in
longitudinal data analysis. *Journal of the Royal Statistical Society:
Series C* 43(1), 49-93.

Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using
simulation studies to evaluate statistical methods. *Statistics in
Medicine* 38(11), 2074-2102.

## License

GPL-3.
