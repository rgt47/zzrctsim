# zzrctsim 0.1.0

First release.

`zzrctsim` simulates longitudinal randomized clinical trials in order
to evaluate design and analysis choices before a trial is run. A study
is expressed as five separable steps, one per element of ADEMP
(Morris, White and Crowther 2019), so that a manuscript's Methods
subsections correspond to code rather than having to be reconstructed
from a script. See `?zzrctsim` for the map and
`vignette("getting-started")` for the walkthrough.

## Data-generating mechanisms

* A DGM may be declared **conditionally**, through random effects
  (`dgm_conditional()`), or **marginally**, through a covariance
  matrix or a covariance function of time (`dgm_marginal()`,
  `cov_ar1()`, `cov_cs()`).
* `as_marginal()` and `as_conditional()` convert between the two
  parameterizations where an exact conversion exists. Where none
  does, `certify()` reports the minimum number of random effects that
  would reach the target covariance, and `reach_error()` quantifies
  the approximation error at a given budget.
* Non-Gaussian outcomes (binomial, Poisson, negative binomial) are
  generated through a GLMM construction; `dgm_tte()` covers
  single and recurrent time-to-event endpoints, including
  piecewise-exponential baselines, delayed effects and frailty;
  `dgm_mixture()` covers latent-subgroup trajectories.

## Design

* `trial_schedule()` spans run-in, treatment and common-close phases;
  `runin_design()` expands it over subjects and builds the model
  columns.
* Staggered accrual is supported end to end: `accrue()` draws
  enrollment times, `close_out()` finds the common close-out date, and
  `realize_schedule()` gives each subject the grid their enrollment
  allows. The resulting ragged design feeds `runin_design()` directly.

## Missingness

* Generated from the Diggle and Kenward (1994) selection model,
  calibrated in closed form to a target dropout proportion, with MCAR,
  MAR and MNAR as nested restrictions of one hazard.
* Mask generation (`dropout_mask()`) and mask application
  (`apply_mask()`) are deliberately separate, so a pattern may be
  drawn, imported from a completed trial (`mask_from_data()`), or
  constructed by hand, and held fixed across replicates under common
  random numbers.
* `reference_based()` supplies post-discontinuation trajectories
  (J2R, CIR, CR, LMCF).

## Simulation and reporting

* `run_simulation()` drives the replicate loop against a declared
  `estimand()`, pinning the L'Ecuyer-CMRG generator and giving each
  replicate its own substream, with the state captured so that any
  single replicate can be replayed in isolation via
  `with_rng_state()`. The calling session's RNG kind and seed are
  restored on exit.
* `compute_performance()` reports the full performance suite -- bias,
  empirical and model-based standard error, relative error, mean
  squared error, coverage, bias-eliminated coverage and rejection
  rate -- each with its Monte Carlo standard error.
  `nsim_for_mcse()` inverts the relationship to choose `B` from a
  precision requirement.
* `sim_power()`, `power_curve()` and `sample_size()` answer the power
  question by simulation. The selected sample size is confirmed by an
  independent run on a separate RNG stream, and it is that confirmed
  power, with its Monte Carlo standard error, that should be quoted.

## Known limitations

* Non-Gaussian families are generated but only one reference fitter
  (`fit_ancova()`) ships; supply your own via the `fit_result()`
  contract.
* `fit_ancova()` reports a single contrast, so call it once per
  contrast in a trial with three or more arms.
* No small-sample denominator degrees-of-freedom adjustment
  (Kenward-Roger, Satterthwaite) is applied.
* `accrue(pattern = "poisson")` sets the *expected* accrual span, so
  a realization may overrun `period`.
