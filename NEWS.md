# zzrctsim (development version)

## Documentation

* The design report formerly kept in `analysis/report/report.Rmd` is
  now part of the package documentation, so it is installed, checked
  and rebuilt with the package rather than rendered separately.

  * New `vignette('design')`: purpose, the layered pipeline and its
    three contracts, what is deliberately out of scope, the Monte
    Carlo conventions of the driver, the exponential-family scope,
    supported families, the three retained scales, which operations
    are Gaussian-only and how they are guarded, and open work.
  * New `vignette('trial-design')`: run-in, treatment and common-close
    phases, the common-close and per-arm hinge conventions, staggered
    accrual with a common close-out, and sample size by simulation
    with and without 30% MAR dropout. The sample-size result is
    precomputed in `inst/extdata/sample-size-dropout.csv` by
    `inst/scripts/sample-size-dropout.R`, because the computation
    takes over a minute.
  * `vignette('missing-data')` gains a section on non-Gaussian
    responses, showing which dropout and reference-based operations
    are refused and why; `vignette('dual-parameterization')` states
    the reachability condition explicitly and its Gaussian-only
    scope; `vignette('time-to-event')` documents the recurrent-event
    and delayed-effect inversions.
* `here` and `pkgload`, used only by the report, are removed from
  Suggests.

# zzrctsim 0.2.0

## Correctness fixes

* **The Monte Carlo standard error of `mean_ci_width` was understated
  when any confidence interval was missing.** The estimate averaged the
  non-missing widths, correctly, but its MCSE divided by the number of
  converged replicates rather than by the number of widths actually
  averaged. With 100 of 500 intervals absent the reported MCSE was 11
  percent below the truth. That matters more than its size suggests:
  the MCSE is what says whether enough replicates have been run, so
  understating it makes a simulation look more precise than its
  replicates support. Every other measure already counted its own
  non-missing values.

## Tests

* `test_morris.R` now checks each performance measure against the
  Morris, White and Crowther (2019) Table 6 formula computed by hand
  in the test, rather than against the implementation's own output:
  bias, empirical SE, model SE, relative error of the model SE, MSE,
  coverage, bias-eliminated coverage, rejection rate and mean interval
  width, plus the MCSEs of bias and empirical SE.
* Suite grows from 311 assertions to 338.

## Verified unchanged

Checked during review and found correct:

* `sim_streams()` and `with_rng_state()` restore both the generator
  kind and the caller's `.Random.seed`, so neither disturbs the calling
  session. Replaying a substream reproduces its draws exactly, distinct
  substreams differ, no two substream states coincide, and across 400
  substreams the pairwise correlations are consistent with
  independence (mean 0.0005) with first draws passing a
  Kolmogorov-Smirnov test against the standard normal.
* `bias_elim_coverage` is written as a shifted comparison but is
  algebraically the coverage of `mean(est)`, which is its definition.
* `mean_ci_width` reports the full interval width, matching its
  documentation, despite being computed via a half-width.

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
