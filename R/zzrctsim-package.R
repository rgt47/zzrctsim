#' @keywords internal
#' @aliases zzrctsim-package
#' @details
#' A simulation study in `zzrctsim` is five steps, kept separate so
#' that each maps onto one element of ADEMP (Morris, White and
#' Crowther 2019) and a manuscript's Methods subsections correspond to
#' code rather than having to be reconstructed from a script.
#'
#' \describe{
#'   \item{1. Design}{When are subjects seen? [trial_schedule()]
#'     defines the visit grid across run-in, treatment and
#'     common-close phases; [runin_design()] expands it over subjects
#'     and builds the model columns. For staggered entry,
#'     [accrue()] draws enrollment times, [close_out()] finds the
#'     common close-out date, and [realize_schedule()] cuts each
#'     subject's grid at that date.}
#'   \item{2. Data-generating mechanism}{What produces the response?
#'     Declare it conditionally through random effects with
#'     [dgm_conditional()], or marginally through a covariance with
#'     [dgm_marginal()] and the [cov_ar1()] / [cov_cs()] builders.
#'     [as_marginal()] and [as_conditional()] convert between the two
#'     where an exact conversion exists; [certify()] reports the rank
#'     needed when it does not. [dgm_tte()] and [dgm_mixture()] cover
#'     time-to-event and latent-subgroup outcomes.}
#'   \item{3. Generation}{Draw one trial with [generate_outcomes()],
#'     or [generate_tte()] for event times. Complete data are always
#'     generated first; missingness is a separate layer, applied with
#'     [dropout_mask()] and [apply_mask()], which is what allows the
#'     cost of dropout to be measured against the complete-data run
#'     and what makes response-dependent (MAR, MNAR) mechanisms
#'     expressible at all. [reference_based()] supplies
#'     post-discontinuation trajectories.}
#'   \item{4. Analysis}{Fit anything, provided it returns a
#'     [fit_result()]. [fit_ancova()] is supplied as a reference
#'     method.}
#'   \item{5. Performance}{[run_simulation()] drives the replicate
#'     loop against a declared [estimand()];
#'     [compute_performance()] summarizes across replicates with a
#'     Monte Carlo standard error on every measure.
#'     [nsim_for_mcse()] chooses `B` from a precision requirement,
#'     and [sim_power()], [power_curve()] and [sample_size()] invert
#'     the power question.}
#' }
#'
#' @section Reproducibility:
#' [run_simulation()] pins the L'Ecuyer-CMRG generator and gives each
#' replicate its own substream via [sim_streams()], capturing the
#' state so that any single replicate can be re-executed in isolation
#' with [with_rng_state()]. The RNG kind and seed of the calling
#' session are restored on exit.
#'
#' @section Where to start:
#' `vignette("getting-started", "zzrctsim")` walks through the five
#' steps end to end. See also
#' `vignette("dual-parameterization", "zzrctsim")`,
#' `vignette("missing-data", "zzrctsim")`,
#' `vignette("time-to-event", "zzrctsim")`, and the
#' `comparison-lme4` and `comparison-simr` vignettes for how this
#' relates to tools you may already use.
#'
#' @references
#' Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using
#' simulation studies to evaluate statistical methods. *Statistics in
#' Medicine* 38(11), 2074-2102. \doi{10.1002/sim.8086}
#'
#' Diggle, P. and Kenward, M. G. (1994). Informative drop-out in
#' longitudinal data analysis. *Journal of the Royal Statistical
#' Society: Series C* 43(1), 49-93. \doi{10.2307/2986113}
"_PACKAGE"
