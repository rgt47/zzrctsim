# zzrctsim: User Interface, Help System, Coding Practices, and CRAN Readiness Evaluation

*2026-08-14 11:24 PDT*

Scope: the exported API surface (39 functions), the roxygen help
system (37 Rd topics), the `R/` sources and `inst/tinytest/` suite,
and a full `R CMD check --as-cran` run (R 4.6.1, aarch64-apple-darwin,
`_R_CHECK_FORCE_SUGGESTS_=false` because `simr` is not installed on
this host). Epistemic status is marked throughout: **verified** means
code was executed and output observed; **inspected** means source was
read directly.

## Verdict

**Not ready for CRAN submission.** `R CMD check --as-cran` returns
2 ERRORs and 4 NOTEs (verified). Beyond the check, this evaluation
found two functional bugs a user would hit in documented workflows,
one statistical-correctness defect in the dropout mechanism labeling,
and a set of help-system gaps concentrated at the package's front
door. The underlying architecture, RNG discipline, and the
statistical quality of the existing test suite are genuinely strong;
the distance to a submittable package is mostly packaging metadata,
one broken example, one test-suite portability bug, and documentation
completion, not redesign.

## 1. R CMD check --as-cran results (verified)

Status: **2 ERRORs, 0 WARNINGs, 4 NOTEs.**

**ERROR 1 — broken shipped example.** `apply_dropout()`'s `@examples`
(`R/generate.R:238-247`) builds a schedule with no run-in
(`trial_schedule(treatment = 5, interval = 3)`), so `runin_design()`
emits no `x_hinge` column, yet the example's `beta` includes
`x_hinge = -0.1`. `generate_outcomes()` correctly rejects it:
"`beta` names not found in `design`: x_hinge". Fix: add
`run_in = 1` (or drop `x_hinge` from the example's beta).

**ERROR 2 — test suite fails against the installed package.**
`test_dgm.R` calls the internal helper `.zgz()` without a
`zzrctsim:::` prefix. Under `pkgload::load_all()` internals are
visible, so the local 239-test run passes; under `R CMD check` the
tests run against the installed namespace and die with "could not
find function \".zgz\"". This is why local testing and check
disagree. Fix: `zzrctsim:::.zgz(...)` (and audit all test files for
other unqualified internal calls).

**NOTEs.**

1. *CRAN incoming feasibility*: new submission; version 0.0.0.9000
   has "large components" — bump to a release version (0.1.0) before
   submitting.
2. *Hidden files*: `.zzcollab`, `.zzcollab-state`, `.devcontainer`
   ship in the tarball. Add to `.Rbuildignore`.
3. *Top-level files*: `LICENSE` is not referenced by the DESCRIPTION
   `License: GPL-3` field — and its content is an MIT-template stub
   ("YEAR / COPYRIGHT HOLDER"), which is wrong for GPL-3 (inspected).
   Delete the file or switch to `GPL-3 + file LICENSE` with correct
   content (deleting is the standard choice). `CITATION.cff` and
   `zzcollab.yaml` are non-standard top-level files — `.Rbuildignore`
   both.
4. *CITATION in a non-standard place*: CRAN expects `inst/CITATION`.
   The current `CITATION.cff` is also stale (verified: version 0.1.0
   vs DESCRIPTION 0.0.0.9000, "Your Institution" placeholder, empty
   ORCID, boilerplate zzcollab abstract).

Additional packaging items check did not flag on this host:
`tools/` contains the zzcollab render wrapper (`render.sh`,
`stamp-render.R`, `stamp.tex`, `README.md`), which is unrelated to
the package build and ships in the tarball — `.Rbuildignore` it.
`R/dgm.R:243` calls `utils::head()` but `utils` is not in
DESCRIPTION Imports (inspected; not flagged by this check run, but
current R-devel checks do flag undeclared `utils`). DESCRIPTION has
no `URL`/`BugReports` fields; now that
`https://github.com/rgt47/zzrctsim` exists, add both.

## 2. Functional bugs found during the audit

These are not check failures; they are defects a user hits on
documented paths.

**2.1 `generate_outcomes()` drops its attribute contract on two of
three paths, breaking `reference_based()` (verified).** The
`beta`/`intercept`/`family` attributes are attached only on the plain
Gaussian path (`R/generate.R:69-72`); the GLMM and mixture branches
return early without them. `reference_based()` on Gaussian mixture
data then errors with "`dat` must carry a `beta` attribute; use
generate_outcomes()" even though the data came from
`generate_outcomes()` exactly as the message instructs. The
documented return ("`design` with two columns added: `mu` and `y`")
is also wrong for those paths (GLMM adds `eta`, `cmean`, a `ranef`
attribute; mixture adds `component`).

**2.2 MAR mislabeling at a first eligible visit (inspected;
statistical correctness).** In both `apply_dropout()`
(`R/generate.R:176`) and `dropout_mask()` (`R/missing.R:173`), the
"previous response" for a subject whose first visit is already
eligible falls back to the *current* response: the hazard then
conditions on an unobserved value under `psi1`, making the mechanism
MNAR at that visit while `dropout_mechanism()` and the `spec`
attribute still report MAR. Fix: contribute nothing (or `NA`
handling) at a first eligible visit, or refuse the configuration.

**2.3 `sample_size()` confirmation is not independent of selection
(inspected).** The confirmation `sim_power()` call reuses the same
`seed` as the curve run (`R/power.R:101-146`), so the "independent
confirmation" the documentation tells the user to quote shares its
randomness with the run that selected `n`. Fix: derive a distinct
confirmation seed, or document the common-random-numbers sharing as
deliberate.

**2.4 Missing-input validation with obscure failure modes
(verified).** `run_simulation()` with a partially named `analyse`
list dies with "arguments imply differing number of rows: 1, 0";
`fit_ancova()` on data lacking a `time` column dies inside `max()`;
`dgm_marginal()` accepts a non-PSD `V` silently and fails B
replicates deep inside `MASS::mvrnorm`; `reference_based()` with a
non-covering mask dies with "missing value where TRUE/FALSE needed"
(the analogous check exists in `apply_mask()` and is simply missing
here). Each has a one-line fix at the entry point.

**2.5 Silent-surprise behaviors (verified/inspected).**
`cov_ar1(sd = c(1, 2))` evaluated on a length-3 grid recycles to
variances 1, 4, 1 with no warning — under staggered accrual with
subject-specific grids this is almost never intended.
`trial_schedule(treatment = 2.5)` silently truncates to 2.
`fit_ancova()` in a 3+ arm design silently reports only the first
non-reference contrast. `accrue(pattern = "poisson")` violates its
own documented `[0, period]` support (observed max 10.89 for
`period = 10`). `apply_dropout()` records `target = 0.25` in its spec
attribute even when calibration was bypassed via `psi0` (the newer
`dropout_mask()` gets this right).

## 3. Help system: extent and accuracy

**Coverage (verified).** All 39 exports are documented across 37 Rd
topics; every `@param` block matches the actual formals — no stale or
missing parameter names anywhere. That baseline is better than most
first submissions.

**The front door is missing.** There is no package-level topic
(`?zzrctsim` resolves to nothing; no `_PACKAGE` roxygen sentinel) and
no root `README.md` — the GitHub landing page published yesterday is
a bare file listing. For "ready access to help when most likely
needed," these are the two highest-leverage additions: a
`zzrctsim-package.R` topic linking the six vignettes and the main
workflow functions, and a README with the getting-started example.

**Examples are sparse: 7 of 39 exports (18%) have `@examples`, and
one of the seven is the broken one in ERROR 1.** CRAN does not
mandate examples on every export, but the main workflow entry points
(`run_simulation()`, `sample_size()`, `dgm_marginal()`,
`trial_schedule()`) deserve them.

**Return-value documentation is incomplete on 15 of 39 exports
(inspected).** The worst pattern is class-only returns: the
`dgm_conditional()`/`dgm_marginal()` docs name the returned class but
describe none of its 8 elements; same for `dgm_tte()` (11 elements),
`dgm_mixture()`, `estimand()`, `fit_result()`. Specific inaccuracies:
`compute_performance()` enumerates 9 measures but emits 11
(`conv_rate` and `mean_ci_width` are undocumented — a user filtering
on the documented list silently misses both);
`run_simulation()`'s `@return` omits the `errors` attribute;
`runin_design()`'s `@param arm` describes an obsolete two-arm 0/1
interface when the implementation supports arbitrary multi-level
factors, and its `@return` omits the `J0 > 0` condition on `x_hinge`
— the precise omission that produced ERROR 1; `generate_tte()`
attaches an undocumented `lp` attribute; `certify()` says
"positive-definite" where the code accepts PSD.

**Vignettes (verified: all six render cleanly).** Coverage of the
main workflows is good. Orphan functionality with no example and no
vignette: the staggered-accrual realization pair
`close_out()`/`realize_schedule()` — a headline DESCRIPTION feature
never demonstrated end-to-end — plus `sim_power()`/`power_curve()`
(reached only through `sample_size()`), and the RNG-replay trio
`sim_streams()`/`with_rng_state()`/`linpred_cov()`. The mixture DGM
has an example but no vignette treatment.

## 4. User interface / API design

**Strengths (inspected).** Naming is uniformly `snake_case` with a
clean noun/verb split: constructors are nouns (`trial_schedule`,
`dgm_*`, `cov_*`, `estimand`), actions are verbs (`generate_*`,
`apply_*`, `certify`, `reconstruct`, `accrue`). Twelve S3 print
methods give informative at-the-console summaries (the
`<sample_size>` print showing target, selected n, and confirmed
power with MCSE is exemplary). The mask-generation/mask-application
separation (`dropout_mask()` / `apply_mask()`) is a genuinely good
design seam. Where error messages were hand-written
(`.check_family_for_hazard`, `cov_at.dgm_marginal`,
`as_conditional`) they name the argument and the remedy.

**Weaknesses.** Argument order is inconsistent between levels of the
same stack: `run_simulation(B, generate, analyse, estimand, ...)`
versus `sim_power(..., estimand, analyse, ...)`. Several exported
entry points fail via bare `stopifnot()` with expression-echo
messages ("nrow(V) == ncol(V) is not TRUE") rather than named
messages — the named-expression `stopifnot("message" = cond)` form
would fix these in place. The `analyse` argument name is British
spelling, inconsistent with the package's otherwise US prose;
renaming is API-breaking, so if it is to change, before first
release is the only cheap moment. Two near-duplicate generations of
the same facility coexist (`apply_dropout()` versus
`dropout_mask()`+`apply_mask()`; the older one carries the spec-
attribute bug in 2.5) — deprecating the older path before first
release would shrink the surface users must learn.

## 5. Coding practices

**Verified sound (run + inspected).** RNG discipline is implemented
as advertised: `sim_streams()` pins L'Ecuyer-CMRG, walks
`parallel::nextRNGStream`, and restores both `RNGkind` and
`.Random.seed` via `on.exit`; `test_morris.R` genuinely tests
restoration, stream replay, determinism, and CRN pairing. No
`sapply`, no `T`/`F`, no `1:n`, no `library()` in `R/`, no
`set.seed()` outside the documented RNG layer, `vapply` with type
templates, `drop = FALSE` discipline throughout, all 12 S3 methods
properly registered, dispatch via `inherits()`. Dependency hygiene is
clean except the `utils` omission. The existing test suite is
behavioral, not structural: type I error and coverage against MCSE
bounds, Weibull shape recovery via Cox fits, `cox.zph` proportional-
hazards violation under `effect_lag`, mixture marginal covariance
against the closed form.

**Gaps.** `R/power.R` — `sim_power()`, `power_curve()`,
`sample_size()` — has zero test coverage, and `sample_size()` has
four distinct branches including two stop paths and a warning path.
These are the top-of-stack functions users call first.
`inst/tinytest/test-basic.R` is a placeholder (`expect_true(TRUE)`)
whose hyphenated name also breaks the suite's `test_*.R` convention;
delete it. Hot loops scan `which(dat$id == i)` per subject inside
per-replicate code paths (O(subjects x rows) in seven functions) and
`generate_tte()` grows its recurrent-event vector by concatenation —
not correctness issues, but the obvious first performance pass.

## 6. Prioritized checklist

### (a) CRAN blockers

1. Fix the `apply_dropout()` example (add `run_in = 1` or drop
   `x_hinge`).
2. Qualify `.zgz()` as `zzrctsim:::.zgz()` in `test_dgm.R`; audit all
   test files for other unqualified internals.
3. Bump Version to 0.1.0.
4. Delete the stub `LICENSE` file (GPL-3 needs none) or make it a
   correctly referenced `+ file LICENSE`.
5. `.Rbuildignore`: `.zzcollab`, `.zzcollab-state`, `.devcontainer`,
   `CITATION.cff`, `zzcollab.yaml`, `tools/`.
6. Add `utils` to Imports (or remove the `utils::head()` call).
7. Move citation metadata to `inst/CITATION` and fix its stale
   content (version, institution, ORCID).
8. Add `URL` and `BugReports` to DESCRIPTION.
9. Re-run `R CMD check --as-cran` with `simr` installed, then on
   win-builder/R-devel before submitting (not done here).

### (b) Bugs to fix before anyone depends on the behavior

10. Attach `beta`/`intercept`/`family` attributes on the GLMM and
    mixture paths of `generate_outcomes()` (unbreaks
    `reference_based()` on those paths).
11. Fix or refuse the first-eligible-visit "previous response"
    fallback so the reported MAR label is true.
12. Decouple `sample_size()`'s confirmation seed from the curve seed,
    or document the sharing.
13. Entry-point validation: fully named `analyse` in
    `run_simulation()`; required-columns check in `fit_ancova()`;
    PSD check in `dgm_marginal()`; mask-coverage check in
    `reference_based()`; no silent `sd` recycling; no silent
    fractional truncation in `trial_schedule()`; warn on dropped
    contrasts in multi-arm `fit_ancova()`; fix `accrue()` "poisson"
    doc or behavior; fix `apply_dropout()`'s spec `target`.

### (c) Help-system completion

14. Add a `_PACKAGE` topic and a root README.
15. Complete the 15 incomplete `@return` sections; correct the
    `compute_performance()` measure list, `runin_design()`'s `arm`
    doc and `x_hinge` condition, `certify()`'s PD/PSD wording.
16. Add `@examples` to the main workflow entry points; add tests for
    `R/power.R`; delete `test-basic.R`.
17. Demonstrate `close_out()`/`realize_schedule()` and the RNG-replay
    workflow in a vignette or examples.

### (d) Design decisions best made before first release

18. Decide `analyse` vs `analyze`; harmonize argument order between
    `run_simulation()` and `sim_power()`; deprecate `apply_dropout()`
    in favor of `dropout_mask()`+`apply_mask()`.

## Not evaluated

Windows/Linux checks, win-builder, R-devel, a check with `simr`
present, spell-check of Rd files, and rhub. All should be run after
the (a)-list lands. External citation DOIs in DESCRIPTION were not
resolved.
