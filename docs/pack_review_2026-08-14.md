# Referee Review: the `zzrctsim` R Package

*2026-08-14 13:10 PDT*

Scope: `R CMD check --as-cran`, four automated tools, the 38-export
help system, the user interface judged from a first-use walkthrough,
and the `R/` sources and tinytest suite. This is an update to
`docs/ui-help-cran-evaluation-2026-08-14.md` (same day, earlier);
its revision history is carried forward in Section 9.

Epistemic marking is used throughout: **verified** means code was
executed and output observed in this pass; **inspected** means source
was read directly. Claims sourced from a subagent are marked with
which of the two that agent reported, and the load-bearing ones were
re-run independently by the reviewer.

## 1. Verdict

**CRAN-ready on the check; not ready as a user-facing release.**

`R CMD check --as-cran`: **0 errors, 0 warnings, 1 NOTE** (verified;
1m 37s; the NOTE is "New submission," unavoidable for a first
submission). Test coverage **84.86%** (verified, `covr`). All six
vignettes render and all 14 example blocks execute (verified).

That is the entire good news about *readiness*, and it is why the
distinction this review's own ground rules insist on matters here:
**passing the check is not being correct.** This pass found a
documented workflow that cannot be executed, a sample-size function
that returns the wrong answer on tied curves, a masking function that
silently fabricates data, and a performance profile that makes the
package's headline use case impractical at realistic sample sizes.
None of these is visible to `R CMD check`.

The three prior remediation rounds today (the (a)–(d) lists of the
earlier review) were real and hold up: every one of their fixes that
this pass re-tested independently is confirmed. The defects below are
newly found, not regressions, with two exceptions noted in Section 9.

## 2. `R CMD check` and tooling results

### 2.1 `R CMD check --as-cran` (verified)

0 errors, 0 warnings, 1 NOTE. The NOTE is the maintainer/new-submission
notice. Both prior runs today with and without `_R_CHECK_FORCE_SUGGESTS_`
agreed; `simr`, `lme4` and `survival` are all installed here, so no
vignette or test was silently skipped in this run.

Manual CRAN-hygiene checks (verified):

- **`NEWS.md` is absent.** Expected for a release; a first-version
  entry should be added.
- **DESCRIPTION lacks a `Language` field.** `spelling` defaults to
  `en-US` and says so; add `Language: en-US` explicitly.
- Tarball contents are clean: `build`, `DESCRIPTION`, `inst`, `man`,
  `NAMESPACE`, `R`, `README.md`, `tests`, `vignettes`. No zzcollab
  scaffolding, no `tools/`, no `CITATION.cff`, no `.devcontainer`.
- Version 0.1.0; `URL`/`BugReports` present; `inst/CITATION` present
  and renders. `.onLoad`/`.onAttach` absent. No `\dontrun{}` or
  `\donttest{}` anywhere. No non-ASCII in `R/`; all files end in a
  newline.
- Every `::` namespace in executable code (`MASS`, `parallel`,
  `stats`, `utils`) is declared in Imports. **A `survival::` match in
  `R/tte.R:214` is inside a roxygen comment, not code** — checked
  specifically because it would otherwise be an undeclared
  dependency; it is a non-issue.

### 2.2 Automated tooling

| Tool | Result |
|---|---|
| `covr::package_coverage()` | **84.86%** (verified) |
| `urlchecker::url_check()` | clean, all URLs correct (verified) |
| `spelling::spell_check_package()` | 2 genuine hits (verified) |
| `lintr::lint_package()` | 414 lints, 69 in `R/` (verified) |
| `codetools::checkUsagePackage()` | 1 genuine finding (verified) |
| `goodpractice::gp()` | **not installed; did not run** |

**Coverage by file** (verified). Weakest: `R/simulate.R` **69.05%** —
which holds `run_simulation()`, `fit_ancova()` and the `fit_result`
contract, the most-used code in the package. Then `R/missing.R`
79.27%, `R/certify.R` 80.00%. Strongest: `R/rng.R` 100%,
`R/performance.R` 97.33%.

**Spelling**: `behaviour` (`vignettes/missing-data.Rmd:156`) and
`centres` (`missing-data.Rmd:68`) are British spellings that survived
today's US-English sweep, which reached `R/` but not the vignette
prose. Every other hit is domain vocabulary (ADEMP, CMRG, Diggle,
Kenward, estimand).

**lintr**: the raw count overstates the problem. 276 of 414 are
`quotes_linter` "Only use double-quotes," all in `inst/tinytest/`,
and they contradict the project's stated single-quote preference —
they are a tool/style disagreement, not defects. `R/` itself carries
only 69, of which 44 are `object_name_linter` firing on the
dot-prefixed internal convention. lintr's `object_usage_linter`
"no visible global function definition" messages are **self-declared
false positives**: it prints "there is no package called 'zzrctsim'
... may lead to false positives" because the package is not
installed. One real hit survives: `widths` assigned and never used.

**codetools**: the "parameter changed by assignment" list is almost
entirely the ordinary `x <- as.matrix(x)` / `rule <- match.arg(rule)`
idiom. One genuine item: `.expected_missing: parameter 'monotone'
may not be used` — deliberate (monotone and intermittent coincide for
"at least one hazard fires") but undocumented at that call site.

## 3. Functional bugs

These are defects a user hits on documented paths. They outrank
everything in Section 2.

### 3.1 The documented staggered-accrual workflow cannot be executed

**Verified** (independently re-run by the reviewer; both subagents
found it separately). `?generate_outcomes` (`R/generate.R:18`) states
the design may come "by `realize_schedule()` passed through
`runin_design()`". That route is impossible:

```
realize_schedule(s, e, cl) |> class()
#> "data.frame"
runin_design(rs, arm)
#> Error: inherits(schedule, "trial_schedule") is not TRUE
```

`realize_schedule()` returns a ragged plain data frame
(`id, enroll, index, time, calendar, phase, h, on_treatment`) with no
`arm` and no `x_*` model columns; `runin_design()` accepts only a
`trial_schedule` and builds only balanced designs
(`R/schedule.R:212,237`).

Staggered accrual with a common close-out is headlined in
`DESCRIPTION:13-14`, in `README.md`, and as stage 1 of the map in
`?zzrctsim`. On the longitudinal path it is unreachable: verified
that `realize_schedule()` output is never fed to `generate_outcomes()`
in any test or vignette, and that `close_out()`/`realize_schedule()`
appear in **zero** vignettes. It is exercised only on the TTE path
(`vignettes/time-to-event.Rmd:60`, via `accrue()` into
`generate_tte()`).

Note the irony: examples for these functions were added earlier today
precisely because they were orphans, and those examples pass — they
stop at `head(rs, 3)` and never reach a DGM, so they demonstrate the
functions without demonstrating the workflow.

**Fix:** either ship a `realized_design()` that attaches `arm` and the
`x_*` columns to ragged output, or delete the claim from
`R/generate.R:18` and `R/zzrctsim-package.R`. The second is honest and
cheap; the first is what the DESCRIPTION promises.

### 3.2 `sample_size()` returns the largest tied size, not the smallest

**Verified** at the exact expression. `?sample_size` promises "the
smallest per-arm sample size attaining `power`" (`R/power.R:141`). On
a tied curve — routine, since power saturates at 1.000 across grid
points —

```
powers 0.6, 0.9, 0.9, 0.9 at n = 50, 100, 150, 200; target 0.9
approx(x = power[o], y = n[o], xout = 0.9, ties = "ordered")$y
#> 200          # 100 already attains it
```

It returns **200 where 100 suffices** — a 2× sample-size overstatement
in the package's flagship deliverable. Separately, on a non-monotone
curve the inversion interpolates *through* the dip silently: powers
`0.70, 0.95, 0.94, 0.99` with target 0.945 gives n = 125, in a region
where power decreases in n.

This is the third distinct defect found in `sample_size()`'s
inversion today; the first two (row/message mismatch in the
upper-bound warning, post-hoc method check) were fixed this morning.
The common cause is that `approx()` assumes a monotone curve that
Monte Carlo error does not deliver.

**Fix:** after interpolating, take
`min(curve$n_per_arm[curve$power >= target])`, and warn when
`is.unsorted(curve$power[order(curve$n_per_arm)])`.

### 3.3 `apply_mask()` silently fabricates an all-`NA` response column

**Verified** (re-run independently):

```
nody <- dat; nody$y <- NULL
apply_mask(nody, mk)      # no error
#> y column CREATED, all NA: TRUE
```

`dat[[response]][m] <- NA_real_` (`R/missing.R:301`) creates the
column when absent. A typo in `response=` therefore yields a fully-`NA`
response that flows downstream and makes `fit_ancova()` return
`null_fit()` — a silent, total loss of the outcome, presenting as
non-convergence rather than as an error.

**Fix:** `if (!response %in% names(dat)) stop(...)` at the top.

### 3.4 `dropout_mask()` is superlinear, making the headline use case impractical

**Verified by timing** (subagent; 8 visits/subject, `psi1 = 0.05`):
0.09 s at n = 100, 0.36 s at n = 400, **2.26 s at n = 1600**.

`.expected_missing()` (`R/missing.R:227`) does `which(dat$id == i)`
per subject, O(n · rows), and `stats::uniroot()` calls it ~40 times
during calibration (`R/missing.R:156`). That is the inner loop of
`sim_power(dropout = ...)`: at n = 1600, B = 1000, roughly **38
minutes of masking alone**, before any model is fitted.

**Fix:** hoist `split(seq_len(nrow(dat)), dat$id)` out of `uniroot`
and precompute the `prev`/`y` vectors so the root-finder only
re-evaluates `plogis(psi0 + const)`.

### 3.5 `run_simulation()` documents warning-trapping it does not do

**Verified** (subagent). `?run_simulation` says "Errors **and
warnings** ... are trapped ... the message is retained in the `errors`
attribute" (`R/simulate.R:206-207`). Only errors are trapped
(`tryCatch(..., error=)`, `R/simulate.R:268,276`). With a three-arm
`fit_ancova` at B = 3, three warnings reach the console and
`attr(res, "errors")` has length 0. At B = 1000 that is 1000 console
warnings, none recorded.

This interacts with a change made earlier today: the multi-arm
warning added to `fit_ancova()` is precisely the warning that now
floods. **Fix:** `withCallingHandlers()` and append to `log$errors`,
or correct the sentence.

### 3.6 `reach_error()` returns silent `NA` and negative variances

**Verified** (re-run independently): `reach_error(V, q = p)` yields
`sigma2 = NA, max_abs_error = NA, rel_frobenius = NA` with no error,
because `(qq+1L):p` counts backwards (`R/certify.R:151`). Unlike its
sibling `certify()`, it validates neither symmetry nor PSD-ness, so a
non-PSD input returns a **negative residual variance** (−1) as if it
were a result.

### 3.7 Smaller silent-surprise defects (all verified by the subagent)

- `as_conditional(dgm, times, q = 5)` on a `q_min = 2` target
  silently returns `q == 2`; `q` is validated and then never used
  (`R/dgm.R:334-346`).
- `generate_tte()` silently overwrites incoming `time`/`status`
  columns (`R/tte.R:258-260`).
- `generate_outcomes()` accepts a `dgm_tte` (it inherits `"dgm"`, the
  documented type) and dies with `no applicable method for 'cov_at'`.
  `dgm_mixture()` guards this case with a clear message; this does
  not.
- `compute_performance()` returns 11 rows per method normally but
  **1 row** when fewer than two replicates converge
  (`R/performance.R:76-78`). `sim_power()` then filters for
  `measure == "rejection"` and gets zero rows, producing a zero-row
  `power_curve` and an obscure failure inside `sample_size()`'s
  `approx`.
- Partial matching defeats the otherwise-good `...` hygiene: `sed=`
  is correctly rejected, but **`s = 99L` silently binds to `seed`**
  through `do.call` (verified independently: results identical to
  `seed = 99L`).

## 4. Help system

**Counts** (verified by the subagent, programmatically):

| Metric | Value |
|---|---|
| Exports | 38 (plus 12 registered S3 methods) |
| Documented | 38/38 |
| `@param` vs formals mismatches | **0** |
| With `@examples` | **14/38 (37%)** |
| `@return` complete and accurate | 32/38 |
| Doc-vs-code mismatches | 5 |
| Orphan exports (no example, no vignette) | **0** |
| Vignettes rendering | 6/6 |
| Example blocks executing | 14/14 |

The `@param`-vs-formals result is genuinely clean and was verified by
parsing `tools::Rd_db()` against `formals()` for all 38, not by
sampling.

**Accuracy defects found** (all verified):

1. **`?fit_result` miscounts its own contract.** "These seven names
   are the whole contract" (`R/simulate.R:70`); the list has **eight**
   elements (`estimate, se, p_value, ci_lower, ci_upper, df,
   converged, level`). This is the package's central contract topic
   and was rewritten today.
2. **`?null_fit` is factually wrong.** "A `fit_result` of `NA`s"
   (`R/simulate.R:96`); actual `df = Inf`, `converged = FALSE`,
   `level = 0.95`.
3. `?mask_from_data` documents no return elements at all, and its
   `spec` differs materially from `dropout_mask`'s (`mechanism =
   "empirical"`, `monotone` auto-detected, `psi*` all NA) — none
   documented.
4. `?dropout_mask` names its `spec` attribute but describes none of
   its 11 elements, including `expected`, the achieved calibration a
   user needs as a diagnostic.
5. `?fit_ancova` never mentions its five `null_fit()` early-return
   branches or the multi-arm warning added today.

**Cross-references are the weakest area** (verified over the Rd link
graph): **one** function-level `@seealso` in the whole package, **zero**
`@family` tags, and **15 of 37 topics have no outbound links** —
including `trial_schedule` and `runin_design`, the two entry points of
stage 1. `dropout_mechanism` and `null_fit` are true islands (no
inbound, no outbound). **Nothing links back to `?zzrctsim`**, so a
user landing on a leaf topic from a search has no route to the
workflow map. Eight exports are absent from the map itself.

## 5. User interface

### 5.1 First-use walkthrough (verified, run verbatim from README)

The README block runs to completion: **12 lines to a first real
result**, 30 to the performance table. Friction points, in the order
a new user meets them:

1. **`beta = c(x_slope = 0.5, x_trt_active = -0.25)`** — these magic
   column names are minted by `runin_design()`
   (`R/schedule.R:256,275`). The README never shows them and never
   prints the design. **A new user cannot write step 3 without
   printing `head(d)` or reading a second topic.** This is the
   single largest friction point in the package.
2. `dgm_conditional(G = diag(c(9, 0.04)))` — the ordering (intercept
   variance, then slope variance) is fixed by the default
   `z = function(t) cbind(1, t)`, which the README never mentions.
3. **Step 4 prints 26 lines of raw list** ending in
   `attr(,"class") [1] "fit_result"`, because there is no
   `print.fit_result` (Section 5.2).
4. The user hand-derives `theta`; nothing checks it.
5. `compute_performance(res)` prints a bare 11-row frame.

**The README's flagship `sample_size()` example is degenerate**
(verified): it warns, and returns an upper bound, because every grid
point sits at power ≥ 0.987 against a target of 0.80. The first
substantial thing a reader runs emits a warning about their own grid.
Fix: `n_grid = c(10, 15, 20, 30)`.

### 5.2 `fit_result` has no print method

**Verified**: `exists("print.fit_result", asNamespace("zzrctsim"))` is
`FALSE`; `methods(class = "fit_result")` is empty. Nine other S3
classes have print methods. The object the entire API is organized
around — "the whole contract" — is the only one that dumps its
internals, and it is what the README displays at step 4. One-line
fix.

### 5.3 Argument naming: the `dat`/`design` claim, re-tested

Earlier today I asserted that `dat` (requires `y`) versus `design`
(no outcome yet) is a principled distinction rather than an
inconsistency, and recorded that in a commit message. **That claim was
re-tested this pass and it does not fully hold** (verified):

- It **fails at `apply_mask()`**, which accepts a frame with no `y`
  and fabricates one (Section 3.3).
- More importantly, **the binary is not a binary**:
  `mask_from_data(observed = ...)` is a *third* name for the
  with-outcome concept, and `generate_tte(subjects = ...)` a *third*
  for the no-outcome concept (verified: it accepts and ignores a `y`
  column).

The distinction is defensible; the vocabulary implementing it is not.
I was too quick to close this. **Fix:** rename `observed` → `dat` and
`subjects` → `design`, keeping the old names as deprecated aliases
while that is still free.

### 5.4 Other interface findings

- **`B` is argument 1 of `run_simulation()` and argument 7 of
  `sim_power()`.** The source comment added today asserts "the two
  functions read the same way" — true for `analyze`→`estimand`, which
  is what was harmonized, but the comment overstates.
- **`print.sample_size()` hard-codes two arms** ("total 100 for two
  arms", `R/power.R:260`) though `design_fn` is arbitrary and
  three-arm designs are supported.
- **Unbounded printing**: `print.trial_schedule` prints 124 lines for
  `trial_schedule(treatment = 120)`; `print.power_curve` and
  `print.sample_size` print whole curves.
- **Verified good**: all nine print methods return `invisible(x)`;
  pipe-friendliness holds (`dat |> apply_mask() |> fit_ancova()`
  works); `cov_at()` is type-stable across all four DGM shapes; the
  export surface has nothing that should be internal.

## 6. Coding practices

**Verified sound.** Idiom is spotless: no `sapply`, no `T`/`F`, no
`1:n` outside prose, no `library()`/`require()` in `R/`, `drop =
FALSE` in all 11 risky places, `vapply()` with explicit `FUN.VALUE`
throughout, `<-` everywhere. All 12 S3 methods registered; dispatch
via `inherits()` with no `class(x) ==` anywhere. All 38 exports have
at least one test. Internals are called with the `zzrctsim:::` prefix
in all six places they appear — the defect class that broke the check
this morning is fully cleared.

**RNG discipline verified correct, including today's fix.**
`sim_streams()` and `with_rng_state()` restore `RNGkind()` *and*
`.Random.seed` bit-identically, with the `on.exit` ordering correct.
The `sample_size()` seed-offset claim was **verified empirically**:
with `seed = 123L` the confirmation reproduced `sim_power(seed =
124L)` exactly (0.355) and not `seed = 123L` (0.360). One blemish:
calling `sim_streams()` in a session with no pre-existing
`.Random.seed` leaves one behind in `globalenv()` (`R/rng.R:44-46`).

**Remaining practice gaps:**

- **41 exported entry points still emit `stopifnot()` expression-echo
  errors** (verified by calling each), e.g.
  `certify(matrix(1:6,2,3))` → `nrow(V) == ncol(V) is not TRUE`.
  Today's validation work fixed the specific functions on the (b)
  list; the sweep was never generalized. File:line list is in the
  subagent record; the concentration is `dgm.R` (5), `tte.R` (4),
  `schedule.R` (3), `mixture.R` (3), `missing.R` (3), `simulate.R`
  (3). The package writes excellent messages elsewhere
  (`R/dgm.R:341`, `R/missing.R:31`), which makes the inconsistency
  more jarring, not less.
- Dead code: `widths` at `R/tte.R:114` (found independently by both
  `lintr` and the help audit).
- Performance smells beyond 3.4: `which(id == i)` scans at
  `R/missing.R:167,385` and `R/generate.R:145`; vector growth at
  `R/tte.R:274`; double subsetting at `R/performance.R:62`.

**Test thin spots** (verified): `close_out(rule = "fixed")` and its
`at = NULL` error path are untested; `reach_error()` has one call
with default `q`, which is why 3.6 survived; no test constructs a
tied or non-monotone power curve, which is why 3.2 survived.

## 7. Prioritized checklist

### (a) CRAN blockers

None. The check is clean. Add `NEWS.md` and `Language: en-US` before
submitting, and run the platform checks in Section 8.

### (b) Bugs to fix before anyone depends on the behavior

1. `sample_size()` tie-breaking and non-monotone inversion (3.2) —
   wrong answer in the flagship deliverable.
2. `apply_mask()` response-column guard (3.3) — silent total data
   loss.
3. `dropout_mask()` performance (3.4) — blocks the headline use case
   at realistic n.
4. Resolve the staggered-accrual claim (3.1): implement the route or
   delete the promise.
5. `run_simulation()` warning trapping (3.5).
6. `reach_error()` validation (3.6); the 3.7 silent surprises.

### (c) Documentation completion

7. Correct `?fit_result` (eight, not seven), `?null_fit`,
   `?mask_from_data`, `?dropout_mask`'s `spec`, `?fit_ancova`'s
   early returns.
8. Fix the README `sample_size` grid; add `head(d)` so the `x_*`
   names are visible before `beta` is written.
9. Add `@family`/`@seealso`; give the 15 terminal topics an outbound
   link and every leaf a route back to `?zzrctsim`.
10. Raise `@examples` coverage above 37%, prioritizing
    `dgm_conditional`, `trial_schedule`, `runin_design`, `estimand`,
    `fit_ancova`, `compute_performance`.
11. `behaviour`/`centres` in `vignettes/missing-data.Rmd`.

### (d) Design decisions still cheapest now

12. Add `print.fit_result` (5.2).
13. Rename `observed` → `dat` and `subjects` → `design` (5.3).
14. Bound the unbounded print methods; drop the two-arm hard-coding.
15. Sweep the remaining 41 `stopifnot()` messages.
16. Reject partial-matched `...` names (`s=` binding to `seed=`).

## 8. Not evaluated

- **Platforms**: only macOS 26.5.1 / aarch64 / R 4.6.1. **No**
  win-builder, R-hub, R-devel, Linux, or older-R check was run.
  `cran-comments.md` still lists these as to-do, correctly.
- **`goodpractice::gp()`** was not installed and did not run.
- **Vignette prose** was checked only for execution, not accuracy.
- `@return` completeness for the 32 passing functions rests on the
  happy path plus selective error branches; not exhaustive.
- No reverse-dependency check (none exist).
- The `.Rd` `\doi{}` targets were not resolved over the network.
- Numerical correctness of the statistical methods themselves (the
  DGM conversions, the certificate, the TTE inversion) was not
  re-derived; it is covered by the existing suite, which passes.

## 9. Revision history

- **2026-08-14 (earlier, `ui-help-cran-evaluation-2026-08-14.md`)**:
  First evaluation. `R CMD check --as-cran` returned 2 ERRORs and 4
  NOTEs. Identified a broken shipped example, a test suite that
  failed under check because of unqualified internals, an MIT-template
  LICENSE stub under a GPL-3 declaration, zzcollab scaffolding in the
  tarball, stale citation metadata, an undeclared `utils` import, plus
  functional bugs (dropped `beta` attribute breaking
  `reference_based()`, MAR/MNAR mislabeling, a non-independent
  confirmation seed) and help-system gaps (no `?zzrctsim`, no README,
  15 incomplete `@return`s, `R/power.R` untested). Verdict: not ready.
- **2026-08-14 (remediation, commits `a0de18f`, `af8622d`, `40fdbb9`,
  `e0e58db`)**: All four lists (a)–(d) completed. Check reached 0/0/1.
  Coverage of `R/power.R` went from none to 45 assertions.
- **2026-08-14 (this document, second full review)**: Independently
  re-verified the prior remediation — the RNG restoration, the
  seed-offset fix, the `generate_outcomes()` attribute fix, the
  `zzrctsim:::` qualification, and the `@param`/formals agreement all
  confirmed by execution, not taken on trust. **Newly found**: a
  documented staggered-accrual workflow that cannot execute (3.1);
  `sample_size()` returning the largest rather than smallest tied size
  (3.2); `apply_mask()` fabricating an all-`NA` response (3.3); a
  superlinear `dropout_mask()` (3.4); undocumented warning behavior
  (3.5); `reach_error()` NA/negative-variance regime (3.6); no
  `print.fit_result` (5.2); partial-matching through `...` (5.7).
  **Corrected from the prior review**: the `dat`/`design` naming
  distinction, which I closed this morning as "principled, not a
  defect," does not fully hold — it fails at `apply_mask()` and there
  are four names for the two concepts, not two (5.3). **Still open
  from the prior review**: the 41 remaining bare `stopifnot()`
  messages, and the `approx()` monotonicity assumption in
  `sample_size()`, which I flagged this morning as a judgment call and
  which has now produced a demonstrable wrong answer. Verdict:
  CRAN-ready on the check, not ready as a user-facing release.
