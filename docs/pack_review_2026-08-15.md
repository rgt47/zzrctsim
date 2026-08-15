# Referee Review: the `zzrctsim` R Package

*2026-08-15 12:10 PDT*

Third review in this series, and the first written after a full
remediation cycle. It supersedes `docs/pack_review_2026-08-14.md`,
which is now substantially stale: that document predates twelve
defect fixes and lists as open a number of items since closed. The
earlier `docs/ui-help-cran-evaluation-2026-08-14.md` remains the
original evaluation. Revision history for all three is carried
forward in Section 9.

Epistemic marking throughout: **verified** means code was executed
and output observed in this pass; **inspected** means source was read
directly. Subagent findings are marked with which of the two the
agent reported, and every load-bearing claim was re-run independently.

## 1. Verdict

**CRAN-ready on the check; close to ready as a user-facing release,
with one documentation defect that should not ship.**

- `R CMD check --as-cran`: **0 errors, 0 warnings, 1 NOTE** (verified;
  47.7s; the NOTE is "New submission").
- Test suite: **325 assertions, 0 failures** (verified), up from 283
  at the start of the cycle.
- Test coverage: **67.74%** (verified, `covr`) — **down** from 84.86%.
  See Section 2.3; this is a real effect and not a regression in
  behavior.
- Every known silent-wrong-result path is closed.

The headline change since the last review is that all six (b)-list
bugs are fixed, plus P1 (the degenerate-shape defect), the
`stopifnot` sweep, and partial matching. Against that, hostile
testing of the newly written staggered-accrual code found **eight**
further defects in a single day, two of which would have produced
wrong scientific results with no diagnostic. That ratio is the most
important thing in this review and is discussed in Section 6.

## 2. `R CMD check` and tooling

### 2.1 Check (verified)

0 errors, 0 warnings, 1 NOTE across every run today. All Suggests
(`simr`, `lme4`, `survival`) are installed here, so no vignette or
test was silently skipped.

Outstanding packaging items, both trivial and both still open
(verified):

- **`NEWS.md` is absent.** A first-release entry is expected.
- **DESCRIPTION has no `Language` field.** `spelling` defaults to
  `en-US` and says so; declare it.

### 2.2 Tooling

| Tool | Result |
|---|---|
| `covr::package_coverage()` | **67.74%** (verified) |
| `R CMD check --as-cran` | 0 / 0 / 1 (verified) |
| `urlchecker` | clean (verified, earlier today) |
| `spelling` | clean of British spellings (verified; the two vignette hits are fixed) |
| `lintr` | `R/` near-clean; the 276 quote lints are all in tests and contradict the project's single-quote standard |
| `codetools` | clean apart from benign `match.arg` reassignments |
| `goodpractice::gp()` | **not installed; did not run** |

### 2.3 Coverage fell 17 points while the suite grew — why

This looks alarming and is not. Coverage went 84.86% → **67.74%**
while assertions went 283 → 325. The cause is that this cycle added a
large volume of *error-handling* code — 15 bare `stopifnot()` calls
became roughly 40 argument-specific `stop()` branches, plus the
ragged-path guards — and **most of those branches are never
triggered by a test**. Uncovered `stop()` lines count against
coverage exactly as uncovered logic does.

Worst-hit files are precisely the ones that gained the most
validation: `R/simulate.R` 69.1% → **52.2%**, `R/rng.R` 100% →
**80.0%**, `R/certify.R` 80.0% → **60.8%**.

Two honest readings, and both are true. The error messages are much
better than they were, which is a real user-facing improvement that
coverage cannot see. And the new branches are genuinely untested, so
a typo inside a `stop()` string — or a guard whose condition is
inverted — would not be caught. The second reading is the actionable
one: **the error paths need tests**, and until they have them the
67.74% is the honest number.

## 3. Functional bugs

**No known silent-wrong-result path remains open.** This section
records what was closed, because the pattern matters more than any
individual fix.

### 3.1 Closed this cycle (all verified by execution)

From the 2026-08-14 review's (b) list:

1. **Staggered accrual could not execute.** `?generate_outcomes`
   documented a route through `runin_design()` that errored.
   Implemented rather than withdrawn, since DESCRIPTION headlines the
   feature: `realize_schedule()` now returns a classed object
   carrying the phase boundaries, and `runin_design()` accepts both
   shapes.
2. **`sample_size()` returned the largest tied size**, a 2×
   overstatement on a saturating curve, against a documented contract
   of "the smallest".
3. **`apply_mask()` fabricated an all-`NA` response** on a mistyped
   `response=`, destroying the outcome silently.
4. **`dropout_mask()` was superlinear** — ~10× faster at n = 1600,
   with output verified bit-identical across six paths.
5. **`run_simulation()` documented warning-trapping it did not do**,
   so a fitter warning emitted B console warnings and recorded none.
6. **`reach_error()` returned silent `NA`s and negative variances**;
   `as_conditional()` ignored `q`; `generate_tte()` overwrote columns.

And, found during this cycle's own re-review:

7. **P1, the degenerate-shape defect.** `compute_performance()`
   returned 1 row instead of 11 when fewer than two replicates
   converged, so `sim_power()`'s filter matched nothing. Two symptoms:
   a crash (`arguments imply differing number of rows: 1, 0`) and,
   worse, a **multi-method curve that silently dropped the degenerate
   method** — a comparison of two fitters could report one. Fixed at
   the root by keeping the row shape stable; the follow-on `if (NA)`
   in `sample_size()` was fixed with it.

### 3.2 Found by attacking code written the same day (all verified)

The staggered-accrual path was new, so it was targeted adversarially.
Eight defects, in the order found:

- **`arm` keyed to row order.** `unique(schedule$id)` returns
  first-appearance order, so reordering rows — *sorting by calendar
  date, the natural thing to do in a staggered design* — silently
  inverted the entire treatment allocation.
- **`sort()` was not enough.** The follow-up fix repaired integer ids
  only: sorting a factor orders by level, and character ids put
  `"S10"` before `"S2"`. The root error was matching by **rank**
  rather than identity; the branch now refuses non-integer ids.
- **`hinge`'s guard was unreachable.** `names(hinge) <- levs` ran
  before the length check, so a wrong-length vector died inside
  `names<-`, and a **right-length unnamed vector was silently matched
  to sorted level order** — `c(TRUE, FALSE, TRUE)` for arms
  low/hi/pbo applied as hi/low/pbo. Silent model misspecification.
- **Hard crash on a documented combination**: a subject enrolled past
  close-out keeps no visits, and `data.frame()` cannot recycle a
  scalar `id` against zero-length columns. Reachable via
  `close_out(rule = "fixed")` and via `accrue(pattern = "poisson")`,
  which is documented as able to overrun `period`.
- Column-subsetting a `realized_schedule` kept the class but dropped
  the attributes, failing with "missing value where TRUE/FALSE
  needed"; `realize_schedule()` validated neither `enroll` nor
  `close`.

Every one of these passed `R CMD check` cleanly.

## 4. Help system

Counts (verified): **38 exports**, all documented; **0 `@param`
vs `formals()` mismatches**; all example blocks execute; all six
vignettes render.

### 4.1 Open: `?runin_design`'s description is garbled (verified)

Still present, and it should not ship. The inline model formula was
parsed as markdown, so `man/runin_design.Rd` renders as:

```
Constructs the three slope columns of the model
`Y_ij = alpha + a_i + (delta (1 - h_j) + beta h_j + gamma h_j g_i
\itemize{
\item b_i) t_j + w_ij`, expanded over subjects.
```

An unbalanced backtick with `+ b_i)` spliced in as a bullet, on the
first paragraph of a central function's help page. Fix: escape the
formula as `\eqn{}` or restructure so `+ b_i)` does not begin a line.

### 4.2 Open: breadth

- **Examples: 14/38 exports (37%)** — unchanged all cycle. The
  fitter/contract layer (`fit_ancova`, `fit_result`, `null_fit`), the
  whole missingness layer, and the DGM constructors have none.
- **Cross-references: 3 `@seealso`, 0 `@family`** (was 1 and 0).
  Improved elsewhere — 77 internal link edges, terminal topics 15 →
  12, islands 2 → 1 — but **no topic links back to `?zzrctsim`**, so
  a user landing on a leaf from a search has no route to the map, and
  8 exports are absent from the map itself.
- `?mask_from_data` still describes no return elements;
  `?dropout_mask`'s `spec` attribute lists none of its 11; `?fit_ancova`
  omits its five `null_fit()` early returns and its multi-arm warning.

### 4.3 Closed

`?fit_result` (eight elements, not seven), `?null_fit` (`df`,
`converged`, `level` are not `NA`), `?generate_outcomes` (no longer
claims an impossible route), `?runin_design`'s `@param`/`@return`
for the ragged case, and `?dropout_mask`'s missing `@details` — all
verified corrected.

## 5. User interface

**The first-use walkthrough is now clean** (verified). The README
runs end to end; `names(d)` shows the `x_*` columns before the reader
must name them in `beta`, which was the largest friction point; and
the `sample_size` example returns a non-degenerate answer — n = 25/arm,
confirmed power 0.8385 (MCSE 0.0082) — with **no warning**, versus
the saturated upper-bound answer it gave this morning.

Closed: `print.fit_result` now exists, so the README's step 4 prints
four readable lines instead of 26 raw ones; partial matching through
`...` is refused with a "did you mean `seed`?" hint; all 15 bare
`stopifnot()` calls are gone, so entry-point errors name the argument
and the fix.

Open:

- **`print.fit_result` respects `getOption("digits")` for the
  estimate only** — `se` is fixed at 3 digits, the CI at 4 — so at
  `digits = 15` the estimate and its SE print at visibly mismatched
  precision. It also emits a meaningless `[NA, NA]` interval when the
  estimate is present but the SE is missing.
- **Four names for two concepts.** `observed` (`mask_from_data`) is a
  third name for the same long frame everything else calls `dat`;
  `subjects` (`generate_tte`) a third for the no-outcome concept. The
  `dat`/`design` split is principled — verified that every `dat`
  function requires `y` — but `observed` is a straight inconsistency.
- `print.sample_size()` hard-codes "for two arms" though `design_fn`
  is arbitrary; `print.trial_schedule` prints 124 lines for a
  120-visit schedule.

## 6. Coding practices, and the pattern worth naming

**Verified sound**: idiom is spotless (no `sapply`, no `T`/`F`, no
`1:n`, `drop = FALSE` throughout, `vapply` with `FUN.VALUE`); all S3
methods registered, dispatch via `inherits()`; internals called with
`zzrctsim:::` everywhere; RNG discipline correct, including the
`sample_size()` seed offset verified empirically rather than read.
Zero bare `stopifnot()` remain in `R/`.

**The pattern.** Across this cycle, eleven defects were found in code
written the same day, and `R CMD check` passed cleanly against every
one. Three rounds of increasingly hostile testing were needed for the
`arm`-mapping bug alone: the first fix (`unique` → `sort(unique)`)
was necessary but insufficient, because ranking ids was itself the
wrong mechanism and sorting merely made it wrong in fewer cases.

Two lessons follow, and they generalize past this package:

1. **The newest code deserves the most hostile testing.** Its happy
   path is the only one anyone has run.
2. **A fix that addresses the symptom can mask the root cause.**
   Sorting made the failure rarer and therefore harder to find. The
   question to ask of any fix is not "does the reproduction pass now"
   but "was that the mechanism, or a special case of it".

## 7. Prioritized checklist

### (a) CRAN blockers

None. Add `NEWS.md` and `Language: en-US`, then run the platform
checks in Section 8.

### (b) Should not ship

1. Fix `?runin_design`'s garbled description (4.1).

### (c) Documentation completion

2. Test the new error branches — this is what the 17-point coverage
   drop is telling you (2.3), and it is the highest-value item here.
3. Raise example coverage above 37%, prioritizing `fit_ancova`,
   `dgm_conditional`, `trial_schedule`, `runin_design`,
   `compute_performance`.
4. Add `@family`; give every leaf topic a route back to `?zzrctsim`;
   link the 8 missing exports from the package topic.
5. Complete `?mask_from_data`, `?dropout_mask`'s `spec`,
   `?fit_ancova`'s early returns.

### (d) Interface polish, cheapest before release

6. Rename `observed` → `dat` (and consider `subjects` → `design`).
7. `print.fit_result`: honor `digits` throughout; suppress the CI
   when the SE is `NA`.
8. Bound the unbounded print methods; drop the two-arm hard-coding.

## 8. Not evaluated

- **Platforms**: macOS 26.5.1 / aarch64 / R 4.6.1 only. **No**
  win-builder, R-hub, R-devel, Linux, or older-R check has ever been
  run for this package. `cran-comments.md` correctly lists these as
  outstanding.
- `goodpractice::gp()` — not installed, never ran.
- Vignette *prose* was checked only for execution, not for accuracy
  against behavior.
- The statistical correctness of the methods themselves (DGM
  conversions, the reachability certificate, TTE inversion) was not
  re-derived; it rests on the existing suite, which passes.
- Coverage of the *new* error branches specifically was not
  itemized; only the package-level and per-file percentages were
  measured.

## 9. Revision history

- **2026-08-14 (`ui-help-cran-evaluation-2026-08-14.md`)**: First
  evaluation. Check returned 2 ERRORs, 4 NOTEs: a broken shipped
  example, a suite that failed under check from unqualified
  internals, an MIT-template LICENSE under a GPL-3 declaration,
  zzcollab scaffolding in the tarball, stale citation metadata, an
  undeclared `utils` import. Plus functional bugs (dropped `beta`
  attribute, MAR/MNAR mislabeling, non-independent confirmation seed)
  and help gaps (no `?zzrctsim`, no README, 15 incomplete `@return`,
  `R/power.R` untested). Verdict: not ready.
- **2026-08-14 (remediation)**: (a)–(d) lists completed across
  `a0de18f`, `af8622d`, `40fdbb9`, `e0e58db`. Check reached 0/0/1.
- **2026-08-14 (`pack_review_2026-08-14.md`)**: Second review.
  Independently re-verified the remediation, then found a documented
  staggered-accrual workflow that could not execute, `sample_size()`
  returning the largest tied size, `apply_mask()` fabricating an
  all-`NA` column, a superlinear `dropout_mask()`, undocumented
  warning behavior, `reach_error()`'s NA regime, no
  `print.fit_result`, and partial matching through `...`. Corrected
  its own prior claim that `dat`/`design` was fully principled.
  Verdict: CRAN-ready on the check, not ready as a release.
- **2026-08-15 (this document)**: All six (b)-list bugs fixed, plus
  P1, the `stopifnot` sweep (15 calls → 0), and partial matching.
  Adversarial testing of the new staggered-accrual code found eight
  further defects, two of them silent wrong-result paths (a reversed
  treatment allocation from `arm` keyed to row order, and `hinge`
  matched to sorted level order); all fixed with regression tests.
  Suite 283 → 325. **Newly measured**: coverage fell 84.86% →
  67.74%, caused by ~40 new untested `stop()` branches — the error
  messages improved, the error *paths* are untested. **Still open**:
  `?runin_design`'s garbled description (should not ship), example
  coverage at 14/38, no inbound links to `?zzrctsim`, `observed`
  naming, `NEWS.md`, `Language`, and no platform check beyond macOS.
  Verdict: CRAN-ready on the check; one documentation defect to fix
  before release.
