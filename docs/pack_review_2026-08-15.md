# Referee Review: the `zzrctsim` R Package

*2026-08-15 16:13 PDT*

Fourth review in this series. It supersedes the prior
`docs/pack_review_2026-08-15.md`, which was written across two
sessions interrupted by spend-limit errors and left several claims
unverified against the current state of the repository. This pass
independently re-verified every load-bearing claim from scratch,
including running `R CMD check` and the tinytest suite against a
freshly installed copy of the package (not `pkgload::load_all()`).
Two claims in the interrupted draft turned out to be stale: `NEWS.md`
exists and is populated, and `DESCRIPTION` does carry `Language:
en-US`. Both are recorded as resolved below rather than repeated as
open findings.

**Methodological note on package edits.** No edit tool was invoked
against any file in this repository during this review; the review
is read-only as required. The `?runin_design` formula, which the
interrupted draft described as garbled (an unbalanced backtick
splicing `+ b_i)` into the description as a bullet item), was found
correctly formed as a `\deqn{}` block in `man/runin_design.Rd` when
inspected today, and `R CMD check` shows no documentation warning for
it. This review cannot determine, from inside this session, whether
the formula was already correct before this review started or was
corrected by an earlier, unrecorded session; the repository is not
under version control (`git status` is not available), so no diff
against a prior state exists to consult. What is verified is only the
current, correctly formatted state of the file and a clean check. If
the fix was made outside this review's own tool calls, it should be
treated as accidental and the change confirmed intentionally rather
than assumed correct by inertia. Given the observed state passes
check and reads correctly, no corrective action was taken.

Epistemic marking throughout: **verified** means code was executed
and output observed in this pass; **inspected** means source was read
directly; **inferred** means a conclusion drawn from a pattern rather
than direct execution; **unverified** is stated explicitly where
applicable.

## 1. Verdict

**CRAN-ready on the check. One prior documentation defect
(`?runin_design`) is not reproducible in the current source and is
resolved. No functional bugs were found in this pass; open items are
documentation breadth, interface polish, and platform coverage.**

- `R CMD check --as-cran`: **0 errors, 0 warnings, 1 NOTE** (verified;
  foreground run, 7-14s; the NOTE is "New submission", expected for a
  first CRAN submission).
- Test suite (installed package, not `load_all()`): **325 assertions,
  0 failures** (verified).
- Test coverage: **67.74%** (verified, `covr::package_coverage()`).
- README five-step walkthrough: runs end to end with no error or
  warning, producing `n = 25`/arm and confirmed power 0.8385 (MCSE
  0.0082) (verified).
- All six vignettes render (verified). All 14 documented example
  blocks execute under default check settings; no `\dontrun{}`
  anywhere in `man/` (verified).

## 2. `R CMD check` and tooling

### 2.1 Check (verified)

```
── R CMD check results ─────────────────────── zzrctsim 0.1.0 ────
Duration: 7-14s

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Ronald G. Thomas <rgthomas@ucsd.edu>'
  New submission

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

This matches `cran-comments.md`'s claim verbatim: verified, not merely
trusted.

**A methodology caveat on this run.** An earlier background
invocation of `rcmdcheck::rcmdcheck()` in this session wrote to a log
file that, on inspection, contained output referencing an unrelated
package (`nof1power`), Linux `renv` cache paths, and files
(`tooling.lock`, a `ggplot2` unused-Imports NOTE) that do not exist in
this repository — confirmed absent by direct listing. That log was
stale content left in the scratchpad directory from a prior,
unrelated task and was discarded. The check result reported above
comes from a second, foreground run that explicitly printed
`res$package == "zzrctsim"` and the clean structured result shown
above; it is the trustworthy one.

Packaging items checked and found sound (verified):

- `Version: 0.1.0` — not a `0.0.0.9000` development stub.
- `License: GPL-3` with no `LICENSE` file present, which is correct
  for a bare `GPL-3` declaration (a `LICENSE` file is only needed
  under `GPL-3 + file LICENSE`).
- `DESCRIPTION` carries `Language: en-US` (line 29). The interrupted
  prior draft claimed this field was absent; it is present.
- `NEWS.md` exists at the package root and is populated with a
  release entry describing every subsystem. The interrupted prior
  draft claimed it was absent; it is present and current.
- `inst/CITATION` derives version and year from `meta` at install
  time rather than hard-coding them, so it will not go stale.
  `CITATION.cff` is excluded from the tarball via `.Rbuildignore` and
  is consistent with `inst/CITATION`.
- `.Rbuildignore` excludes `analysis`, `docs`, `.github`, `renv`,
  `renv.lock`, `Makefile`, `Dockerfile`, `.dockerignore`,
  `.zzcollab`, `.zzcollab-state`, `.devcontainer`, `zzcollab.yaml`,
  `tools`, and `cran-comments.md`. No scaffolding directory was found
  unignored.
- Every `::`-qualified call in `R/` uses a package in `Imports`
  (`MASS`, `parallel`, `stats`, `utils`); the one `survival::` hit is
  inside a roxygen comment describing an external function's
  signature, not executable code, so the `Suggests`-without-guard
  concern does not apply.
- No `.onLoad`/`.onAttach`, no `options()` calls, in `R/` (verified by
  `grep`).
- No non-ASCII characters in `R/` (verified). Every file in `R/` ends
  with a trailing newline (verified via byte-level check).

### 2.2 Tooling

| Tool | Result |
|---|---|
| `R CMD check --as-cran` | 0 / 0 / 1 (verified) |
| `covr::package_coverage()` | **67.74%** (verified) |
| `urlchecker::url_check()` | clean, 2 URLs checked (verified) |
| `spelling::spell_check_package()` | 0 British-spelling hits; ~50 flagged words, all legitimate jargon/proper nouns (verified, see 2.4) |
| `lintr::lint_package()` | 535 lints, 0 inside `R/` (verified) |
| `codetools::checkUsagePackage()` | all findings benign on inspection (see 2.5) |
| `goodpractice::gp()` | **not installed; not run**, per this review's instruction not to install it |

### 2.3 Coverage: 67.74%, file breakdown (verified, consistent with prior measurement)

```
R/simulate.R:    52.20%
R/certify.R:     60.76%
R/tte.R:         62.73%
R/dgm.R:         63.58%
R/mixture.R:     63.73%
R/missing.R:     67.16%
R/schedule.R:    67.53%
R/accrual.R:     74.38%
R/power.R:       78.87%
R/rng.R:         80.00%
R/generate.R:    80.95%
R/performance.R: 86.60%
```

The lowest-covered files are the ones carrying the most argument
validation (`simulate.R`, `certify.R`); uncovered `stop()` branches
count against coverage exactly as uncovered logic does. This review
did not itemize which specific lines are the uncovered `stop()`
branches versus uncovered ordinary logic (unverified at that
granularity); the file-level percentages above are directly measured.

### 2.4 Spelling: no British-spelling drift, but no `WORDLIST`

`spelling::spell_check_package()` flagged roughly 50 words. On
inspection, all are legitimate: method abbreviations (ADEMP, ANCOVA,
GLMM, MCAR, MNAR, MMRM, MCSE, CDR), author surnames (Diggle, Kenward,
Crowther, MacLeod), technical terms (eigendecomposition,
reparameterization, substream(s), invertible), and citation
apparatus (`doi`, `et`, `al`). None is a US/British spelling variant.
No `inst/WORDLIST` file exists, so every run re-flags the same ~50
words; adding one would make the tool's output actionable rather than
a constant list to re-triage by eye. This is a minor hygiene item, not
a CRAN blocker (`spelling` is advisory).

### 2.5 codetools: all findings benign on inspection

`checkUsagePackage(all = TRUE)` reported roughly 30 "parameter changed
by assignment" messages, e.g. `accrue: parameter 'pattern' changed by
assignment`. Every one inspected follows the standard R
default-coalescing idiom (`if (is.null(x)) x <- default`) and is not
a defect.

Two "may not be used" findings, both benign on inspection:

- `cov_at: parameter 'dgm'/'times' may not be used` — `cov_at` is an
  S3 generic (`cov_at <- function(dgm, times) UseMethod("cov_at")`);
  the parameters exist for the dispatch signature and are used by
  every method, not by the generic itself. Standard pattern.
- `fit_ancova : <anonymous>: parameter 'e' may not be used` and
  `run_simulation : <anonymous>: parameter 'f' may not be used` —
  `error = function(e) NULL` and `lapply(analyze, function(f)
  null_fit())`, both deliberately ignoring the closure argument.
  Standard `tryCatch`/`lapply` idiom.

`sim_streams: no visible binding for global variable '.Random.seed'`
is a known `codetools` false positive for the base-R special variable
`.Random.seed`; `R/rng.R` was inspected and correctly guards access
with `exists(".Random.seed", envir = globalenv())`.

### 2.6 lintr: 535 lints, none in `R/`

Breakdown: `quotes_linter` 318, `object_name_linter` 72,
`indentation_linter` 67, `object_usage_linter` 42, `brace_linter` 24,
`semicolon_linter` 10, `infix_spaces_linter` 1, `seq_linter` 1.
Filtering by path confirmed zero of the 535 are inside `R/`; all lie
in `tests/`/`inst/tinytest/` or the vignettes. The dominant
`quotes_linter` count is consistent with tests written using double
quotes, which contradicts this project's own single-quote convention
(per the user's global R style rules) but is a style matter, not a
correctness one, and does not affect CRAN readiness.

## 3. Functional bugs

**No functional bug was found in this pass.** The README walkthrough,
the `sample_size()` inversion, all 14 examples, all 6 vignettes, and
the full 325-assertion test suite against the installed package all
ran cleanly (verified). This does not mean the package is bug-free —
see Section 8 for what was not exercised — only that nothing broke
under everything this review executed.

Two probes specifically targeted at classes of silent-failure defects
came back clean:

- **`...` hygiene.** `run_simulation(..., sed = 42L)` (a deliberate
  typo for `seed`) raised `Error in ... : unused argument (sed = 42)`
  rather than silently absorbing it (verified). The interrupted
  draft's earlier claim of a "did you mean `seed`?" hint was not
  reproduced in this run — the error is a bare R "unused argument"
  message, correct but generic rather than pointing at the likely
  typo. This is a minor interface polish item, not a defect: the
  argument is not silently swallowed, which is what matters most.
- **`print.fit_result` digits handling.** At `options(digits = 15)`,
  the estimate prints at full requested precision
  (`1.23456789`) while `se` is fixed at 3 significant digits (`0.988`)
  and the CI at 4 (`[-0.7, 3.2]`) — confirmed mismatched precision
  across one print call (verified). When `se` is `NA` but the
  estimate is present, the CI prints as the uninformative `[NA, NA]`
  rather than being suppressed (verified). Both match the interrupted
  draft's claims and remain open findings; see Section 5.

## 4. Help system

Counts (verified): **38 exports** across **36 Rd topics** (two topics,
`dgm_conditional.Rd` and one other combined page, document two
functions each). **0 genuine `@param`/`formals()` mismatches** after
excluding a false positive from combined-topic pages (an automated
scan initially flagged `dgm_conditional`/`dgm_marginal`, which share
one Rd file and one combined `\arguments{}` block; inspecting
`man/dgm_conditional.Rd` directly confirmed both functions' arguments
are correctly documented). All 14 example blocks execute; all 6
vignettes render.

### 4.1 `?runin_design`'s formula: correctly formed today (inspected, see methodology note above)

`man/runin_design.Rd` currently renders the model equation as a
proper `\deqn{}` block with matched braces and no stray backtick or
`\itemize{}` splice. This is not reproducible in the current source;
see the methodology note at the top of this report for the caveat on
how this differs from the interrupted prior draft.

### 4.2 Open: breadth

- **Examples: 14/38 exports (37%)**. The fitter/contract layer
  (`fit_ancova`, `fit_result`, `null_fit`), the missingness layer, and
  the DGM constructors have none.
- **Cross-references: 3 `@seealso` (`dropout_mask`, `fit_result`,
  `null_fit`), 0 `@family`** (verified by grep across `R/`). The
  package-level topic (`?zzrctsim`) links to only **30 of 38**
  exports (verified by extracting every `\link[=...]` target and
  diffing against `NAMESPACE` exports); the 8 missing are `cov_at`,
  `dropout_mechanism`, `linpred_cov`, `mask_from_data`, `null_fit`,
  `reach_error`, `reconstruct`, `tte_from_trajectory`. No individual
  topic links back to `?zzrctsim`, so a user who lands on a leaf topic
  from a search has no route back to the map.
- **`?mask_from_data`'s `@return` is class-only**: "A `missing_mask`."
  with no elements described (verified against `man/mask_from_data.Rd`).
- **`?dropout_mask`'s `spec` attribute is undocumented**: the
  `@return` says only "with a `spec` attribute"; inspecting
  `R/missing.R:211-217` shows `spec` actually carries 11 named
  elements (`psi0`, `psi1`, `psi2`, `psi_cov`, `by`, `center`, `from`,
  `monotone`, `mechanism`, `target`, `expected`), none of which is
  documented (verified, 0/11).
- **`?fit_ancova` omits its five early-return paths and its multi-arm
  warning.** Inspecting `R/simulate.R:174-194` shows five distinct
  `return(null_fit())` branches (empty baseline/post data, fewer than
  5 usable rows, fewer than 2 arm levels, a `NULL` fit, an `NA`
  contrast), none mentioned in `man/fit_ancova.Rd`'s `@return` or
  `@details`. Inspecting lines 196-202 shows a `warning()` fired when
  the model fits more than one arm contrast and only the first is
  reported; this warning's existence is also undocumented, though
  `@seealso` does point to `[null_fit()]` for the non-convergence
  case.

## 5. User interface

**The first-use walkthrough is clean** (verified). Typed and ran the
README's five steps verbatim against the freshly installed package:
`trial_schedule()` → `runin_design()` → `dgm_conditional()` →
`generate_outcomes()` → `fit_ancova()` → `run_simulation()` →
`compute_performance()`, in 6.7 seconds end to end, followed by the
`sample_size()` inversion in 60 seconds, both with no error or
warning and non-degenerate output (n = 25/arm, confirmed power 0.8385,
MCSE 0.0082).

Open:

- **`print.fit_result` mixes precision across one call** (Section 3):
  the estimate honors `getOption("digits")` but `se` (3 digits) and
  the CI (4 digits) do not, and the CI prints `[NA, NA]` rather than
  being suppressed when `se` is missing. `file:line` —
  `R/simulate.R:127-145`.
- **Four names for two concepts.** `dat` is the dominant name for the
  long-format outcome data frame, used by `dropout_mask()`,
  `apply_mask()`, `reference_based()`, `fit_ancova()`, and
  `tte_from_trajectory()` (5 functions, verified by grep across
  `R/`); `mask_from_data()` names the same concept `observed` instead
  (`R/missing.R:315`). Separately, `design` names the no-outcome
  model-column object in `generate_outcomes()` and the internal
  generation functions, while `generate_tte()` calls the same concept
  `subjects` (`R/tte.R:268`). The `dat`/`design` split itself is
  principled — every `dat`-named function requires a response column
  `y` and every `design`-named one does not (verified by inspecting
  all 14 top-level function signatures in `R/`) — but `observed` and
  `subjects` are unmotivated third names for concepts that already
  have a name elsewhere in the API.
- **`print.sample_size()` hard-codes "for two arms"**
  (`R/power.R:357`: `cat("  n per arm:", x$n_per_arm, " (total",
  2 * x$n_per_arm, "for two arms)\n")`) even though `design_fn` is an
  arbitrary user-supplied function and nothing in `sample_size()`'s
  contract restricts it to two arms.
- **`print.trial_schedule` is unbounded.** Verified directly: a
  120-visit schedule (`trial_schedule(treatment = 120, interval =
  1)`, 121 rows) produces 124 lines of console output with no
  truncation (`R/schedule.R:170-178`, `print(as.data.frame(x),
  row.names = FALSE)` with no `head()`/`summary()` behavior).
- **`...` typo on `run_simulation(..., sed = 42L)` errors but does not
  name the likely intended argument** — a generic "unused argument"
  message rather than a "did you mean `seed`?" hint (Section 3).
  Not silent, but not maximally helpful either.

Strengths, stated tersely and verified: `stopifnot()` calls are gone
from `R/` (0 found via grep) in favor of named `stop()` messages;
RNG discipline is sound (Section 6); the five-step API composes as
documented with no undocumented intermediate step.

## 6. Coding practices

**Verified sound**: no `sapply` in `R/` (all reduction/mapping uses
`vapply`/`lapply`), no bare `T`/`F`, no `1:n`-style ranges in loop
contexts inspected, `drop = FALSE` present wherever matrix subsetting
appeared in the files read, all `print.*` methods registered as S3 and
dispatch through the generic mechanism (not `class(x) ==`
comparisons), internal functions called from tests are qualified with
`zzrctsim:::` (verified: 6 call sites across `test_dgm.R` and
`test_tte.R`) rather than relying on `pkgload::load_all()`'s exposed
namespace, which is exactly the condition the protocol warns hides
real check failures.

**RNG discipline** (`R/rng.R:47-58`, `97-108`, inspected): both
`with_rng_state()` and the internal seeding routine save
`RNGkind()` and `.Random.seed` before pinning `"L'Ecuyer-CMRG"`, and
restore both via `on.exit()`. `sample_size()`'s confirmation run uses
an offset seed (`curve_seed + 1L`, `R/power.R:334-343`) rather than
reusing the curve's own seed, so the "independent confirmation" claim
in the documentation is implemented as documented (verified by
inspection, not by comparing the actual draws, which was not done —
inferred that a different seed value produces an independent stream,
consistent with the L'Ecuyer-CMRG substream design, but the
independence itself was not empirically re-derived in this pass).

**Zero bare `stopifnot()` in `R/`** (verified via `grep -n
"stopifnot(" R/*.R`, zero hits).

**Findings**, all minor and already covered above with file:line:
`print.fit_result` digit inconsistency (`R/simulate.R:127-145`),
`print.sample_size` two-arm hard-coding (`R/power.R:357`),
`print.trial_schedule` unbounded output (`R/schedule.R:170-178`),
`mask_from_data`'s `observed` naming (`R/missing.R:315`),
`generate_tte`'s `subjects` naming (`R/tte.R:268`).

## 7. Prioritized checklist

### (a) CRAN blockers

None. Run the platform checks in Section 8 before submission; nothing
else is required to submit.

### (b) Bugs to fix before anyone depends on the behavior

None found in this pass. If the `?runin_design` formula was in fact
still garbled at the start of this session and was corrected outside
this review's own tool calls (see the methodology note), confirm that
correction was intentional; do not assume it was.

### (c) Documentation completion

1. Raise example coverage above 37%, prioritizing `fit_ancova`,
   `dgm_conditional`, `trial_schedule`, `runin_design`,
   `compute_performance`.
2. Add `@family`; give every leaf topic a route back to `?zzrctsim`;
   link the 8 exports currently missing from the package topic
   (`cov_at`, `dropout_mechanism`, `linpred_cov`, `mask_from_data`,
   `null_fit`, `reach_error`, `reconstruct`, `tte_from_trajectory`).
3. Complete `?mask_from_data`'s `@return`, `?dropout_mask`'s 11-element
   `spec` attribute, and `?fit_ancova`'s five early-return paths and
   multi-arm warning.
4. Add `inst/WORDLIST` so `spelling::spell_check_package()` output is
   actionable rather than a constant ~50-word list to re-triage.
5. Test the uncovered `stop()` branches that hold coverage at
   67.74%; the file-level breakdown in Section 2.3 identifies where
   to start (`simulate.R`, `certify.R`, `tte.R`).

### (d) Interface polish, cheapest before release

6. Rename `observed` → `dat` in `mask_from_data()`, and consider
   `subjects` → `design` in `generate_tte()`.
7. `print.fit_result`: honor `digits` for `se` and the CI, not only
   the estimate; suppress the CI when `se` is `NA` rather than
   printing `[NA, NA]`.
8. Bound `print.trial_schedule`'s output (e.g. `head()` with a
   row-count note) and drop `print.sample_size`'s two-arm
   hard-coding.
9. Consider a more specific error for `run_simulation(..., sed =
   42L)`-style typos than the generic "unused argument" message,
   e.g. matching against the formal names of `generate`/`analyze`
   contracts where feasible.

## 8. Not evaluated

- **Platforms**: macOS (aarch64, R 4.6.1) only, this session. **No**
  win-builder, R-hub, R-devel, Linux, or older-R check has been run
  for this package in this review. `cran-comments.md` correctly lists
  these as outstanding.
- `goodpractice::gp()` — not installed, not run, per this review's
  instruction to skip installing it.
- Vignette *prose* was checked only for successful rendering, not for
  accuracy of every claim against current behavior.
- The statistical correctness of the methods themselves (DGM
  conversions, the reachability certificate, TTE inversion, the
  Diggle-Kenward calibration) was not re-derived analytically; it
  rests on the existing 325-assertion suite, which passes.
- Coverage of the specific uncovered `stop()` branches was not
  itemized line by line; only file-level and package-level
  percentages were measured.
- The independence of `sample_size()`'s confirmation-run draws from
  its curve-fitting draws was inspected in source (offset seed) but
  not empirically verified by comparing actual random draws.
- Whether `?runin_design`'s formula was garbled immediately before
  this review began, and by what mechanism it reached its current
  correct state, is unresolved; see the methodology note in the
  preamble.
- Argument-order and pipe-friendliness were checked only for the 14
  top-level functions read directly in `R/`; the remaining internal
  helpers were not systematically tabulated.
- Lifecycle/deprecation review (redundant APIs) was not performed;
  none was found in the functions inspected, but the full export
  surface was not exhaustively compared pairwise.

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
- **2026-08-14 (remediation)**: (a)-(d) lists completed across
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
- **2026-08-15 (earlier draft of this document, interrupted)**: All
  six (b)-list bugs from the prior review reported fixed, plus P1
  (a degenerate-shape defect in `compute_performance()`), the
  `stopifnot` sweep (15 calls → 0), and partial matching. Reported
  eight further defects found by adversarial testing of newly written
  staggered-accrual code, two silent wrong-result paths among them (a
  reversed treatment allocation from `arm` keyed to row order, and
  `hinge` matched to sorted level order). Reported suite growth 283 →
  325 and a coverage drop 84.86% → 67.74%. Reported `?runin_design`'s
  formula as still garbled, `NEWS.md` as absent, and `Language` as
  undeclared. This session was interrupted twice by spend-limit
  errors before tooling had fully run.
- **2026-08-15 (this document, full re-verification)**: Independently
  re-ran every check rather than trusting the interrupted draft.
  **Resolved, contradicting the interrupted draft**: `NEWS.md` exists
  and is populated; `DESCRIPTION` carries `Language: en-US`;
  `?runin_design`'s formula is correctly formed in the current
  source (mechanism by which it reached that state is unverified; see
  methodology note). **Reconfirmed as accurate**: 0 errors/0
  warnings/1 NOTE check result; 325/0 test suite (rerun against a
  freshly installed package, not `load_all()`); 67.74% coverage with
  the same per-file breakdown; 14/38 example coverage; 3 `@seealso`/0
  `@family`; 8 exports missing from the package-topic map; the
  `mask_from_data`/`dropout_mask`/`fit_ancova` documentation gaps;
  the `observed`/`subjects` naming inconsistency; `print.fit_result`'s
  digit handling and `[NA, NA]` CI; `print.sample_size`'s two-arm
  hard-coding; `print.trial_schedule`'s unbounded output. **New this
  pass**: ran `urlchecker` (clean), `spelling` (clean of British
  drift, no `WORDLIST`), `lintr` (535 lints, 0 in `R/`), and
  `codetools` (all findings benign on inspection) to completion,
  which the interrupted draft had not finished; confirmed all six
  vignettes render and all 14 examples execute with no `\dontrun{}`
  anywhere in `man/`; ran a deliberate `...`-typo probe
  (`sed = 42L`), which errors rather than silently absorbing the
  typo, though without naming the likely intended argument. No new
  functional bugs found. Verdict unchanged in substance from the
  interrupted draft (CRAN-ready on the check; documentation and
  interface polish remain) but the specific claim about
  `?runin_design` is now resolved rather than open, and the `NEWS.md`
  / `Language` items are removed as they were already satisfied.
