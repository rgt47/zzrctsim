# zzrctsim: Requirements and Design Rationale
*2026-08-10 12:19 PDT*

## 1. Purpose

`zzrctsim` is a single trial-simulation engine intended to serve the
simulation needs of the research compendia under `~/prj/res/`. It
absorbs the existing `trialsim` package (`res/28-trialsim`) rather than
sitting alongside it.

The design criteria stated by the author are **accurate, flexible,
fast**, in that order of importance. Section 8 records what each of
those means operationally.

## 2. Evidence base

A survey of the 35 compendia under `~/prj/res/` on 2026-08-10 found
**21 that are trial-simulation dependent**. Most have an empty `R/`
directory: the manuscripts are scaffolded and the simulation code is
not yet written. This is the central practical opportunity. A shared
engine introduced now prevents 21 divergent reimplementations rather
than having to reconcile them afterwards.

The compendia already converge informally. Their reports share a
section template -- Introduction, Present study, Methods (ADEMP
structure, Data-generating model, Parameter values, Analysis models,
Simulation procedure), Results -- and at least three (`13`, `14`, `18`)
carry byte-identical comments pinning `RNGkind` and a single
`set.seed()` per Morris et al. (2019) §4.1. The convention exists; what
is missing is the code that implements it once.

## 3. The core design thesis

The survey overturned the assumption that motivated `trialsim`.

`trialsim` generates **marginally**: it draws the full response vector
from a multivariate normal with a directly specified covariance. The
portfolio, however, is dominated by **conditional** (random-effects)
data-generating models:

| Compendium | Data-generating model |
|---|---|
| 05-optimal-visit-placement | random intercept + random slope |
| 07-multicenter-rct | random site intercept + random treatment x site |
| 08-mmrm-linearity-robust | random intercept + slope, nonlinear effect `g(t, k)` |
| 09-summary-stats-efficacy | random intercept + slope |
| 11-site-covariate-analysis | site intercept + treatment x site + subject effect |
| 17-mixed-r2 | two-level LMM, components solved from target R^2 and ICC |
| 31-long-power | random intercept + slope |
| 30-missforest | marginal MVN + covariate |
| 28-trialsim | marginal MVN over (baseline, change scores) |
| 02-adaptive-alloc | Weibull proportional hazards + covariates |
| 13-mmrm-vs-survival-ad | latent continuous + time-to-event, jointly |

`trialsim` is the outlier, not the template.

More consequentially, the **relationship between the two
parameterizations is itself the subject matter** of several
manuscripts. `18-mmrm-vs-lm-simple` asks when random-effects models
reduce to simpler forms. `04-runin-power-analysis` uses the Woodbury
identity to move between conditional and marginal variance
expressions. `09-summary-stats-efficacy` and `31-long-power` both turn
on the same mapping.

**Therefore the central object of `zzrctsim` is a data-generating
specification that can be declared in either parameterization, with
exact conversion in both directions where one exists, and a computable
certificate reporting when it does not.**

Concretely:

- Declare conditionally as `(Z, G, sigma^2)` and obtain the induced
  marginal `V = Z G Z' + sigma^2 I` exactly.
- Declare marginally as `V` and obtain the minimal conditional
  representation, or a certificate that none exists at the requested
  rank.
- The certificate is the eigenvalue argument documented in
  `res/28-trialsim/docs/rogers-process-reconstruction-whitepaper.md`
  §13: a target `V` is reachable from `q` random effects if and only if
  some `sigma^2 >= 0` makes `V - sigma^2 I` positive semi-definite with
  rank at most `q`. Take `sigma^2` as the smallest eigenvalue of `V`;
  the resulting rank is the minimum `q`. Compound symmetry certifies at
  `q = 1`; AR(1) with `p = 6` certifies at `q = 5`.

No existing package does this. `simr` and `lme4` are conditional only;
`mmrm` and `nlme` are marginal only.

### 3.1 What the certificate is, and what it is not

The certificate is **constructive**. Eigendecompose the target `V`, set
`sigma^2` to its smallest eigenvalue, and count the strictly positive
eigenvalues of `V - sigma^2 I`; that count is `q_min`. The
corresponding eigenvectors supply `Z` and the shifted eigenvalues
supply `G`. Verified reconstruction to machine precision:

```
compound symmetry      p=6  q_min=1  sigma^2=0.5000  max|V-(ZGZ'+s2 I)|=8.9e-16
AR(1)                  p=6  q_min=5  sigma^2=0.2654  max|V-(ZGZ'+s2 I)|=1.8e-15
HC hetero SD x AR(1)   p=6  q_min=5  sigma^2=6.5151  max|V-(ZGZ'+s2 I)|=7.5e-14
```

**`q_min` is necessary, not sufficient for a particular model.** The
certificate permits `Z` to be arbitrary, namely the eigenvectors. A
real random-effects specification fixes `Z` to design columns --
intercept `1`, slope `t`, and so on. Compare leading eigenvectors:

```
CS  :  0.408  0.408  0.408  0.408  0.408  0.408   flat, = 1/sqrt(6)
AR1 : -0.324 -0.419 -0.469 -0.469 -0.419 -0.324   bowed
```

Compound symmetry's is exactly flat, which is why a random intercept
(forcing `Z = 1`) reproduces it exactly. AR(1)'s is bowed, so even the
optimal one-dimensional direction is unavailable to a random intercept.
The certificate therefore proves impossibility cleanly -- no model with
`q < q_min` can represent the target -- but does not prove possibility
for a given `Z`. Any use of it must state that asymmetry.

### 3.2 Novelty: a caution

The underlying linear algebra is classical. Writing
`V = Lambda Lambda' + sigma^2 I` and solving by eigendecomposition is
**probabilistic PCA** (Tipping and Bishop, 1999), and more generally
factor analysis with isotropic uniquenesses. There is no new
mathematics here, and a referee will recognise it immediately. This
must not be presented as a theoretical contribution.

What is defensible is narrower:

1. Packaging the check as an **automatic step** in a trial-simulation
   workflow, so that a user who supplies an empirical covariance is
   told immediately whether the random-effects analysis they intend can
   represent it at all.
2. Using it as the **objective axis of a software comparison**,
   replacing "expressible / not expressible" judgment calls with a
   computed integer that is identical in form for every package in the
   `lme4` family.
3. The accompanying **error-magnitude table** (28% relative Frobenius
   error at `q = 1` for AR(1) with `p = 6`), which converts
   "unreachable" into a statement of how far wrong a practitioner would
   be who tried anyway.

That is a methods *section* of the software paper, or a short applied
note. It is not a standalone theoretical result.

## 4. Endpoint families

Five, per the author's decision:

1. **Gaussian longitudinal** -- the continuous repeated-measures case.
   Required by 04, 05, 07, 08, 09, 11, 12, 14, 17, 18, 21, 28, 30, 31.
2. **Time-to-event** -- implemented. Weibull and piecewise-exponential
   baseline hazards, single and recurrent events, Gaussian frailty, a
   delayed-effect option giving non-proportional hazards, staggered
   accrual with administrative censoring, and an endpoint derivable
   from a generated longitudinal trajectory. Required by 02 (survival
   arm), 10 (log-rank), 13, 20. Not implemented: accelerated failure
   time, competing risks, cure models, and baselines other than
   Weibull or piecewise exponential.
3. **Binary** -- required by 10 (exact conditional power at small event
   rates), 19.
4. **Count** -- Poisson and negative binomial, including recurrent
   events. Required by 10 (counts arm).
5. **Finite mixtures of multivariate normals** -- a component label is
   drawn per subject and the response from that component's mean and
   covariance. Required by 29 (Gaussian mixture models for MRI-based
   progression), and useful wherever a population is a blend of
   progressor types rather than one homogeneous group.
6. **Mixtures of alpha-stable laws**, of which the Gaussian mixture is
   the `alpha = 2` case. See §4.1.

Families 5 and 6 were previously listed as one item. Separating them
matters because they are not comparably difficult. A normal mixture
has finite moments, so the covariance machinery, the reachability
certificate and the ADEMP performance measures all continue to apply
unchanged; it is an additive feature. A stable mixture with
`alpha < 2` has no variance, which withdraws the foundation those
three rest on (§4.2); it is a structural change. Grouping them made
the easy case look as blocked as the hard one.

**Implementation status as of 2026-08-11.** Families 1 to 5 are built:
Gaussian, the exponential-family set (binomial, Poisson, negative
binomial) via `family=` on the conditional generator, time-to-event
via `dgm_tte()`, and normal mixtures via `dgm_mixture()`. Family 6 is
not built, and remains a requirement rather than a feature. Compendia
25 and 26 cannot migrate until it is; 29 now can.

### 4.0 Should the families be unified as an exponential family?

The natural suggestion is to collapse Gaussian, binary and count into
a single generator with a `family` and `link` argument, as `glm()` and
`glmmTMB()` do. The conclusion is that this is worth doing for part of
the problem but must not become the organizing principle.

**What it unifies, and what it does not.** The exponential family
unifies *univariate marginal* distributions. It supplies no
multivariate analogue: there is no canonical multivariate binomial or
multivariate Poisson corresponding to the multivariate normal. For
Gaussian data the marginal family and the dependence structure are the
same object, the covariance matrix. For every other family they
decouple, and a separate dependence mechanism must be chosen. That
choice, not the marginal law, is the design decision.

**Two routes, mirroring the dichotomy of §3.**

*Conditional (GLMM).* Draw `b_i ~ N(0, G)`, then
`y_ij | b_i ~ ExpFam(g^-1(x'beta + z'b_i))`, conditionally independent
given `b_i`. This generalizes cleanly to any exponential family and is
what `lme4`, `glmmTMB` and `simr` do. Two costs: the dependence
reachable is exactly what random effects can induce, so the constraint
of §3 applies with no marginal escape hatch; and under a non-identity
link the marginal mean is not `g^-1(x'beta)`, so conditional and
marginal treatment effects differ and the estimand must say which is
intended.

*Marginal (copula, NORTA).* Specify arbitrary marginals plus a
correlation matrix and push a latent multivariate normal through. This
is the `simstudy::genCorGen` and `SimCorMultRes` approach. The
difficulty is that the requested correlation is not the achieved one.
Measured, latent Gaussian correlation against realized Pearson
correlation, `n = 2e5`:

```
marginal                 latent achieved retained
Gaussian                   0.60    0.600     100%
Poisson mu=3               0.60    0.578      96%
Bernoulli p=0.5            0.60    0.410      68%
Bernoulli p=0.10           0.60    0.321      53%
Bernoulli p=0.02           0.60    0.209      35%
```

Counts survive nearly intact. Binary does not, and the loss worsens as
events become rarer: at `p = 0.02` only 35% of the requested
correlation is realized. That is precisely the regime of compendium
10, exact conditional power at small event rates. An interface that
accepted a correlation matrix for binary outcomes and silently
delivered a third of it would be a defect of the same character as the
uncentred dropout logit of §4.6. If a copula route is offered, the
**achieved** correlation must be computed and returned, not assumed.

**Coverage.** The abstraction reaches three of the five families:

| Family | Exponential family? |
|---|---|
| Gaussian | yes |
| Binary | yes |
| Count (Poisson, negative binomial) | yes |
| Time-to-event | no: censoring, hazards and risk sets are a different structure |
| alpha-stable mixtures | no: no closed-form density, not an exponential family |

Two of the five committed families fall outside it, so it cannot be
the top-level abstraction.

**What it would break.** `cov_at()` presupposes second moments and a
covariance-parameterized joint law; for binary outcomes the variance
is determined by the mean and the two cannot be specified
independently. `certify()` and `as_conditional()` are Gaussian
results. The residual-preserving construction in `reference_based()`
assumes additive residuals. The centred hazard in `dropout_mask()`
reads the response on a continuous scale and is meaningless for
`y` in `{0, 1}`.

**Decision.**

1. Add `family` and `link` to the *conditional* generator. Standard,
   matches `glmmTMB`, covers three families at low cost.
2. Offer a copula layer for the marginal route, reporting the achieved
   correlation alongside the requested one.
3. Do **not** make the exponential family the central abstraction.
   Treat the endpoint families as sibling generators behind one
   contract.

The justification for (3) is that the surrounding machinery is already
distribution-free. `trial_schedule()`, `accrue()`,
`realize_schedule()`, `apply_mask()`, `run_simulation()`,
`compute_performance()`, `sim_power()` and `sample_size()` make no
distributional assumption whatever. The reuse across families lives in
that spine, which is built, not in a unified generator, which would
cover only three of five cases. The correct framing is **five
generators behind one contract**, not one generator with five
settings.

### 4.1 Mixtures of stable laws

Taking the fifth family to be *mixtures* of stable laws rather than a
single stable law is a substantial generalisation, and it unifies four
compendia rather than one:

| Compendium | Connection |
|---|---|
| 26-ecf-estimation | stable-law parameters via the empirical characteristic function |
| 25-fourier-gmm | Fourier density estimation applied to mixtures |
| 29-mixnormalmri | Gaussian mixture models for MRI-based progression |
| all Gaussian-endpoint repos | the `alpha = 2` degenerate case |

The stable family `S(alpha, beta, gamma, delta)` has no closed-form
density except in three cases -- Gaussian (`alpha = 2`), Cauchy
(`alpha = 1, beta = 0`), and Levy (`alpha = 1/2, beta = 1`) -- but its
characteristic function is closed-form throughout. That is precisely
why 26 estimates through the ECF, and why simulation is the practical
route to everything else.

Two distinct constructions are needed and must not be conflated.

**(a) Finite mixtures.** `sum_k w_k S(alpha_k, beta_k, gamma_k,
delta_k)`. Gaussian mixtures are the case `alpha_k = 2` for all `k`,
and a single Gaussian is the degenerate one-component case. This is the
marginal/univariate construction and serves 25, 26, and 29 directly.

**(b) Sub-Gaussian alpha-stable, for multivariate dependence.**
`X = sqrt(W) . Z` where `Z ~ N(0, Sigma)` and `W` is positive stable
with index `alpha/2`. This is a scale mixture of normals, it is
elliptical, and critically it **retains a shape matrix `Sigma`** that
plays the structural role a covariance plays in the Gaussian case.
This is the tractable route to longitudinal dependence under heavy
tails; the general multivariate stable law requires a spectral measure
on the unit sphere and is not tractable for routine simulation.

Verified numerically (`alpha = 1.5`, `p = 4`, AR(1) shape `rho = 0.6`):

```
  sample variance per visit : 85.0  50.9  46.6  56.9   <- non-convergent
  target shape row 1        : 1.000 0.600 0.360 0.216
  robust proxy row 1        : 1.000 0.569 0.340 0.186
```

The co-variation structure is preserved in the shape matrix even though
the variance does not exist.

### 4.2 Architectural consequence: moments may not exist

This is the most important constraint the stable family imposes, and it
reaches well beyond the generator.

A stable law has finite variance **only** when `alpha = 2`, and finite
mean only when `alpha > 1`. For `alpha < 2`, therefore:

- The covariance machinery of §3 is undefined. `V = D R D`, the
  reachability certificate, and any conversion between conditional and
  marginal parameterisations all presuppose second moments. Under
  heavy tails the analogous object is a **shape matrix**, and the
  engine must keep the two concepts distinct rather than overloading
  one name.
- MMRM, GLS, and REML are not applicable as specified.
- **The ADEMP performance measures break.** Bias, empirical standard
  error, and mean squared error are undefined without moments, and
  their Monte Carlo standard errors are meaningless. Median bias,
  interquartile-range-based dispersion, and coverage -- which survives,
  being a probability -- are the required substitutes. `compute_performance()`
  as inherited from `trialsim` must therefore dispatch on whether
  moments exist, and must refuse rather than silently return a
  finite-looking number.

Verification metrics need the same treatment: the Wishart z-scores and
the Stein/LRT statistic in `test_dgm_verification.R` are Gaussian
results and do not transfer. Characteristic-function-based checks are
the natural analogue, and are already the estimation tool in repo 26.

### 4.3 Parameterisation traps

Two, both verified, both capable of silently corrupting a study.

**The `alpha = 2` scale convention.** Drawing with `alpha = 2,
gamma = 3, delta = 5` yields a sample SD of 4.246, not 3:

```
  gamma = 3          sqrt(2) * gamma = 4.243        observed sd = 4.246
```

So `S(2, beta, gamma, delta) = N(delta, 2 * gamma^2)`. A user who sets
`gamma` to their intended standard deviation gets data inflated by
`sqrt(2)`. The engine must expose a Gaussian-native parameterisation
and perform this conversion internally.

**Nolan's parameterisations.** `stabledist` takes a `pm` argument
selecting among the 0, 1, and 2 conventions, which differ in how the
location parameter is centred and are not interchangeable near
`alpha = 1`. The convention must be pinned once, documented, and
asserted in tests.

### 4.4 Available tooling

A CRAN scan on 2026-08-10 found 16 packages referencing stable
distributions. The relevant ones:

| Package | Version | Relevance |
|---|---|---|
| `stabledist` | 0.7-2 | Density, CDF, quantile, RNG. Installed. The baseline dependency. |
| `mixSSG` | 2.1.1 | Clustering via **mixtures of sub-Gaussian stable** distributions -- construction (b) above |
| `MixStable` | 0.1.0 | Parameter estimation for stable distributions **and their mixtures** |
| `StableEstim` | 2.4 | Four-parameter estimation, including ECF methods -- directly relevant to repo 26 |
| `libstable4u` | 1.0.5 | Fast C implementations |
| `alphastable`, `FMStable`, `TempStable` | | Inference, finite-moment and tempered variants |

### 4.5 Build versus reuse: assessment of `mixSSG` and `MixStable`

Both were examined on 2026-08-10 (function inventories from the CRAN
reference manuals, metadata from `tools::CRAN_package_db()`).

**The decisive observation is that these are estimation packages.**
`zzrctsim` needs *generation*. In every one of them the generation
surface is a small corner of a large inferential package, and the
generation itself is the easy part.

#### `mixSSG` 2.1.1 (Teimouri, 2022-09-11)

Eight exported objects. Two generate:

- `rpstable()` -- positive stable random variables.
- `rssg()` -- **skewed sub-Gaussian stable random vectors**. This is
  exactly construction (b) of §4.1, and it supplies skewness, which
  the symmetric construction does not.

The remainder estimate or cluster (`fitmssg`, `fitBayes`, `stoch`) or
approximate the density (`dssg`). Dependencies are light: `ars`,
`MASS`, `rootSolve`. Risk: single maintainer, no release since
September 2022.

#### `MixStable` 0.1.0 (Manou-Abi, Najib, Slaoui, 2025-11-03)

Roughly 150 exported objects. Generation is confined to `rstable()`,
`generate_alpha_stable_mixture()`, `generate_mixture_data()`,
`generate_synthetic_data()`, and `simulate_mixture()`, all
**univariate**. The bulk is estimation: quantile/McCulloch, maximum
likelihood, SEM/EM, Bayesian Gibbs, Metropolis-Hastings, and importance
sampling.

Four concerns argue against depending on it:

1. **Version 0.1.0, released November 2025.** No track record.
2. **Package hygiene.** The export list includes `mock_gibbs_sampling`,
   `mock_lookup_alpha_beta`, `clip`, `Re`, `Im`, `safe_integrate`, and
   `unpack_params`. Exported mock objects and leaked internal helpers
   indicate the namespace was not curated.
3. **Irrelevant heavy dependencies.** `openxlsx` and `jsonlite` -- an
   Excel writer and a JSON parser -- in a distributions package, plus
   `mixtools`, `nortest`, `e1071`, `libstable4u`.
4. **Embedded unrelated domain code.** `compute_serial_interval()`,
   `plot_effective_reproduction_number()`, `est_r0_ml()`,
   `empirical_r0()`, and a `serial_interval` dataset. This appears to
   be an epidemiological reproduction-number analysis released as a
   general package.

**However, `MixStable` is important prior art for compendium 26.** Its
ECF machinery -- `ecf_empirical()`, `ecf_regression()`,
`ecf_estimate_all()`, `ecf_components()`, `fit_stable_ecf()`,
`robust_ecf_regression()`, `estimate_stable_kernel_ecf()`,
`estimate_stable_recursive_ecf()` -- overlaps directly with what
`26-ecf-estimation` proposes to do. That manuscript's state-of-the-field
section must engage with it. This is a finding for repo 26, not a
dependency decision for `zzrctsim`.

#### `alphastable` 0.2.1 (Teimouri et al., 2019)

Provides `urstab()`, `urstab.trunc()`, `mrstab()`, and
`mrstab.elliptical()` for multivariate elliptically contoured stable
generation, plus `ufitstab.cauchy.mix()`. Same lead author as `mixSSG`.
Unmaintained since 2019.

#### Decision

1. **Import `stabledist`** as the univariate baseline. Maintained by
   the Rmetrics group, minimal, focused.
2. **Build the sub-Gaussian construction rather than depend on it.**
   `X = sqrt(W) . Z`, with `W` positive stable of index `alpha/2` and
   `Z ~ N(0, Sigma)`, is roughly five lines and was verified working
   during this analysis. Taking a dependency on an unmaintained
   single-maintainer package to avoid five lines is the wrong trade,
   and the construction must in any case be integrated with the design,
   dropout, and analysis layers rather than called in isolation.
3. **Suggest, do not Import, `mixSSG` and `alphastable`**, and use
   `rssg()` / `mrstab.elliptical()` in the test suite as independent
   cross-validation of our own generator. This is the highest-value use
   of both: an oracle, not a dependency.
4. **Do not depend on `MixStable`** for the reasons above. Record it as
   prior art for compendium 26.
5. **Revisit if skewness is required.** The symmetric sub-Gaussian
   construction is straightforward; the *skewed* variant implemented in
   `rssg()` is not. If skewed sub-Gaussian endpoints turn out to be
   needed, reassess reuse at that point rather than pre-emptively
   reimplementing.

The general principle this case illustrates, and which should be
applied to the other endpoint families: **for a simulation engine, the
generation half of a distributional package is usually small and worth
owning; the estimation half is large and worth deferring to.**

**Joint generation is required, not merely parallel support.**
`13-mmrm-vs-survival-ad` generates a latent continuous severity process
and derives a time-to-event endpoint from it (time to CDR threshold
crossing). The engine must therefore let a survival endpoint be a
*function of* a generated longitudinal trajectory, not an independently
drawn variable.

## 4.6 Dropout generation: the standard, and current conformance

### The standard

There is a settled convention for generating dropout in longitudinal
simulation studies, and it is the **selection-model factorization**:
write the joint distribution of the response `Y` and the missingness
indicator `R` as `f(Y, R) = f(Y) f(R | Y)`, generate complete data from
`f(Y)`, then generate missingness from `f(R | Y)`. This is why
complete data are generated first and missingness applied as a
separate layer.

Within that factorization the canonical parameterization is the
discrete-time dropout hazard of Diggle and Kenward (1994,
<doi:10.2307/2986113>). For a subject still under observation at visit
`j`, the probability of dropping out at that visit is

```
logit P(drop at j | history) = psi_0 + psi_1 * y_{j-1} + psi_2 * y_j
```

with the mechanism determined by which coefficients are non-zero:

| Mechanism | Restriction |
|---|---|
| MCAR | `psi_1 = psi_2 = 0` |
| MAR | `psi_2 = 0`, `psi_1` free |
| MNAR | `psi_2` free |

The essential point is that the three mechanisms are **nested
restrictions of one model**, not three separate models. MNAR extends
MAR by adding dependence on the current, unobserved response; it does
not replace the dependence on history.

Two further conventions are near-universal in practice:

- **The response entering the hazard is centred**, or the model is
  written on a standardized scale, so that `psi_0` retains its meaning
  as a baseline log-odds.
- **`psi_0` is calibrated** to achieve a target marginal dropout
  proportion, because the realized proportion depends on `psi_1`,
  `psi_2`, and the response distribution jointly. A per-visit hazard
  parameter is not itself the realized dropout rate unless
  `psi_1 = psi_2 = 0`.

The alternative factorization, `f(Y, R) = f(R) f(Y | R)`, is the
pattern-mixture model. It is the natural home for reference-based
assumptions such as jump-to-reference and copy-increments-in-reference
(Carpenter, Roger and Kenward 2013), and it is what `rbmi` implements
for analysis. It is a legitimate alternative generator and should
eventually be offered, but the selection model is the default
convention for simulation.

### Does `apply_dropout()` meet it?

**Partially. It conforms in structure and deviates in two respects,
one of which is a defect.**

**Conforms.** The selection-model factorization is respected: complete
data are generated first and `apply_dropout()` applies missingness
afterwards, so the hazard may depend on realized responses. Dropout is
monotone and absorbing, and a discrete-time logistic hazard is used.

**Deviation 1: MNAR is not the Diggle-Kenward form.** The
implementation treats the mechanisms as mutually exclusive rather than
nested:

```
mar  : plogis(psi_0 + slope * y_{j-1})     # psi_2 = 0        correct
mnar : plogis(psi_0 + slope * y_j)         # psi_1 = 0        NOT standard
```

The standard MNAR retains `psi_1 * y_{j-1}` and *adds* `psi_2 * y_j`.
The current MNAR is the restricted special case `psi_1 = 0`, which is
not what the literature means by MNAR and does not nest MAR. Any
sensitivity analysis contrasting MAR against MNAR under this
implementation is varying two things at once: it removes the
history dependence at the same moment it adds the current-value
dependence.

**Deviation 2: the response is not centred, so `rate` is not the
dropout rate.** Measured on a simulated trial with mean response 22.2,
requesting a per-visit hazard of `0.12`:

| Mechanism | `slope` | Subjects ever dropping |
|---|---|---|
| MCAR | any | 45.5% |
| MAR | 0.05 | 81.3% |
| MAR | 0.25 | **100.0%** |
| MNAR | 0.05 | 82.2% |
| MNAR | 0.25 | **100.0%** |

Because `plogis(qlogis(0.12) + 0.25 * 20) = 0.95`, the hazard
saturates at clinically realistic ADAS-Cog or CDR-SB values and every
subject drops out at the first eligible visit. The `rate` argument is
therefore meaningless whenever `slope` is non-zero, and silently so:
no warning is issued and the returned data look plausible.

Note that even under MCAR the realized 45.5% is not 12%: `rate` is a
per-visit hazard, and over five eligible visits the cumulative
proportion is `1 - (1 - 0.12)^5 = 47%`. That is arithmetically correct
but the argument name invites the wrong reading.

### Required changes, and their status

All seven are now implemented.

1. **Done.** Nested Diggle-Kenward form with separate `psi1` and
   `psi2`. The mechanism is derived from the coefficients by
   `dropout_mechanism()` rather than selected by a string, so MNAR
   nests MAR by construction instead of replacing its history term.
2. **Done.** The response entering the hazard is centred, by default
   at the complete-data mean, and the centring value is recorded in
   the returned specification.
3. **Done.** `rate` is replaced by `target`, the proportion of
   subjects with any missing observation, and `psi0` is obtained by
   root-finding. Because the hazards depend only on complete-data
   responses, the expected proportion is available in closed form as
   `1 - prod(1 - p_j)` averaged over subjects, so calibration needs no
   inner Monte Carlo and is exact to the solver tolerance.
4. **Done.** `dropout_mask(target = c(placebo = 0.15, active = 0.35),
   by = "arm")` calibrates `psi0` separately within each group;
   `psi_cov` adds covariate-dependent terms. Measured realized
   proportions 0.123 and 0.290 against targets 0.15 and 0.35.
5. **Done.** Mask generation and mask application are now separate
   layers: `dropout_mask()` produces a `missing_mask` object,
   `mask_from_data()` imports an empirical pattern from a completed
   trial, and `apply_mask()` applies either. A mask is therefore
   saveable, inspectable, and reusable across replicates, which is what
   permits the missingness pattern to be held fixed under common
   random numbers.
6. **Done.** `monotone = FALSE` gives intermittent missingness. At a
   common target of 0.30 the monotone mask removed 535 observations
   from 169 subjects and the intermittent mask 204 from the same 169,
   which is the expected relationship: the same subjects are affected,
   but absorbing dropout costs far more data.
7. **Done.** `reference_based()` implements jump to reference, copy
   increments in reference, copy reference, and last mean carried
   forward as *generators*, modifying the post-discontinuation mean
   while preserving each subject's realized residual so that
   within-subject correlation is retained. Measured mean shifts on the
   treated arm: J2R `+2.660`, CIR `+1.840`, LMCF `-1.104`, CR
   `+2.660`. The ordering `0 < CIR < J2R` is the expected one, since
   copying increments retains the accumulated benefit that jumping to
   reference discards.

`apply_dropout()` is retained as a convenience wrapper over
`dropout_mask()` plus `apply_mask()` for the common single-group case.

### A distinction the implementation makes explicit

`reference_based()` alters the trajectory; `apply_mask()` hides
observations. These are deliberately separate because discontinuation
and missingness are not the same event. Under a treatment-policy
estimand, post-discontinuation values are observed and the
reference-based trajectory is the truth being estimated; under a
hypothetical estimand they are unobserved. Fusing the two would make
the treatment-policy case inexpressible, and it is the case ICH E9(R1)
puts first.

## 5. Design-structure requirements

Traced to the compendia that impose them.

| Requirement | Imposed by |
|---|---|
| Multilevel nesting (patients in sites, clusters) | 07, 11, 16, 17, 36 |
| Randomization procedures: covariate-adaptive, matched-tuple, blocked | 02, 16 |
| Factorial allocation | 16 |
| Multi-arm with dose-response and multiplicity adjustment | 19, 28 |
| Nonlinear mean trajectories | 08, 30 |
| Flexible visit schedules, unequal spacing | 05 |
| Run-in and common-close phases (`J0`, `J1`, `J2`) | 04 |
| Staggered accrual and partial enrollment at interim | 02, 13, 28 |
| Dropout: MCAR, MAR, MNAR; monotone and intermittent | 14, 05, 30 |
| Interim analysis, conditional power, sequential monitoring | 10, 20, 21, 28 |
| Baseline covariates with specified correlation to outcome | 02, 11, 12, 30 |
| Solving variance components from target R^2 / ICC | 17 |
| ADEMP performance measures with Monte Carlo SE | all 21 |
| Single-stream RNG with per-replicate state capture | all 21; explicit in 13, 14, 18 |

Two of these deserve emphasis because no surveyed package provides
them:

- **Run-in / common-close phase structure** (04). The design is
  parameterized as `J0` pre-randomization observations, `J1`
  randomized-phase observations, `J2` common-close observations. This
  is not expressible in any existing simulator's design vocabulary.
- **Target-driven parameter solving** (17). The user specifies a
  marginal R^2 and an ICC, and the variance components are solved
  analytically to hit them. `mlmpower` does this for multilevel models
  and is the reference implementation to study.

## 6. What is absorbed from `trialsim`

Carry forward, having been verified correct this session:

- `build_corr()` / `build_covariance()` -- separation strategy
  `V = D R D`, correctly scaled.
- `mrm_covariance()` -- between/within ICC decomposition with AR(1),
  correctly using `cov2cor()`.
- `joint_covariance()` -- attaches the baseline block to the
  change-score block. This joint `(baseline, change_1..change_p)`
  parameterization is a genuine asset and is not offered by any
  surveyed package.
- `hc_change_sds()`, `hc_baseline_change_cor()` -- Homocysteine
  empirical calibration.
- `compute_performance()` -- Morris Table 6 suite with MCSE.
- Per-replicate RNG state capture in `run_simulation()`.
- The DGM verification suite (`test_dgm_verification.R`, 26
  assertions): structural invariants, Hotelling and Wishart moment
  recovery, Stein/LRT global covariance test, dropout mechanism-label
  tests.

Do **not** carry forward:

- The `# Source:` provenance headers citing `c001/rogers.R` by line
  range. Their origin is a pharmaceutical sponsor's internal code
  (Pfizer protocol B0341002, compound PF-04494700); see the white paper
  §2. `zzrctsim` must be written from the process reconstruction, not
  from that source.
- `generate_trial()`'s row-wise assembly. See section 8.
- The claim of Kenward-Roger degrees of freedom, which `nlme::gls`
  does not provide.

## 7. Architecture

A four-layer separation, mirroring ADEMP so that manuscript sections
map onto code modules:

```
  design   ->  specify arms, visits, phases, accrual, randomization
  dgm      ->  specify mean structure + (conditional | marginal) covariance
               + endpoint family; convert between parameterizations
  generate ->  draw complete data, then apply missingness as a layer
  analyze  ->  pluggable fitters with a uniform result contract
  measure  ->  ADEMP performance with MCSE
```

Three contracts fix the seams:

1. **The DGM object** carries both parameterizations and knows which is
   canonical for a given specification. Conversion is explicit and
   certified, never silent.
2. **Complete data are always generated first**; missingness is a
   separate layer. This permits the information loss from dropout to be
   isolated by comparison, and it is what makes response-dependent
   (MAR, MNAR) mechanisms expressible at all.
3. **Fitters return a uniform shape.** Note that `trialsim` claims this
   but does not achieve it: `fit_ancova()` returns a per-contrast list
   with no confidence interval, contradicting its own documentation.
   The contract must be tested, not asserted.

## 8. Performance requirement

Measured on `trialsim` 0.1.0, two arms, six visits:

```
n_per_arm=  100   0.057 s/call  ->   56.8 s for B=1000
n_per_arm=  400   0.197 s/call  ->  197.0 s for B=1000
n_per_arm= 1600   0.804 s/call  ->  804.4 s for B=1000

MVN draw only (2 arms x 400): 0.0002 s
```

The multivariate normal draw accounts for **0.1%** of runtime. The
remaining 99.9% is one `data.frame()` allocation per subject followed
by `do.call(rbind, ...)`. Vectorized long-format assembly is therefore
worth roughly a thousandfold, and it is the single highest-value change
to make.

Target: generation cost dominated by the draw itself, so that a
power curve over five sample sizes at `B = 1000` completes in seconds
rather than a day.

Secondary performance requirements:

- Draw all `B` replicates in one call where the design is fixed, rather
  than looping.
- Parallel execution must use `L'Ecuyer-CMRG` streams. Naive per-worker
  `set.seed()` breaks the per-replicate reproducibility that is a
  stated requirement of every compendium.
- Store per-replicate estimates, not only summaries, so new performance
  measures can be computed without re-running.

## 9. Non-goals

- Not a group-sequential design package. `rpact` and `gsDesign` are
  comprehensive; `zzrctsim` simulates designs, it does not derive
  boundaries.
- Not a general Monte Carlo harness. `SimDesign` fills that role, and
  building on it should be evaluated rather than reimplemented.
- Not a fitting package. `mmrm`, `nlme`, `lme4`, and `survival` are the
  analysis engines; `zzrctsim` supplies the uniform contract around
  them.
- Not a reference-based imputation package. `rbmi` is the standard for
  MNAR sensitivity analysis under reference-based assumptions.

## 10. Open questions

1. **Dependency posture.** `trialsim` imports only `MASS`, `nlme`,
   `stats`. Covering five endpoint families plausibly requires
   `survival`, and alpha-stable generation requires `stabledist` or an
   equivalent. Decide whether these are Imports or Suggests with
   guarded use.
2. **Build on `SimDesign`?** It already supplies the replicate loop,
   error trapping, seed management, and parallelism that section 8
   requires. The counterargument is dependency weight and idiom
   lock-in. This should be decided deliberately, since a JOSS reviewer
   will ask.
3. **Naming.** `TrialSimulator` (CRAN, 1.20.1, Han Zhang) already
   exists with a closely related purpose. `zzrctsim` avoids collision;
   confirm that the `zz` prefix is intended for a package that may go
   to CRAN.
4. **Relationship to the manuscript.** The reachability certificate is
   not a standalone theoretical result (§3.2) and should not carry a
   statistics-journal submission on its own. It belongs as a methods
   section of the software paper, or as a short applied note on
   covariance specification in simulation studies. The open question is
   therefore what the *manuscript's* contribution is, given that the
   software contribution is the engine itself. Two candidates worth
   weighing:

   - A **software paper** (JOSS for the package; JSS or the R Journal
     if a fuller treatment is wanted), in which the certificate is one
     designed-in feature among several.
   - A **comparative evaluation** across the simulation packages
     surveyed in
     `res/28-trialsim/docs/rogers-process-reconstruction-whitepaper.md`
     §§9-14, using the certificate as the objective axis and the DGM
     verification metrics (§13) as the fidelity axis. This is a
     genuine empirical contribution and is not currently in the
     literature, but it is a different paper from the software paper
     and should not be conflated with it.

## 11. Implemented features

As of 2026-08-10 the package exports 34 functions and passes 153
tinytest assertions. This section documents what is built. Everything
described here was verified by execution; measured figures are quoted
rather than asserted.

### 11.1 Design layer

**`trial_schedule()`** defines the observation schedule across up to
three phases, following the formulation of compendium 04:

```
  run-in         j = -J0..-1   all untreated, t_j < 0
  randomization  j = 0         t_j = 0
  treatment      j = 1..J1
  common close   j > J1        all arms off treatment
```

Two constructors: phase counts with a common `interval`, or explicit
`times` for unequally spaced designs (compendium 05). Returns one row
per visit with `index`, `time`, `phase`, the phase indicator `h`, and
`on_treatment`, carrying `J0`/`J1`/`J2` as attributes. A `print`
method summarizes.

**`runin_design()`** expands a schedule over subjects and builds the
fixed-effect columns. Two design decisions are exposed rather than
buried:

- `common_close = "revert"` (default) sets the treatment column to
  zero during the common close, per compendium 04; `"retain"` holds it
  at its last on-treatment value. Verified: `x_trt` is
  `0 0 0 3 6 9 12 0 0` and `0 0 0 3 6 9 12 12 12` respectively. The
  first is a discontinuity in the mean; the second is the usual
  clinical reading. Neither is silently assumed.
- `hinge` controls, **per arm**, whether the post-randomization slope
  differs from the run-in slope. An unhinged reference arm constrains
  `beta = delta`, dropping one parameter, which is what makes run-in
  observations directly informative about the control slope.

The parameterization is `delta` common to all visits plus per-arm
post-randomization increments. Verified to span the identical column
space as compendium 04's `delta`/`beta`/`gamma` form
(`rank(cbind(orig, new)) == rank(orig) == 3`), so the closed-form
variance derivations there transfer unchanged.

**`accrue()`, `close_out()`, `realize_schedule()`** implement
staggered accrual with a common close-out. Enrolment may be uniform,
linear or Poisson. Under the last-subject-last-visit rule, follow-up
is a decreasing function of enrolment time and early enrollees are
observed *past* the nominal schedule. Verified with five linearly
accrued subjects and a nominal six visits:

```
  visits per subject:  10  9  8  7  6
  follow-up duration:  24 21 18 15 12
```

Nothing occurs after close-out; the last enrollee receives exactly the
nominal schedule. This is the mechanism that determines the
information fraction at an interim analysis.

### 11.2 Data-generating mechanism layer

**Dual parameterization.** `dgm_conditional(G, sigma2, z)` declares
random effects; `dgm_marginal(V, times)` declares the response
covariance directly, either as a matrix on a fixed grid or as a
function of time. `cov_at(dgm, times)` is the single polymorphic
operation both answer.

The function form matters: staggered accrual makes the observed grid
subject-specific, so a marginal specification given as a fixed matrix
cannot be evaluated on an extended grid. `cov_at()` on a matrix-valued
marginal DGM takes the appropriate submatrix for a subset of the grid
and **refuses** off-grid evaluation with a message pointing at the
function form, rather than extrapolating.

**Continuous-time structures.** `cov_ar1()` and `cov_cs()` return
covariance functions. AR(1) uses `rho^|t_i - t_j|`, the exponential
form, so it extends to unequal spacing naturally: verified that
`V[1,4]` equals `sd^2 * 0.6^7` on the grid `c(0, 1, 2.5, 7, 11)`.
Standard deviations may be scalar, vector, or a function of time.

**Conversion and the certificate.** `as_marginal()` is always exact.
`as_conditional()` uses `certify()` and refuses an insufficient
budget rather than silently approximating. `certify()` returns
`q_min`, the minimum number of random effects, together with the `Z`
and `G` attaining it; `reconstruct()` rebuilds the target to machine
precision. `reach_error()` tabulates the best attainable error at each
budget, converting "unreachable" into a magnitude.

Verified: compound symmetry certifies at `q_min = 1` with a flat
leading eigenvector, which is why a random intercept reproduces it
exactly; AR(1) at `p = 6` certifies at `q_min = 5` with a bowed
leading eigenvector, and a single random effect incurs 28% relative
Frobenius error.

### 11.3 Generation layer

**`generate_outcomes()`** draws from the DGM given a design and a
named fixed-effect vector, returning the data with `mu` and `y` added
and `beta`/`intercept` retained as attributes for downstream use.
Subjects sharing a time vector share a covariance and are drawn in one
call, so balanced designs collapse to a single `mvrnorm()`; only
genuinely ragged grids need more.

**Missingness in two composable layers.** A *mask* says which
observations are missing; a *filter* applies one. Keeping them
separate is what allows an externally supplied pattern to be used and
what allows the pattern to be held fixed across replicates.

- `dropout_mask()` generates from the Diggle-Kenward selection model
  `logit P(drop at j) = psi0 + psi1 (y_{j-1} - c) + psi2 (y_j - c)`.
  MCAR, MAR and MNAR are the nested restrictions `psi1 = psi2 = 0`,
  `psi2 = 0`, and `psi2` free, so MNAR extends MAR rather than
  replacing its history term. `dropout_mechanism()` names the implied
  mechanism.
- The response is centred, and `psi0` is calibrated to a **target
  proportion of subjects with any missing observation** by
  root-finding. Because hazards depend only on complete-data
  responses, the expectation is available in closed form as
  `1 - prod(1 - p_j)` averaged over subjects, so calibration needs no
  inner Monte Carlo. Verified to hit targets of 0.10, 0.25 and 0.50
  across MCAR, MAR and MNAR and across a threefold range of `psi1`.
- `by =` gives per-arm calibration; `psi_cov =` gives
  covariate-dependent missingness. Verified realized proportions 0.123
  and 0.290 against per-arm targets 0.15 and 0.35.
- `monotone = FALSE` gives intermittent missingness. At a common
  target of 0.30 the monotone mask removed 535 observations from 169
  subjects and the intermittent mask 204 from the same 169.
- `mask_from_data()` imports an empirical pattern from a completed
  trial and detects whether it is monotone.
- `apply_mask()` applies any mask, and errors if the mask does not
  cover every row.

**`reference_based()`** implements jump-to-reference, copy-increments-
in-reference, copy-reference and last-mean-carried-forward as
*generators*. The post-discontinuation mean is modified while each
subject's realized residual is preserved, so within-subject
correlation survives. Measured mean shifts on the treated arm: J2R
`+2.660`, CIR `+1.840`, LMCF `-1.104`, CR `+2.660`; the ordering
`0 < CIR < J2R` is the expected one.

Discontinuation and missingness are deliberately separate operations,
because under a treatment-policy estimand post-discontinuation values
are observed and under a hypothetical estimand they are not.

### 11.4 Analysis and inference layer

**`fit_result()`** is the contract every fitter must satisfy:
`estimate`, `se`, `p_value`, `ci_lower`, `ci_upper`, `df`,
`converged`, `level`. Limits are derived from a t reference when `df`
is finite. `null_fit()` is the non-converged sentinel.

**`fit_ancova()`** is the reference implementation: change from
baseline at a nominated visit, regressed on arm and baseline value.
The analysis layer is otherwise a deliberate non-goal; `mmrm`, `nlme`,
`lme4` and `survival` are the engines, and the contract is the seam.

**`estimand()`** declares the target quantity and its true value, so
bias and coverage have a referent that is recorded rather than passed
around implicitly.

### 11.5 Simulation driver and Monte Carlo discipline

**`sim_streams()` and `with_rng_state()`** implement Morris section
4.1. `L'Ecuyer-CMRG` substreams are used rather than per-replicate
`set.seed()`, because independence is then guaranteed by construction
and the same substreams replay across design cells, which is what
common random numbers require. Both functions restore the caller's
generator kind and state on exit, so using them has no side effect on
surrounding code.

**`run_simulation()`** drives the replicate loop and stores the RNG
state that produced each replicate. Errors in generation or analysis
are trapped: the replicate is recorded as non-converged, the message
is retained in an `errors` attribute, and the run continues. Passing a
named list of fitters analyses the identical data set with each, which
is the paired comparison Morris recommends; verified that two fitters
sharing data produce positively correlated estimates.

Verified: the study is reproducible from its master seed, and any
individual replicate can be reproduced exactly from its stored stream.

### 11.6 Performance measures

**`compute_performance()`** implements Morris Table 6 with a Monte
Carlo standard error on every measure: `n_converged`, `conv_rate`,
`bias`, `emp_se`, `mod_se`, `rel_error_mod_se`, `mse`, `coverage`,
`bias_elim_coverage`, `rejection`, `mean_ci_width`.

`rel_error_mod_se` is the diagnostic of variance estimation, and
`bias_elim_coverage` separates interval-width problems from bias.
Verified that `MCSE(bias) = emp_se / sqrt(B)` and that
`MCSE(coverage)` takes the binomial form, and that the reference
fitter is unbiased with nominal coverage and nominal type I error
under a null data-generating mechanism.

**`nsim_for_mcse()`** inverts the MCSE formula so the number of
replicates can be justified rather than conventional, per Morris
section 5.3.

### 11.7 Power and sample size

**`sim_power()`** returns simulated power with its MCSE at one sample
size, optionally applying a dropout specification to each replicate.
**`power_curve()`** sweeps a grid. **`sample_size()`** inverts a
monotone interpolation of the curve and then **confirms** the selected
size with a longer independent run, reporting the achieved power and
its MCSE. Where the target lies outside the range of the curve the
function says so explicitly rather than extrapolating.

## 12. Sample size: `zzrctsim` contrasted with `simr`

`simr` is the established simulation-based power package for mixed
models and is the natural comparator. The contrast is instructive
because the two tools solve adjacent but different problems.

### Workflow

In `simr`, sample size is obtained indirectly. One begins from a
fitted `merMod`, or constructs one without data via
`makeLmer(formula, fixef, VarCorr, sigma, data)`, extends the design
with `extend(model, along = "subject", n = ...)`, evaluates
`powerCurve(fit, along = "subject", breaks = ...)`, and reads off the
break at which power crosses the target. There is no inverse function;
the analyst interpolates by eye or by hand.

In `zzrctsim`, `sample_size(target = 0.80, n_grid = ...)` performs the
inversion and then confirms the answer with an independent longer run
at the selected size, returning the achieved power with its Monte
Carlo standard error. Measured:

```
<sample_size> target power 0.8
  n per arm: 104  (total 208 for two arms)
  confirmed power: 0.8187  (MCSE 0.00995)

 n_per_arm method  power        mcse
        40 ancova 0.4725 0.024962159
        70 ancova 0.6900 0.023124662
       110 ancova 0.8200 0.019209373
       160 ancova 0.9625 0.009499178
```

### What each can express

| | `simr` | `zzrctsim` |
|---|---|---|
| Starting point | fitted `merMod`, or `makeLmer` | design plus DGM specification |
| Covariance | `VarCorr`, so constrained to `Z G Z' + sigma^2 I` | either parameterization; arbitrary marginal `V` |
| Sample size | `powerCurve` read off by hand | `sample_size()` with inversion and confirmation |
| Dropout | must be imposed on the design *before* generation, so MCAR or covariate-dependent only | Diggle-Kenward MCAR/MAR/MNAR, per-arm, calibrated to target |
| Run-in, common close | no design vocabulary | native |
| Staggered accrual | no | native |
| Non-Gaussian | GLMMs: binomial, Poisson | planned, not yet built |
| Output | rejection rate with a binomial interval | full Morris Table 6, MCSE on every measure |
| RNG | `set.seed` | `L'Ecuyer-CMRG` substreams, per-replicate state stored |

The decisive structural difference is the one established in §3:
because `simr` generates conditionally from a random-effects
specification, the reachable marginal covariance is
`Z G Z' + sigma^2 I`. An empirically estimated covariance with
heterogeneous visit variances is not expressible at any parameter
setting. For a study whose *subject* is the covariance, that is
disqualifying; for a study where the covariance is a nuisance, it is
irrelevant.

The second difference is dropout. Because `simr` generates the
response from a fixed design matrix, missingness must be imposed
before generation and therefore cannot depend on the response. MAR
conditioning on an earlier observed value, and MNAR, are structurally
unreachable. For compendium 14, whose entire subject is dropout
adjustment, this rules `simr` out.

### Where `simr` is preferable

Honestly stated, because a comparison in which the local package wins
on every axis is not a comparison:

- **Pilot data.** `simr` takes a fitted model directly. If a pilot
  study exists, `simr` needs no covariance elicitation at all, whereas
  `zzrctsim` requires the DGM to be specified.
- **Generalized outcomes.** `simr` handles binomial and Poisson GLMMs
  today; `zzrctsim` is Gaussian-only until the endpoint families of §4
  are built.
- **Maturity.** `simr` is widely used, well documented, and its
  results are recognized by reviewers. A new package carries the
  burden of proof.
- **API.** `extend(along =, within =, n =)` is a clean abstraction for
  growing a design along any factor, and `zzrctsim` has no equivalent.

The recommendation for the portfolio is therefore not "replace `simr`
everywhere" but "use `zzrctsim` where the covariance structure, the
dropout mechanism, or the trial phase structure is the object of
study, and `simr` where a fitted pilot model is the natural starting
point and the outcome is generalized".

## 13. Migration plan for `~/prj/res`

### 13.1 Smoke test

Before proposing any migration, the engine was validated against a
**closed-form** answer rather than against another simulator, since
agreement between two simulators is weaker evidence than agreement
with an analytic result.

Design: two arms, four post-baseline visits at three-month spacing,
random intercept and slope with residual variance, ANCOVA of change
from baseline at the final visit adjusting for baseline. The exact
power follows from the induced marginal covariance:

```
  var_adj = a' V a  -  (a' V e_base)^2 / V_base,   a = y_last - y_base
  power   = Phi( |delta| / sqrt(2 var_adj / n)  -  z_{1-alpha/2} )
```

Result, `B = 2000` per row:

```
n/arm        closed  simulated       mcse        z
60           0.6058     0.5925     0.0110    -1.21
120          0.8832     0.8810     0.0072    -0.31
200          0.9825     0.9785     0.0032    -1.23
320          0.9993     0.9995     0.0005     0.44
```

Agreement is within Monte Carlo error at every point, across a power
range of 0.59 to 0.9995. The engine is sound for the Gaussian
longitudinal case. **This validates the generator and the ANCOVA
fitter only**; it says nothing about the endpoint families of §4,
which are not built, nor about the accrual and reference-based layers,
which have no closed-form comparator and are covered by unit tests
instead.

### 13.2 Migration inventory

A scan of `~/prj/res` on 2026-08-10 found 35 compendia, of which 21
are trial-simulation dependent. Existing simulation code is
concentrated in eight:

| Compendium | `R/` files | Status |
|---|---|---|
| 02-adaptive-alloc | 8 | substantial; survival plus MMRM arms |
| 14-dropout-adjustment | 2 | substantial; archive tree also present |
| 16-factorial-matched-random | 4 | substantial |
| 19-multiple-comparisons | 3 | moderate |
| 28-trialsim | 8 | the package being absorbed |
| 30-missforest | 12 | substantial |
| 31-long-power | 14 | closed-form focus, light simulation |
| 36-pmsimstats-ng | 10 | substantial; N-of-1 designs |

The remaining thirteen have an empty `R/`: the manuscripts are
scaffolded and the simulation is not yet written. **This is the
central fact of the migration.** Most of the work is not replacement
but first provision, and the cost of acting now is far lower than the
cost of reconciling thirteen divergent implementations later.

### 13.3 Sequencing

**Phase 0 — do not migrate anything yet.** `zzrctsim` is Gaussian-only
and has no fitters beyond ANCOVA. Migrating a survival or count
compendium today is not possible. Phase 0 is completing the endpoint
families of §4 and at least an MMRM fitter.

**Phase 1 — greenfield adoption, no migration risk.** Use `zzrctsim`
for compendia whose `R/` is empty and whose needs are Gaussian
longitudinal: 05, 07, 08, 09, 11, 12, 18, 21. Nothing is replaced, so
nothing can regress. This also exercises the API against eight
independent studies before any working code depends on it.

**Phase 2 — absorb 28-trialsim.** The package being absorbed by
decision. Port `hc_change_sds()`, `hc_baseline_change_cor()` and
`joint_covariance()` as the Alzheimer's calibration layer. Retire the
`# Source:` headers citing the vendored Pfizer-origin code, per the
white paper. Verify by reproducing the three RDS outputs in
`analysis/data/derived_data/` to within Monte Carlo error.

**Phase 3 — migrate compendia with existing code, one at a time,
each behind an equivalence test.** Order by tractability: 19
(multiplicity, Gaussian), then 14 (dropout; the missingness layer is
built for it), then 30 (imputation; needs `mask_from_data()`), then 16
(factorial; needs a randomization layer), then 36 (N-of-1), then 02
(needs survival).

**The equivalence test is the gate.** For each compendium, run the
existing code and the `zzrctsim` implementation under matched
parameters and compare headline performance measures with their
MCSEs. A difference exceeding roughly three combined MCSEs blocks the
swap and must be explained before proceeding. Discrepancies are as
likely to reveal a defect in the original as in the replacement; the
covariance defect found in the vendored Rogers code is precedent.

### 13.4 Risks

- **Single point of failure.** Twenty-one manuscripts depending on one
  package means a defect propagates everywhere. Mitigation: the
  verification suite of `res/28-trialsim/inst/tinytest/test_dgm_verification.R`
  should be ported and run against `zzrctsim`, and the closed-form
  smoke test above should become a permanent regression test.
- **API churn.** Compendia adopting early will be broken by later
  design changes. Mitigation: freeze the design and DGM layers before
  Phase 1, and version thereafter.
- **Reviewer unfamiliarity.** A manuscript resting on an unpublished
  local package invites the question of whether the engine is correct.
  Mitigation: the closed-form validation of §13.1 belongs in the
  supplementary material of every compendium that uses it.
- **Loss of provenance.** Migrating replaces code that produced
  already-reported results. No published or submitted analysis should
  be migrated; migration applies to work in progress only.
