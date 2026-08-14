library(tinytest)

# Dual parameterization, the reachability certificate, and generation.

p <- 6L
tt <- seq_len(p)
CS <- cov_cs(sd = 1, rho = 0.5)(tt)
AR1 <- cov_ar1(sd = 1, rho = 0.6)(tt)

# ---- the certificate ------------------------------------------------

cs_cert <- certify(CS)
expect_equal(cs_cert$q_min, 1L,
             info = 'compound symmetry needs exactly one random effect')
expect_equal(cs_cert$sigma2, 0.5, tolerance = 1e-10,
             info = 'CS residual variance is 1 - rho')

ar_cert <- certify(AR1)
expect_equal(ar_cert$q_min, 5L,
             info = 'AR(1) at p = 6 needs p - 1 random effects')

# The construction is exact, not approximate.
expect_equal(reconstruct(cs_cert), CS, tolerance = 1e-10,
             info = 'certificate reconstructs CS exactly')
expect_equal(reconstruct(ar_cert), AR1, tolerance = 1e-10,
             info = 'certificate reconstructs AR(1) exactly')

# CS has a flat leading eigenvector; AR(1) does not. This is why a
# random intercept reproduces CS and cannot reproduce AR(1).
expect_true(diff(range(abs(cs_cert$Z[, 1]))) < 1e-8,
            info = 'CS leading eigenvector is flat (all-ones direction)')
expect_true(diff(range(abs(ar_cert$Z[, 1]))) > 0.05,
            info = 'AR(1) leading eigenvector is bowed')

expect_error(certify(matrix(c(1, 2, 3, 4), 2, 2)),
             info = 'asymmetric input is rejected')

# Degradation is monotone in the budget and hits zero at q_min.
err <- reach_error(AR1)
expect_true(all(diff(err$rel_frobenius) < 0),
            info = 'error decreases as the budget grows')
expect_true(err$rel_frobenius[err$q == 5L] < 1e-8,
            info = 'error vanishes at q = q_min')
expect_true(err$rel_frobenius[err$q == 1L] > 0.25,
            info = 'a single random effect is badly wrong for AR(1)')

# ---- dual parameterization round trips ------------------------------

# Conditional -> marginal is always exact.
G <- matrix(c(4, 0.5, 0.5, 0.25), 2, 2)
dc <- dgm_conditional(G = G, sigma2 = 1.5)
Vc <- cov_at(dc, tt)
expect_true(isSymmetric(Vc), info = 'induced marginal is symmetric')
dm <- as_marginal(dc, tt)
expect_equal(cov_at(dm, tt), Vc, tolerance = 1e-12,
             info = 'as_marginal is exact')

# Marginal -> conditional succeeds at q_min and refuses below it.
dm_cs <- dgm_marginal(cov_cs(1, 0.5), times = NULL)
dc_cs <- as_conditional(dm_cs, tt)
expect_equal(cov_at(dc_cs, tt), CS, tolerance = 1e-10,
             info = 'CS round-trips through the conditional form')

dm_ar <- dgm_marginal(cov_ar1(1, 0.6))
expect_error(as_conditional(dm_ar, tt, q = 2),
             info = 'conversion refuses an insufficient budget')
dc_ar <- as_conditional(dm_ar, tt)
expect_equal(cov_at(dc_ar, tt), AR1, tolerance = 1e-10,
             info = 'AR(1) round-trips at q_min = 5')

# A random intercept induces compound symmetry and nothing else.
ri <- dgm_conditional(G = matrix(4), sigma2 = 2,
                      z = function(t) matrix(1, length(t), 1))
Vri <- cov_at(ri, tt)
expect_equal(certify(Vri)$q_min, 1L,
             info = 'random intercept certifies at q = 1')
off <- Vri[upper.tri(Vri)]
expect_true(diff(range(off)) < 1e-10,
            info = 'random intercept gives constant off-diagonal (CS)')

# ---- covariance functions extend to irregular grids ------------------

t_irr <- c(0, 1, 2.5, 7, 11)
Vi <- cov_ar1(sd = 2, rho = 0.6)(t_irr)
expect_equal(dim(Vi), c(5L, 5L), info = 'AR(1) evaluates on any grid')
expect_equal(Vi[1, 2], 4 * 0.6^1, tolerance = 1e-12,
             info = 'continuous-time AR(1) uses |t_i - t_j|')
expect_equal(Vi[1, 4], 4 * 0.6^7, tolerance = 1e-12,
             info = 'unequal spacing handled correctly')

# A matrix-valued marginal DGM cannot leave its grid, and says so.
dm_fixed <- dgm_marginal(AR1, times = tt)
expect_equal(cov_at(dm_fixed, c(2, 4)), AR1[c(2, 4), c(2, 4)],
             info = 'submatrix taken for a subset of the grid')
expect_error(cov_at(dm_fixed, c(2, 99)),
             info = 'off-grid evaluation is refused, not extrapolated')

# ---- generation ------------------------------------------------------

set.seed(42)
s <- trial_schedule(run_in = 1, treatment = 4, interval = 3)
arm <- factor(rep(c("placebo", "active"), each = 400),
              levels = c("placebo", "active"))
d <- runin_design(s, arm, reference = "placebo")
gd <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
b <- c(x_slope = 0.5, x_trt_active = -0.25)
out <- generate_outcomes(d, gd, beta = b, intercept = 20)

expect_true(all(c("mu", "y") %in% names(out)),
            info = 'generate_outcomes adds mu and y')
expect_equal(nrow(out), nrow(d), info = 'no rows added or lost')
expect_true(all(!is.na(out$y)), info = 'complete data before dropout')

# The fixed-effect mean is recovered.
fit <- stats::lm(y ~ x_slope + x_trt_active, data = out)
est <- stats::coef(fit)[names(b)]
se <- summary(fit)$coefficients[names(b), "Std. Error"]
expect_true(all(abs(est - b) < 4 * se),
            info = 'fixed effects recovered within 4 SE')

expect_error(generate_outcomes(d, gd, beta = c(nope = 1)),
             info = 'unknown beta names are rejected')

# ---- dropout ---------------------------------------------------------

set.seed(7)
dd <- apply_dropout(out, target = 0.30, from = 0)
expect_true(any(is.na(dd$y)), info = 'dropout removes observations')
expect_true(all(!is.na(dd$y[dd$time <= 0])),
            info = 'run-in and baseline visits are retained')

# Mechanisms are nested restrictions of one model.
expect_equal(dropout_mechanism(0, 0), "MCAR", info = 'psi1 = psi2 = 0')
expect_equal(dropout_mechanism(0.2, 0), "MAR", info = 'psi2 = 0')
expect_equal(dropout_mechanism(0.2, 0.1), "MNAR", info = 'psi2 free')
expect_equal(dropout_mechanism(0, 0.1), "MNAR",
             info = 'psi2 alone is still MNAR')

# Calibration: the realized proportion must hit the target, and must
# keep hitting it as psi1 and psi2 vary. This is the property the
# previous `rate` argument did not have.
for (spec in list(c(0, 0), c(0.15, 0), c(0.15, 0.10), c(0, 0.10))) {
  for (tg in c(0.15, 0.35)) {
    set.seed(21)
    z <- apply_dropout(out, target = tg,
                       psi1 = spec[1], psi2 = spec[2], from = 0)
    sp <- attr(z, "dropout_spec")
    expect_equal(sp$expected, tg, tolerance = 1e-6,
                 info = paste0('expected dropout matches target ', tg,
                               ' at psi1=', spec[1], ' psi2=', spec[2]))
    realized <- mean(tapply(z$y, z$id, function(v) any(is.na(v))))
    expect_true(abs(realized - tg) < 0.06,
                info = paste0('realized dropout near target ', tg,
                              ' at psi1=', spec[1], ' psi2=', spec[2]))
  }
}

# Centring is what makes calibration stable on an untransformed scale.
sp <- attr(apply_dropout(out, target = 0.25, psi1 = 0.15), "dropout_spec")
expect_equal(sp$center, mean(out$y), tolerance = 1e-10,
             info = 'default centring is the complete-data mean')

# Monotone: once dropped, dropped thereafter.
mono <- tapply(seq_len(nrow(dd)), dd$id, function(k) {
  k <- k[order(dd$time[k])]
  m <- is.na(dd$y[k])
  !any(diff(m) < 0)
})
expect_true(all(unlist(mono)), info = 'dropout is monotone')

# Mechanism-label tests. The hazard reads the COMPLETE-data previous
# response, so the check must regress realized dropout on `out$y`, not
# on the post-dropout column and not on the deterministic mean.
risk_set <- function(complete, realized, from = 0) {
  do.call(rbind, lapply(unique(complete$id), function(i) {
    k <- which(complete$id == i)
    k <- k[order(complete$time[k])]
    y_c <- complete$y[k]
    na <- is.na(realized$y[k])
    tv <- complete$time[k]
    m <- seq_len(length(k) - 1L)
    keep <- tv[m + 1L] > from & !na[m]
    data.frame(prev = y_c[m][keep],
               drop_next = (na[m + 1L] & !na[m])[keep])
  }))
}

set.seed(11)
mar <- apply_dropout(out, target = 0.30, psi1 = 0.20, from = 0)
pm <- risk_set(out, mar)
fit_mar <- stats::glm(drop_next ~ prev, data = pm,
                      family = stats::binomial())
expect_true(stats::coef(fit_mar)[2] > 0,
            info = 'MAR dropout increases with the prior response')
# and the association must be of the right size, not merely positive:
# the previous check passed even when 100% of subjects dropped out.
expect_true(mean(tapply(mar$y, mar$id, function(v) any(is.na(v)))) < 0.45,
            info = 'MAR dropout stays near its target, not saturated')

set.seed(11)
mc <- apply_dropout(out, target = 0.30, from = 0)
pc <- risk_set(out, mc)
fit_mc <- stats::glm(drop_next ~ prev, data = pc,
                     family = stats::binomial())
ci <- suppressMessages(stats::confint(fit_mc, "prev", level = 0.999))
expect_true(ci[1] <= 0 && ci[2] >= 0,
            info = 'MCAR dropout does not depend on the prior response')

# MNAR nests MAR: with psi1 held fixed, adding psi2 must not remove the
# history dependence. This is the defect the old parameterization had.
set.seed(11)
mnar <- apply_dropout(out, target = 0.30, psi1 = 0.20, psi2 = 0.10,
                      from = 0)
sp_mnar <- attr(mnar, "dropout_spec")
expect_equal(sp_mnar$psi1, 0.20,
             info = 'MNAR retains the psi1 history term')
expect_equal(sp_mnar$mechanism, "MNAR", info = 'psi2 makes it MNAR')

# ---- family / link on the conditional generator ----------------------

set.seed(4242)
s2 <- trial_schedule(treatment = 4, interval = 3)
arm2 <- factor(rep(c("placebo", "active"), each = 500),
               levels = c("placebo", "active"))
d2 <- runin_design(s2, arm2, reference = "placebo")
bb <- c(x_slope = 0.02, x_trt_active = -0.04)

# Binomial with a logit link.
gb <- dgm_conditional(G = matrix(0.5), sigma2 = NULL,
                      z = function(t) matrix(1, length(t), 1),
                      family = stats::binomial())
ob <- generate_outcomes(d2, gb, beta = bb, intercept = 0)
expect_true(all(ob$y %in% c(0, 1)),
            info = 'binomial family yields a binary response')
expect_true(all(ob$cmean > 0 & ob$cmean < 1),
            info = 'conditional mean is a probability')
expect_equal(ob$cmean, stats::plogis(ob$eta), tolerance = 1e-12,
             info = 'cmean is the inverse logit of eta')
expect_true(abs(mean(ob$y) - mean(ob$cmean)) < 0.02,
            info = 'realized rate tracks the conditional mean')
expect_equal(dim(attr(ob, "ranef")), c(1000L, 1L),
             info = 'random effects retained, one row per subject')

# Poisson with a log link.
gp <- dgm_conditional(G = matrix(0.2), sigma2 = NULL,
                      z = function(t) matrix(1, length(t), 1),
                      family = stats::poisson())
op <- generate_outcomes(d2, gp, beta = bb, intercept = 1)
expect_true(all(op$y >= 0 & op$y == floor(op$y)),
            info = 'poisson family yields non-negative counts')
expect_equal(op$cmean, exp(op$eta), tolerance = 1e-10,
             info = 'cmean is exp(eta) under the log link')

# A probit link is accepted and differs from logit.
gpr <- dgm_conditional(G = matrix(0.5), sigma2 = NULL,
                       z = function(t) matrix(1, length(t), 1),
                       family = stats::binomial("probit"))
opr <- generate_outcomes(d2, gpr, beta = bb, intercept = 0)
expect_equal(opr$cmean, stats::pnorm(opr$eta), tolerance = 1e-12,
             info = 'probit link applied')

# Dependence is induced by the shared random effect: repeated measures
# on a subject must be positively associated.
wide <- matrix(ob$y, ncol = nrow(s2), byrow = TRUE)
expect_true(stats::cor(wide[, 2], wide[, 3]) > 0.05,
            info = 'shared random intercept induces within-subject dependence')

# ---- guards on Gaussian-only machinery -------------------------------

expect_error(cov_at(gb, s2$time),
             info = 'cov_at refuses a non-Gaussian family')
expect_error(dgm_conditional(matrix(1), family = stats::gaussian()),
             info = 'sigma2 is required for the Gaussian family')
expect_error(dgm_conditional(matrix(1), sigma2 = 1,
                             family = stats::Gamma()),
             info = 'unsupported families are rejected')

# linpred_cov is defined for every family and equals Z G Z'.
lp <- linpred_cov(gb, s2$time)
expect_equal(dim(lp), rep(nrow(s2), 2L),
             info = 'linpred_cov returns a covariance on the eta scale')
expect_true(all(abs(lp - 0.5) < 1e-12),
            info = 'random intercept gives constant eta covariance')

# Gaussian behaviour is unchanged by the new argument.
set.seed(11)
g_old <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
v1 <- cov_at(g_old, s2$time)
expect_equal(v1, zzrctsim:::.zgz(g_old, s2$time) + 4 * diag(nrow(s2)),
             info = 'Gaussian cov_at is still Z G Z\' + sigma^2 I')
