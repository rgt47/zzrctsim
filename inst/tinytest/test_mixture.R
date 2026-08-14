library(tinytest)

# Finite mixtures of multivariate normals. A component is drawn per
# subject and held across that subject's visits, so the mixture is a
# statement about latent subpopulations rather than visit-level noise.

tt <- 1:5
slow <- dgm_conditional(G = diag(c(4, 0.01)), sigma2 = 2)
fast <- dgm_conditional(G = diag(c(16, 0.09)), sigma2 = 2)

# ---- construction ---------------------------------------------------

mx <- dgm_mixture(list(slow, fast), weights = c(0.7, 0.3))
expect_inherits(mx, "dgm_mixture", info = 'constructor returns a dgm')
expect_inherits(mx, "dgm", info = 'a mixture is a dgm')
expect_equal(mx$K, 2L, info = 'component count recorded')
expect_equal(sum(mx$weights), 1, info = 'weights normalized')

expect_error(dgm_mixture(list(slow), weights = 1),
             info = 'a mixture needs at least two components')
expect_error(dgm_mixture(list(slow, fast), weights = c(0.5, 0.7)),
             info = 'weights must sum to one')
expect_error(dgm_mixture(list(slow, "not a dgm")),
             info = 'components must be dgm objects')
expect_error(dgm_mixture(list(slow, dgm_tte())),
             info = 'a time-to-event component cannot be mixed')

# Default weights are equal; numeric mean_shift is accepted.
expect_equal(dgm_mixture(list(slow, fast))$weights, c(0.5, 0.5),
             info = 'default weights are equal')
mx_num <- dgm_mixture(list(slow, fast), mean_shift = c(0, 3))
expect_equal(mx_num$mean_shift[[2]](tt), rep(3, length(tt)),
             info = 'numeric mean_shift becomes a constant function')

# ---- the marginal covariance is the within-plus-between decomposition

mx2 <- dgm_mixture(list(slow, fast), weights = c(0.6, 0.4),
                   mean_shift = list(function(t) 0 * t,
                                     function(t) 0.5 * t))
Vmix <- cov_at(mx2, tt)
w <- c(0.6, 0.4)
shifts <- list(0 * tt, 0.5 * tt)
mbar <- w[1] * shifts[[1]] + w[2] * shifts[[2]]
within <- w[1] * cov_at(slow, tt) + w[2] * cov_at(fast, tt)
between <- w[1] * outer(shifts[[1]] - mbar, shifts[[1]] - mbar) +
  w[2] * outer(shifts[[2]] - mbar, shifts[[2]] - mbar)
expect_equal(Vmix, within + between, tolerance = 1e-12,
             info = 'cov_at is exactly within + between')
expect_true(isSymmetric(Vmix), info = 'mixture covariance is symmetric')
expect_true(all(eigen(Vmix, only.values = TRUE)$values > 0),
            info = 'mixture covariance is positive definite')

# A mean shift that differs across components inflates the covariance
# relative to the within-component average. That is the between term.
expect_true(all(diag(Vmix) >= diag(within) - 1e-10),
            info = 'between-component spread adds variance')

# With identical shifts the between term vanishes.
mx_same <- dgm_mixture(list(slow, fast), weights = c(0.6, 0.4))
expect_equal(cov_at(mx_same, tt), within, tolerance = 1e-12,
             info = 'equal means give a pure within-component mixture')

# ---- generation ------------------------------------------------------

set.seed(6001)
s <- trial_schedule(treatment = 4, interval = 3)
arm <- factor(rep(c("placebo", "active"), each = 1500),
              levels = c("placebo", "active"))
d <- runin_design(s, arm, reference = "placebo")
b <- c(x_slope = 0.4, x_trt_active = -0.2)

gm <- dgm_mixture(list(slow, fast), weights = c(0.7, 0.3),
                  mean_shift = list(function(t) 0 * t,
                                    function(t) 0.6 * t))
out <- generate_outcomes(d, gm, beta = b, intercept = 12)

expect_true(all(c("component", "y", "mu") %in% names(out)),
            info = 'component label returned alongside the response')
expect_equal(nrow(out), nrow(d), info = 'no rows added or lost')
expect_true(all(!is.na(out$y)), info = 'complete data')

# The component is constant within a subject: this is what makes it a
# subpopulation rather than visit-level contamination.
per_subject <- tapply(out$component, out$id, function(v) length(unique(v)))
expect_true(all(unlist(per_subject) == 1L),
            info = 'component is fixed within a subject')

# Realized mixing proportions match the weights.
comp1 <- attr(out, "component")
expect_equal(mean(comp1 == 1L), 0.7, tolerance = 0.03,
             info = 'realized weights match the specification')

# The empirical covariance matches cov_at() -- the decomposition is
# the true covariance of the generated data, not a nominal quantity.
Y <- matrix(out$y, ncol = nrow(s), byrow = TRUE)
emp <- stats::cov(Y)
target <- cov_at(gm, s$time)
expect_true(max(abs(emp - target) / target) < 0.12,
            info = 'empirical covariance matches the mixture covariance')

# Components differ in trajectory: the shifted component ends higher.
last <- out[out$time == max(out$time), ]
expect_true(mean(last$y[last$component == 2L]) >
              mean(last$y[last$component == 1L]),
            info = 'the shifted component progresses further')

# The fixed effects are still recoverable, and the treatment contrast
# is unaffected by the mixture, because components are assigned
# independently of arm.
fit <- stats::lm(y ~ x_slope + x_trt_active, data = out)
est <- stats::coef(fit)[names(b)]
se <- summary(fit)$coefficients[names(b), "Std. Error"]
expect_true(abs(est[["x_trt_active"]] - b[["x_trt_active"]]) <
              4 * se[["x_trt_active"]],
            info = 'treatment effect recovered under a mixture')

# ---- the mixture composes with the rest of the machinery ------------

# Missingness applies unchanged, because a normal mixture has moments.
mk <- dropout_mask(out, target = 0.25)
expect_equal(as.numeric(attr(mk, "spec")$expected), 0.25,
             tolerance = 1e-6,
             info = 'dropout calibration works on mixture data')

msk <- apply_mask(out, mk)
expect_true(any(is.na(msk$y)), info = 'mask applies to mixture data')

# And the driver runs it end to end.
set.seed(6002)
res <- run_simulation(
  B = 40,
  generate = function() generate_outcomes(d, gm, beta = b,
                                          intercept = 12),
  analyze = list(ancova = function(z)
    fit_ancova(z, reference = "placebo")),
  estimand = estimand("final contrast",
                      b[["x_trt_active"]] * max(s$time)),
  seed = 3L)
perf <- compute_performance(res)
bias <- perf$estimate[perf$measure == "bias"]
bias_mcse <- perf$mcse[perf$measure == "bias"]
expect_true(abs(bias) < 4 * bias_mcse,
            info = 'ANCOVA remains unbiased under a normal mixture')
