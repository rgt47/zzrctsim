library(tinytest)

# Missingness as two composable layers: mask generation, mask
# application, and pattern-mixture post-discontinuation trajectories.

set.seed(101)
s <- trial_schedule(treatment = 5, interval = 3)
arm <- factor(rep(c("placebo", "active"), each = 250),
              levels = c("placebo", "active"))
d <- runin_design(s, arm, reference = "placebo")
gd <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
b <- c(x_slope = 0.5, x_hinge = -0.1, x_trt_active = -0.25)
out <- generate_outcomes(d, gd, beta = b, intercept = 20)
out$armf <- arm[out$id]

# ---- item 4: per-arm and covariate-dependent psi0 -------------------

m_by <- dropout_mask(out, target = c(placebo = 0.15, active = 0.35),
                     by = "armf")
sp <- attr(m_by, "spec")
expect_equal(as.numeric(sp$expected[["placebo"]]), 0.15,
             tolerance = 1e-6, info = 'per-arm target met, placebo')
expect_equal(as.numeric(sp$expected[["active"]]), 0.35,
             tolerance = 1e-6, info = 'per-arm target met, active')
expect_true(sp$psi0[["active"]] > sp$psi0[["placebo"]],
            info = 'higher target implies higher psi0')

ever <- tapply(m_by$missing, m_by$id, any)
by_arm <- tapply(ever, arm[as.integer(names(ever))], mean)
expect_true(by_arm[["active"]] > by_arm[["placebo"]],
            info = 'realized dropout is differential by arm')

# A single target still works and is unnamed-safe.
m1 <- dropout_mask(out, target = 0.20)
expect_equal(as.numeric(attr(m1, "spec")$expected), 0.20,
             tolerance = 1e-6, info = 'single target met')

# Covariate dependence.
m_cov <- dropout_mask(out, target = 0.25, psi_cov = c(x_slope = 0.02))
expect_equal(as.numeric(attr(m_cov, "spec")$expected), 0.25,
             tolerance = 1e-6,
             info = 'calibration holds with a covariate term')
expect_error(dropout_mask(out, psi_cov = c(nope = 1)),
             info = 'unknown psi_cov column is rejected')

# ---- item 6: intermittent missingness --------------------------------

set.seed(5)
m_mono <- dropout_mask(out, target = 0.30, monotone = TRUE)
m_int <- dropout_mask(out, target = 0.30, monotone = FALSE)

is_mono <- function(mk) {
  all(vapply(split(mk, mk$id), function(z) {
    v <- z$missing[order(z$time)]
    !any(diff(v) < 0)
  }, logical(1)))
}
expect_true(is_mono(m_mono), info = 'monotone mask is absorbing')
expect_false(is_mono(m_int),
             info = 'intermittent mask has interior gaps')
expect_true(sum(m_int$missing) < sum(m_mono$missing),
            info = 'intermittent loses fewer observations than dropout')
expect_equal(as.numeric(attr(m_int, "spec")$expected), 0.30,
             tolerance = 1e-6,
             info = 'intermittent calibrates to the same target')

# ---- item 5: mask application and import -----------------------------

masked <- apply_mask(out, m_mono)
expect_equal(sum(is.na(masked$y)), sum(m_mono$missing),
             info = 'apply_mask hides exactly the masked rows')
expect_true(all(!is.na(masked$y[masked$time <= 0])),
            info = 'baseline retained by default `from`')

# A mask imported from observed data reproduces itself.
m_back <- mask_from_data(masked)
expect_equal(m_back$missing, m_mono$missing,
             info = 'mask_from_data round-trips')
expect_true(attr(m_back, "spec")$monotone,
            info = 'monotonicity detected on import')
expect_equal(attr(m_back, "spec")$mechanism, "empirical",
             info = 'imported masks are labelled empirical')

# The same mask applies to an independent replicate of the same design.
set.seed(202)
out2 <- generate_outcomes(d, gd, beta = b, intercept = 20)
masked2 <- apply_mask(out2, m_mono)
expect_equal(is.na(masked2$y), is.na(masked$y),
             info = 'a mask is reusable across replicates (CRN-safe)')

expect_error(apply_mask(out[1:10, ], m_mono[1:3, ]),
             info = 'incomplete mask coverage is an error')

# ---- item 7: reference-based trajectories ----------------------------

set.seed(303)
m_ref <- dropout_mask(out, target = 0.40, by = "armf",
                      monotone = TRUE)

j2r <- reference_based(out, m_ref, method = "j2r")
cir <- reference_based(out, m_ref, method = "cir")
lmcf <- reference_based(out, m_ref, method = "lmcf")
cr <- reference_based(out, m_ref, method = "cr")

expect_true("discontinued" %in% names(j2r),
            info = 'discontinuation flag added')
expect_equal(sum(j2r$discontinued), sum(m_ref$missing),
            info = 'discontinuation starts at the first masked visit')

# Pre-discontinuation values are untouched by every method except CR,
# which copies the reference throughout by construction.
pre <- !j2r$discontinued
expect_equal(j2r$y[pre], out$y[pre],
             info = 'J2R leaves pre-discontinuation data unchanged')
expect_equal(cir$y[pre], out$y[pre],
             info = 'CIR leaves pre-discontinuation data unchanged')

# The residual is preserved, so within-subject correlation survives.
expect_equal(j2r$y - j2r$mu, out$y - out$mu,
             info = 'J2R preserves realized residuals')

# The methods are genuinely different from each other and from MAR.
post <- j2r$discontinued & (arm[j2r$id] == "active")
expect_true(any(abs(j2r$mu[post] - out$mu[post]) > 1e-8),
            info = 'J2R shifts the treated post-withdrawal mean')
expect_false(isTRUE(all.equal(j2r$mu[post], cir$mu[post])),
             info = 'J2R and CIR differ')
expect_false(isTRUE(all.equal(cir$mu[post], lmcf$mu[post])),
             info = 'CIR and LMCF differ')

# Ordering: on a worsening scale the treated arm does better than
# reference, so jumping to reference must raise the mean, and copying
# increments must raise it by less than jumping.
d_j2r <- mean(j2r$mu[post] - out$mu[post])
d_cir <- mean(cir$mu[post] - out$mu[post])
expect_true(d_j2r > 0,
            info = 'J2R moves the treated mean toward the control arm')
expect_true(d_cir > 0 && d_cir < d_j2r,
            info = 'CIR retains accumulated benefit, so shifts less than J2R')

# CR modifies the whole trajectory for discontinuing subjects.
disc_ids <- unique(j2r$id[j2r$discontinued])
k <- cr$id %in% disc_ids & !cr$discontinued & arm[cr$id] == "active"
expect_true(any(abs(cr$mu[k] - out$mu[k]) > 1e-8),
            info = 'CR alters pre-discontinuation means too')

expect_error(reference_based(out, m_ref, method = "j2r",
                             trt_cols = character(0)),
             info = 'missing treatment columns is an error')
