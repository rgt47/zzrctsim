# Performance measures with Monte Carlo standard errors.
#
# Formulas follow Morris, White and Crowther (2019) Table 6. Every
# measure is reported with its MCSE, because a performance measure
# without one cannot be distinguished from Monte Carlo noise, and
# reporting it is what justifies the choice of B.

#' Performance measures for a simulation study
#'
#' @param results A `sim_results` object from [run_simulation()], or a
#'   data frame with `estimate`, `se`, `ci_lower`, `ci_upper`,
#'   `p_value` and `converged`.
#' @param theta Numeric. True value of the estimand. Taken from the
#'   `estimand` attribute when absent.
#' @param alpha Numeric. Significance level for the rejection rate.
#' @return A data frame with columns `method`, `measure`, `estimate`
#'   and `mcse`, and one row per method and measure: eleven rows per
#'   method, or a single `n_converged` row for any method with fewer
#'   than two converged replicates. `mcse` is `NA` for the two measures
#'   that are counts rather than estimates.
#' @details
#' Measures, in the order emitted, with `B` the number of converged
#' replicates:
#' \describe{
#'   \item{`n_converged`}{count of usable replicates: those flagged
#'     converged with a non-missing estimate. No MCSE.}
#'   \item{`conv_rate`}{`n_converged` divided by the number of
#'     replicates attempted for that method. No MCSE. A rate below one
#'     means the remaining measures are conditional on convergence and
#'     may be selectively biased.}
#'   \item{`bias`}{`mean(est) - theta`; MCSE is the standard error of
#'     the mean}
#'   \item{`emp_se`}{`sd(est)`; MCSE `emp_se / sqrt(2(B-1))`}
#'   \item{`mod_se`}{`sqrt(mean(se^2))`, the average model-based
#'     standard error}
#'   \item{`rel_error_mod_se`}{percentage by which `mod_se` over- or
#'     under-states `emp_se`; a key diagnostic of variance estimation}
#'   \item{`mse`}{`mean((est - theta)^2)`}
#'   \item{`coverage`}{proportion of intervals containing `theta`}
#'   \item{`bias_elim_coverage`}{coverage of intervals recentered on
#'     `mean(est)`, which separates interval width problems from bias}
#'   \item{`rejection`}{proportion with `p < alpha`; power under an
#'     alternative, type I error under the null}
#'   \item{`mean_ci_width`}{`mean(ci_upper - ci_lower)`, with MCSE the
#'     standard error of that mean. Read alongside `coverage`: an
#'     interval can reach nominal coverage by being needlessly wide,
#'     and the pair distinguishes that from a well-calibrated one.}
#' }
#' @export
compute_performance <- function(results, theta = NULL, alpha = 0.05) {
  if (is.null(theta)) {
    es <- attr(results, "estimand")
    if (is.null(es)) {
      stop("`theta` is required when `results` carries no estimand.")
    }
    theta <- es$value
  }
  meths <- unique(results$method)
  out <- lapply(meths, function(m) {
    z <- results[results$method == m & results$converged &
                   !is.na(results$estimate), , drop = FALSE]
    .perf_one(z, theta, alpha, nrow(results[results$method == m, ]))
  })
  res <- do.call(rbind, Map(function(m, d) cbind(method = m, d),
                            meths, out))
  rownames(res) <- NULL
  res
}

.perf_one <- function(z, theta, alpha, n_total) {
  B <- nrow(z)
  add <- function(measure, estimate, mcse) {
    data.frame(measure = measure, estimate = estimate, mcse = mcse,
               stringsAsFactors = FALSE)
  }
  if (B < 2L) {
    return(add("n_converged", B, NA_real_))
  }

  est <- z$estimate
  se <- z$se
  m_est <- mean(est)
  v_est <- stats::var(est)

  bias <- m_est - theta
  bias_mcse <- sqrt(v_est / B)

  emp <- sqrt(v_est)
  emp_mcse <- emp / sqrt(2 * (B - 1))

  ok_se <- !is.na(se)
  mod <- sqrt(mean(se[ok_se]^2))
  mod_mcse <- if (sum(ok_se) > 1L) {
    sqrt(stats::var(se[ok_se]^2) / (4 * sum(ok_se) * mod^2))
  } else NA_real_

  rel <- 100 * (mod / emp - 1)
  rel_mcse <- 100 * (mod / emp) *
    sqrt(mod_mcse^2 / mod^2 + 1 / (2 * (B - 1)))

  sq <- (est - theta)^2
  mse <- mean(sq)
  mse_mcse <- sqrt(stats::var(sq) / B)

  cov_ind <- z$ci_lower <= theta & z$ci_upper >= theta
  cov_ind <- cov_ind[!is.na(cov_ind)]
  cvg <- mean(cov_ind)
  cvg_mcse <- sqrt(cvg * (1 - cvg) / length(cov_ind))

  half <- (z$ci_upper - z$ci_lower) / 2
  be_ind <- (z$ci_lower - m_est + theta) <= theta &
    (z$ci_upper - m_est + theta) >= theta
  be_ind <- be_ind[!is.na(be_ind)]
  be <- mean(be_ind)
  be_mcse <- sqrt(be * (1 - be) / length(be_ind))

  rej_ind <- z$p_value < alpha
  rej_ind <- rej_ind[!is.na(rej_ind)]
  rej <- mean(rej_ind)
  rej_mcse <- sqrt(rej * (1 - rej) / length(rej_ind))

  rbind(
    add("n_converged", B, NA_real_),
    add("conv_rate", B / n_total, NA_real_),
    add("bias", bias, bias_mcse),
    add("emp_se", emp, emp_mcse),
    add("mod_se", mod, mod_mcse),
    add("rel_error_mod_se", rel, rel_mcse),
    add("mse", mse, mse_mcse),
    add("coverage", cvg, cvg_mcse),
    add("bias_elim_coverage", be, be_mcse),
    add("rejection", rej, rej_mcse),
    add("mean_ci_width", mean(2 * half, na.rm = TRUE),
        stats::sd(2 * half, na.rm = TRUE) / sqrt(B))
  )
}

#' Replicates required for a target Monte Carlo standard error
#'
#' Morris et al. section 5.3: the number of replicates should be
#' justified by the precision required of the headline performance
#' measure, not chosen by convention.
#'
#' @param measure Character. `"proportion"` for coverage, power or type
#'   I error; `"bias"` for a mean.
#' @param mcse Numeric. Target Monte Carlo standard error.
#' @param p Numeric. Anticipated proportion, for `"proportion"`. The
#'   worst case is `0.5`.
#' @param sd Numeric. Anticipated standard deviation of the estimator,
#'   for `"bias"`.
#' @return Required number of replicates.
#' @examples
#' # Type I error to within 0.005 at a nominal 0.05
#' nsim_for_mcse("proportion", mcse = 0.005, p = 0.05)
#' @export
nsim_for_mcse <- function(measure = c("proportion", "bias"),
                          mcse, p = 0.5, sd = 1) {
  measure <- match.arg(measure)
  if (!is.numeric(mcse) || length(mcse) != 1L) {
    stop("`mcse` must be a single number, but it is a ",
         class(mcse)[1], " of length ", length(mcse),
         ". Supply the target Monte Carlo standard error, ",
         "for example 0.005.",
         call. = FALSE)
  }
  if (is.na(mcse) || mcse <= 0) {
    stop("`mcse` must be a positive number; received ", mcse,
         ". It is the target Monte Carlo standard error on the ",
         "scale of the measure, so a smaller value asks for more ",
         "replicates.",
         call. = FALSE)
  }
  n <- switch(measure,
              proportion = p * (1 - p) / mcse^2,
              bias = (sd / mcse)^2)
  ceiling(n)
}
