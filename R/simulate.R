# Estimands, the fitter contract, and the replicate driver.
#
# ADEMP separates the estimand (what is being estimated), the methods
# (how), and the performance measures (how well). These are kept as
# distinct objects so that a simulation can be re-summarized under a
# different estimand or a new performance measure without re-running.

#' Declare an estimand
#'
#' @param name Character. Label, used in output.
#' @param value Numeric. The true value of the target quantity under
#'   the data-generating mechanism.
#' @param description Character. Optional prose definition. Recording
#'   this matters: an estimand is a statement about a population
#'   quantity, not about an estimator.
#' @return An object of class `estimand`, a list with
#'   \describe{
#'     \item{name}{the label, as supplied}
#'     \item{value}{the true value; [compute_performance()] takes
#'       `theta` from here when none is passed}
#'     \item{description}{the prose definition, or `NULL`}
#'   }
#' @export
estimand <- function(name, value, description = NULL) {
  stopifnot(is.character(name), length(name) == 1L,
            is.numeric(value), length(value) == 1L)
  structure(list(name = name, value = value,
                 description = description),
            class = "estimand")
}

#' @export
print.estimand <- function(x, ...) {
  cat("<estimand>", x$name, "=", format(x$value, digits = 6), "\n")
  if (!is.null(x$description)) cat("  ", x$description, "\n")
  invisible(x)
}

#' The fitter result contract
#'
#' Every analysis method must return this shape, so that competing
#' methods are interchangeable downstream and performance measures need
#' no knowledge of the fitter.
#'
#' @param estimate Numeric point estimate of the estimand.
#' @param se Numeric standard error.
#' @param p_value Numeric p-value for the null that the estimand is
#'   zero.
#' @param ci_lower,ci_upper Numeric confidence limits. Computed from
#'   `estimate`, `se` and `df` when not supplied.
#' @param df Numeric denominator degrees of freedom, or `Inf` for a
#'   normal reference.
#' @param converged Logical.
#' @param level Numeric confidence level used when limits are derived.
#' @return An object of class `fit_result`, a list with
#'   \describe{
#'     \item{estimate}{the point estimate}
#'     \item{se}{the standard error}
#'     \item{p_value}{the p-value}
#'     \item{ci_lower, ci_upper}{the confidence limits; when either was
#'       not supplied, **both** are recomputed as
#'       `estimate +/- crit * se`, with `crit` a `t` quantile on `df`
#'       degrees of freedom, or a normal quantile when `df` is
#'       infinite}
#'     \item{df}{the denominator degrees of freedom used}
#'     \item{converged}{logical convergence flag; non-converged
#'       replicates are excluded by [compute_performance()]}
#'     \item{level}{the confidence level}
#'   }
#'   These seven names are the whole contract: [run_simulation()] reads
#'   `estimate`, `se`, `p_value`, `ci_lower`, `ci_upper` and
#'   `converged` from every fitter, and nothing else.
#' @export
fit_result <- function(estimate = NA_real_, se = NA_real_,
                       p_value = NA_real_,
                       ci_lower = NULL, ci_upper = NULL,
                       df = Inf, converged = TRUE, level = 0.95) {
  if (is.null(ci_lower) || is.null(ci_upper)) {
    crit <- if (is.finite(df)) {
      stats::qt(1 - (1 - level) / 2, df = df)
    } else {
      stats::qnorm(1 - (1 - level) / 2)
    }
    ci_lower <- estimate - crit * se
    ci_upper <- estimate + crit * se
  }
  structure(list(estimate = estimate, se = se, p_value = p_value,
                 ci_lower = ci_lower, ci_upper = ci_upper,
                 df = df, converged = converged, level = level),
            class = "fit_result")
}

#' A failed fit
#'
#' @param converged Logical, normally `FALSE`.
#' @return A `fit_result` of `NA`s.
#' @export
null_fit <- function(converged = FALSE) {
  fit_result(NA_real_, NA_real_, NA_real_,
             ci_lower = NA_real_, ci_upper = NA_real_,
             converged = converged)
}

#' ANCOVA of change from baseline at a nominated visit
#'
#' A reference fitter satisfying the contract. Regresses the change
#' from baseline at `visit_time` on arm and the baseline value.
#'
#' @param dat Data with `id`, `time`, `arm`, `y`.
#' @param visit_time Numeric. Visit at which to test. Defaults to the
#'   last.
#' @param reference The reference arm level.
#' @param level Confidence level.
#' @return A `fit_result` for the contrast of the non-reference arm
#'   against the reference.
#' @export
fit_ancova <- function(dat, visit_time = NULL, reference = NULL,
                       level = 0.95) {
  if (!is.data.frame(dat)) {
    stop("`dat` must be a data frame of one row per subject-visit.")
  }
  need <- c("id", "time", "arm", "y")
  absent <- setdiff(need, names(dat))
  if (length(absent)) {
    stop("`dat` is missing required column(s): ",
         paste(absent, collapse = ", "), ".")
  }
  if (is.null(visit_time)) visit_time <- max(dat$time)
  base <- dat[abs(dat$time) < .Machine$double.eps^0.5, ]
  post <- dat[abs(dat$time - visit_time) < .Machine$double.eps^0.5, ]
  if (!nrow(base) || !nrow(post)) return(null_fit())

  m <- merge(post[, c("id", "arm", "y")],
             base[, c("id", "y")], by = "id",
             suffixes = c("", "_base"))
  m <- m[stats::complete.cases(m$y, m$y_base), ]
  if (nrow(m) < 5L) return(null_fit())

  m$arm <- factor(m$arm)
  if (!is.null(reference)) m$arm <- stats::relevel(m$arm, reference)
  if (nlevels(m$arm) < 2L) return(null_fit())
  m$chg <- m$y - m$y_base

  fit <- tryCatch(stats::lm(chg ~ arm + y_base, data = m),
                  error = function(e) NULL)
  if (is.null(fit)) return(null_fit())

  cf <- summary(fit)$coefficients
  arm_trms <- grep("^arm", rownames(cf), value = TRUE)
  trm <- arm_trms[1]
  if (is.na(trm)) return(null_fit())
  # A `fit_result` carries a single contrast. With three or more arms
  # the model fits several, and returning the first without saying so
  # would silently discard the rest.
  if (length(arm_trms) > 1L) {
    warning("`fit_ancova()` returns a single contrast, but ",
            length(arm_trms), " were fitted (",
            paste(arm_trms, collapse = ", "), "). Reporting ", trm,
            " only; call it once per contrast, with `reference` set, ",
            "to obtain the others.")
  }

  fit_result(estimate = cf[trm, "Estimate"],
             se = cf[trm, "Std. Error"],
             p_value = cf[trm, "Pr(>|t|)"],
             df = fit$df.residual, converged = TRUE, level = level)
}

#' Run a simulation study
#'
#' Drives the replicate loop, storing per-replicate estimates together
#' with the RNG state that produced them.
#'
#' @param B Integer. Number of replicates.
#' @param generate Function of no arguments returning one simulated
#'   data set.
#' @param analyse Function of a data set returning a [fit_result()], or
#'   a named list of such functions to compare methods on the *same*
#'   data (common random numbers).
#' @param estimand An [estimand()].
#' @param seed Integer. Master seed.
#' @return An object of class `c("sim_results", "data.frame")` with
#'   `B * length(analyse)` rows, one per replicate and method, and
#'   columns
#'   \describe{
#'     \item{rep}{replicate index, `1` to `B`}
#'     \item{method}{the name of the element of `analyse` that produced
#'       the row}
#'     \item{estimate, se, p_value, ci_lower, ci_upper, converged}{the
#'       corresponding elements of that method's [fit_result()]}
#'   }
#'   carrying the attributes
#'   \describe{
#'     \item{estimand}{the [estimand()] supplied, from which
#'       [compute_performance()] takes the true value}
#'     \item{streams}{the list of `B` RNG states, one per replicate, for
#'       replay with [with_rng_state()]}
#'     \item{seed}{the master seed}
#'     \item{B}{the number of replicates}
#'     \item{errors}{a character vector of trapped error messages, each
#'       tagged with its replicate number and whether it arose in
#'       `generate` or `analyse`. Empty when every replicate succeeded,
#'       and reported by the `print` method when not.}
#'   }
#' @details
#' Errors and warnings inside `generate` or `analyse` are trapped: the
#' replicate is recorded as non-converged rather than aborting the run,
#' and the message is retained in the `errors` attribute. Because the
#' RNG state for every replicate is stored, a failing replicate can be
#' reproduced exactly with [with_rng_state()].
#'
#' Passing several methods analyses the identical data set with each,
#' which is the paired comparison Morris et al. recommend.
#' @examples
#' sch <- trial_schedule(treatment = 4, interval = 3)
#' arm <- factor(rep(c("placebo", "active"), each = 30),
#'               levels = c("placebo", "active"))
#' d <- runin_design(sch, arm, reference = "placebo")
#' g <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
#' bt <- c(x_slope = 0.5, x_trt_active = -0.25)
#'
#' res <- run_simulation(
#'   B = 20,
#'   generate = function() {
#'     generate_outcomes(d, g, beta = bt, intercept = 20)
#'   },
#'   analyse = fit_ancova,
#'   estimand = estimand("difference in change at month 12", -3))
#' res
#' head(compute_performance(res), 4)
#' @export
run_simulation <- function(B, generate, analyse, estimand,
                           seed = 20260810L) {
  stopifnot(B >= 1, is.function(generate))
  if (is.function(analyse)) analyse <- list(method = analyse)
  if (!is.list(analyse) || !length(analyse)) {
    stop("`analyse` must be a function or a non-empty list of ",
         "functions, one per analysis method.")
  }
  # Every element must be named, not merely some: results are collected
  # by name, so an unnamed element is silently dropped and surfaces
  # later as a confusing row-count mismatch.
  nms <- names(analyse)
  unnamed <- if (is.null(nms)) seq_along(analyse) else which(!nzchar(nms))
  if (length(unnamed)) {
    stop("`analyse` must be a fully named list; element(s) ",
         paste(unnamed, collapse = ", "), " have no name.")
  }
  if (!all(vapply(analyse, is.function, logical(1)))) {
    stop("Every element of `analyse` must be a function; ",
         "element(s) ",
         paste(names(analyse)[!vapply(analyse, is.function,
                                      logical(1))],
               collapse = ", "), " are not.")
  }
  stopifnot(inherits(estimand, "estimand"))

  streams <- sim_streams(B, seed)
  rows <- vector("list", B)
  # An environment rather than `<<-`: the trapping code is evaluated in
  # this frame, and `<<-` searches from the *parent* environment, so it
  # would skip a local and silently write into the namespace.
  log <- new.env(parent = emptyenv())
  log$errors <- character(0)

  for (b in seq_len(B)) {
    res <- with_rng_state(streams[[b]], {
      dat <- tryCatch(generate(), error = function(e) e)
      if (inherits(dat, "error")) {
        log$errors <- c(log$errors,
                        paste0("rep ", b, " generate: ",
                               conditionMessage(dat)))
        lapply(analyse, function(f) null_fit())
      } else {
        lapply(analyse, function(f) {
          out <- tryCatch(f(dat), error = function(e) e)
          if (inherits(out, "error")) {
            log$errors <- c(log$errors,
                            paste0("rep ", b, " analyse: ",
                                   conditionMessage(out)))
            null_fit()
          } else out
        })
      }
    })
    rows[[b]] <- do.call(rbind, lapply(names(res), function(nm) {
      r <- res[[nm]]
      data.frame(rep = b, method = nm, estimate = r$estimate,
                 se = r$se, p_value = r$p_value,
                 ci_lower = r$ci_lower, ci_upper = r$ci_upper,
                 converged = r$converged,
                 stringsAsFactors = FALSE)
    }))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "estimand") <- estimand
  attr(out, "streams") <- streams
  attr(out, "seed") <- seed
  attr(out, "B") <- B
  attr(out, "errors") <- log$errors
  class(out) <- c("sim_results", "data.frame")
  out
}

#' @export
print.sim_results <- function(x, ...) {
  es <- attr(x, "estimand")
  cat("<sim_results> B =", attr(x, "B"),
      " methods:", paste(unique(x$method), collapse = ", "), "\n")
  cat("  estimand:", es$name, "=", format(es$value, digits = 6), "\n")
  cat("  converged:", sum(x$converged), "/", nrow(x), "\n")
  if (length(attr(x, "errors"))) {
    cat("  errors:", length(attr(x, "errors")), "\n")
  }
  invisible(x)
}
