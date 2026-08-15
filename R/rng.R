# Random number discipline.
#
# Morris, White and Crowther (2019) section 4.1: pin the generator
# kind, set the seed once, and store enough state that any individual
# replicate can be re-executed in isolation for diagnosis.
#
# L'Ecuyer-CMRG substreams are used rather than per-replicate
# `set.seed()`. Independent substreams are guaranteed by construction,
# and the same substreams can be replayed across design cells, which is
# what common random numbers require. Naive `set.seed(i)` gives streams
# with no such guarantee and silently desynchronizes when the amount of
# random number consumption changes between cells.

#' Generate independent RNG substreams
#'
#' @param n Integer. Number of substreams, normally the number of
#'   replicates.
#' @param seed Integer. Master seed.
#' @return A list of `n` RNG state vectors, each usable with
#'   [with_rng_state()].
#' @details
#' Pins `RNGkind("L'Ecuyer-CMRG")` for the duration of the call and
#' restores the previous kind on exit, so calling this does not change
#' the caller's generator.
#' @examples
#' streams <- sim_streams(3, seed = 20260810L)
#' length(streams)
#'
#' # Each substream is independent, and replaying one reproduces its
#' # draws exactly.
#' a <- with_rng_state(streams[[2]], rnorm(3))
#' b <- with_rng_state(streams[[2]], rnorm(3))
#' identical(a, b)
#' identical(a, with_rng_state(streams[[3]], rnorm(3)))
#' @export
sim_streams <- function(n, seed) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("`n` must be a single number of substreams of at least 1, ",
         "but it is ", paste(class(n), collapse = "/"),
         " of length ", length(n),
         if (is.numeric(n) && length(n) == 1L) {
           paste0(" with value ", format(n))
         } else "",
         ". Pass the number of replicates to be run, e.g. n = 1000.",
         call. = FALSE)
  }
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  s <- .Random.seed
  out <- vector("list", n)
  for (i in seq_len(n)) {
    out[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  out
}

#' Evaluate an expression under a stored RNG state
#'
#' Restores the caller's generator state afterwards, so a replicate can
#' be re-run in isolation without disturbing anything else.
#'
#' @param state An RNG state vector from [sim_streams()].
#' @param expr Expression to evaluate.
#' @return The value of `expr`.
#' @examples
#' # The replay workflow. `run_simulation()` stores one RNG state per
#' # replicate, so any single replicate can be regenerated on its own
#' # for diagnosis, without re-running the whole study.
#' sch <- trial_schedule(treatment = 4, interval = 3)
#' arm <- factor(rep(c("placebo", "active"), each = 20),
#'               levels = c("placebo", "active"))
#' d <- runin_design(sch, arm, reference = "placebo")
#' g <- dgm_conditional(G = diag(c(9, 0.04)), sigma2 = 4)
#' bt <- c(x_slope = 0.5, x_trt_active = -0.25)
#' gen <- function() generate_outcomes(d, g, beta = bt, intercept = 20)
#'
#' res <- run_simulation(B = 10, generate = gen, analyze = fit_ancova,
#'                       estimand = estimand("change diff", -3))
#'
#' # Re-execute replicate 3 in isolation and recover its estimate.
#' rep3 <- with_rng_state(attr(res, "streams")[[3]], gen())
#' all.equal(fit_ancova(rep3)$estimate, res$estimate[3])
#' @export
with_rng_state <- function(state, expr) {
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", state, envir = globalenv())
  force(expr)
}
