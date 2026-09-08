#' Combine multiple estimates of the same QAP model
#'
#' Pools several runs of the SAME model on the SAME data -- for example two
#' runs of 500 permutations each, giving a 1000-permutation result. Runs are
#' weighted by how many permutations actually contributed, so they need not
#' all use the same number. For a model fitted with \code{comparison=}, each
#' comparison is pooled separately and all of them are kept.
#'
#' @section What must match:
#' Pooling is only meaningful for repeated runs of one model on one dataset,
#' so the runs are checked before they are combined: the class, the null
#' hypothesis, the family, the estimator, the outcome mode, the diagonal and
#' data-shape flags, the comparisons, the coefficient names, and the
#' coefficients themselves. The baseline fit does not depend on the
#' permutations, so identical data implies identical coefficients -- and
#' coefficients that differ prove the runs are not poolable, whatever their
#' settings say.
#'
#' @param res list; either a single fitted model of class \code{QAP}, or a
#'   \code{list} of several such models that should all be combined.
#' @param res2 list; if \code{res} is a single model output, \code{res2} is
#'   the other model output that should be combined with the first.
#'
#' @returns Returns a \code{list} of the \code{class} of the combined model
#'   output. \code{reps} is the total across runs, and \code{null_dist}
#'   holds every run's permutation draws stacked, so \code{confint()} on the
#'   result uses all of them.
#' @export
#'
#' @examples
#' set.seed(1)
#' n  <- 10
#' x1 <- matrix(rnorm(n^2), n, n)
#' d  <- list(y = 0.5 * x1 + matrix(rnorm(n^2), n, n), x1 = x1)
#'
#' # Two independent runs of 25 permutations...
#' a <- QAP(y ~ x1, data = d, reps = 25, seed = 1)
#' b <- QAP(y ~ x1, data = d, reps = 25, seed = 2)
#'
#' # ...pooled into a single 50-permutation result.
#' ab <- combine_qap_estimates(a, b)
#' ab$reps

combine_qap_estimates <- function(res, res2 = NULL) {
  # A fit is itself a list, so a single one passed alone would otherwise be
  # read as a list of models -- one per element of the fit -- and get well
  # past the length check before failing incomprehensibly.
  if (inherits(res, "QAP")) {
    if (is.null(res2))
      stop("combine_qap_estimates() needs at least two model outputs: ",
           "pass the second as res2, or pass a list of fits as res.")
    res <- list(res, res2)
  }

  n_res <- length(res)
  if (n_res < 2)
    stop("combine_qap_estimates() needs at least two model outputs.")

  not_fits <- !vapply(res, inherits, logical(1), what = "QAP")
  if (any(not_fits))
    stop("Every element of res must be a QAP fit; element(s) ",
         paste(which(not_fits), collapse = ", "), " are not.")

  return_res <- res[[1]]
  ref_comp   <- names(return_res$comp)

  # --- compatibility checks --------------------------------------------
  # Settings that must be identical for the runs to describe one model.
  # A mismatch here means the p-values being averaged answer different
  # questions, which no amount of weighting can reconcile.
  settings <- c("nullhyp", "family", "estimator", "mode", "diag",
                "css", "multi_mode", "robust_se", "reference")

  for (i in seq_len(n_res)[-1]) {
    ri <- res[[i]]

    if (!identical(class(ri), class(return_res)))
      stop("All models must be of the same class; model ", i, " is a '",
           class(ri)[1], "' but model 1 is a '", class(return_res)[1], "'.")

    for (s in settings) {
      if (!identical(ri[[s]], return_res[[s]])) {
        shown <- function(v) if (is.null(v)) "NULL" else format(v)
        stop("All models must be fitted with the same ", s, "; model ", i,
             " used '", shown(ri[[s]]), "' but model 1 used '",
             shown(return_res[[s]]), "'.")
      }
    }

    if (!identical(names(ri$comp), ref_comp))
      stop("All models must use the same comparisons.")
    if (!identical(dimnames(ri$abs), dimnames(return_res$abs)))
      stop("All models must have the same coefficients in the same order.")

    # The decisive check. The baseline fit is computed on the unpermuted
    # data, so two runs of one model on one dataset must agree on it
    # exactly; if they do not, the runs are of different data and pooling
    # their p-values would silently invent a result.
    if (!isTRUE(all.equal(ri$coefficients, return_res$coefficients)))
      stop("Model ", i, " has different coefficients from model 1, so the ",
           "two were not fitted to the same data. Only repeated runs of ",
           "the same model on the same data can be pooled.")
  }

  # --- pooling ----------------------------------------------------------
  # Weight by how many permutations CONTRIBUTED, not by how many were
  # requested: a run that lost permutations to non-convergence has a
  # smaller denominator, and weighting by reps would over-weight it.
  weight_of <- function(f) {
    if (!is.null(f$n_valid)) return(f$n_valid)
    warning("Model fitted before n_valid was recorded; weighting by reps. ",
            "If any permutations failed, the pooled p-values are slightly ",
            "off. Refit to remove this warning.")
    f$reps
  }

  # n_valid is a scalar under qapy and one entry per predictor under
  # qapspp, where each column of the result was computed from its own set
  # of permutations. Weight column by column so a predictor that lost
  # permutations in one run does not drag the others.
  #
  # A column with no weight -- the intercept under qapspp, or a predictor
  # that was skipped -- pools to NA, which is what its p-values already are.
  pool <- function(a, b, wa, wb) {
    k  <- ncol(a)
    cn <- colnames(a)
    # Column count comes from ncol(), not from the names: an unnamed matrix
    # would otherwise spread to length zero and take the whole result with it.
    spread <- function(w) if (length(w) > 1L && !is.null(cn)) {
      unname(w[cn])
    } else {
      rep(unname(w[1L]), k)
    }
    w_a <- spread(wa)
    w_b <- spread(wb)

    num <- a * rep(w_a, each = nrow(a)) + b * rep(w_b, each = nrow(b))
    num / rep(w_a + w_b, each = nrow(a))
  }

  # Stack the permutation draws so confint() on the pooled object uses
  # every run's, not just the first's.
  bind_draws <- function(a, b) {
    if (is.null(a) || is.null(b)) return(NULL)
    list(b = rbind(a$b, b$b), t = rbind(a$t, b$t))
  }

  for (i in seq_len(n_res - 1)) {
    nxt   <- res[[i + 1]]
    w_old <- weight_of(return_res)
    w_new <- weight_of(nxt)

    if (is.null(return_res$comp)) {
      return_res$lower  <- pool(return_res$lower,  nxt$lower,  w_old, w_new)
      return_res$larger <- pool(return_res$larger, nxt$larger, w_old, w_new)
      return_res$abs    <- pool(return_res$abs,    nxt$abs,    w_old, w_new)
      return_res$null_dist <- bind_draws(return_res$null_dist, nxt$null_dist)
    } else {
      for (com in ref_comp) {
        return_res$lower[[com]]  <- pool(return_res$lower[[com]],
                                         nxt$lower[[com]],  w_old, w_new)
        return_res$larger[[com]] <- pool(return_res$larger[[com]],
                                         nxt$larger[[com]], w_old, w_new)
        return_res$abs[[com]]    <- pool(return_res$abs[[com]],
                                         nxt$abs[[com]],    w_old, w_new)
        return_res$null_dist[[com]] <- bind_draws(return_res$null_dist[[com]],
                                                  nxt$null_dist[[com]])
      }
    }

    return_res$reps    <- return_res$reps + nxt$reps
    return_res$n_valid <- w_old + w_new
  }

  return_res
}
