#' Combine multiple estimates of the same QAP model
#'
#' Pools several runs of the SAME model on the SAME data -- for example two
#' runs of 500 permutations each, giving a 1000-permutation result. Runs are
#' weighted by their \code{reps}, so they need not all use the same number.
#' For a model fitted with \code{comparison=}, each comparison is pooled
#' separately and all of them are kept.
#'
#' @param res list; either a single fitted model of class \code{QAPCSS},
#'   \code{QAPRegression} or \code{QAPGLM}, or a \code{list} of several such
#'   models that should all be combined.
#' @param res2 list; if \code{res} is a single model output, \code{res2} is the other model output that should be combined with the first.
#'
#' @returns Returns a \code{list} of the \code{class} of the combined model output.
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
  if (inherits(res, c("QAPCSS", "QAPRegression", "QAPGLM")) &&
      !is.null(res2)) {
    res <- list(res, res2)
  }
  n_res <- length(res)
  if (n_res < 2)
    stop("combine_qap_estimates() needs at least two model outputs.")

  return_res <- res[[1]]

  # --- compatibility checks --------------------------------------------
  ref_comp <- names(return_res$comp)
  for (i in seq_len(n_res)[-1]) {
    if (!identical(class(res[[i]]), class(return_res)))
      stop("All models must be of the same class; model ", i, " is a '",
           class(res[[i]])[1], "' but model 1 is a '",
           class(return_res)[1], "'.")
    if (!identical(res[[i]]$nullhyp, return_res$nullhyp))
      stop("All models must use the same nullhyp; model ", i, " uses '",
           res[[i]]$nullhyp, "' but model 1 uses '",
           return_res$nullhyp, "'.")
    if (!identical(names(res[[i]]$comp), ref_comp))
      stop("All models must use the same comparisons.")
    if (!identical(dimnames(res[[i]]$abs), dimnames(return_res$abs)))
      stop("All models must have the same coefficients in the same order.")
  }

  # Reps-weighted running mean of the permutation proportions.  For models
  # fitted with comparison=, $lower/$larger/$abs are *lists* of K matrices
  # (one per comparison) stored as fit$lower[[com]], so the pooling has to
  # index into that list -- the K comparisons stay separate throughout.
  pool <- function(a, b, wa, wb) (a * wa + b * wb) / (wa + wb)

  for (i in seq_len(n_res - 1)) {
    nxt  <- res[[i + 1]]
    w_old <- return_res$reps
    w_new <- nxt$reps

    if (is.null(return_res$comp)) {
      return_res$lower  <- pool(return_res$lower,  nxt$lower,  w_old, w_new)
      return_res$larger <- pool(return_res$larger, nxt$larger, w_old, w_new)
      return_res$abs    <- pool(return_res$abs,    nxt$abs,    w_old, w_new)
    } else {
      for (com in ref_comp) {
        return_res$lower[[com]]  <- pool(return_res$lower[[com]],
                                         nxt$lower[[com]],  w_old, w_new)
        return_res$larger[[com]] <- pool(return_res$larger[[com]],
                                         nxt$larger[[com]], w_old, w_new)
        return_res$abs[[com]]    <- pool(return_res$abs[[com]],
                                         nxt$abs[[com]],    w_old, w_new)
      }
    }

    return_res$reps <- w_old + w_new
  }

  return(return_res)
}
