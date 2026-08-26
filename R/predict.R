# ============================================================
# Prediction, and a picture of the permutation distribution.
#
# Both read the baseline fit -- the one to the unpermuted data. The
# permutation machinery tests coefficients; it plays no part in prediction.
# ============================================================

#' Predictions from a QAP fit
#'
#' @section Getting a network back:
#' \code{type = "matrix"} returns the predictions in the shape of the
#' network rather than as a vector, with \code{NA} in every cell that did
#' not enter the model: the diagonal, the redundant triangle in undirected
#' mode, and any cell dropped for missingness. This is usually what you
#' want in order to plot a predicted network or compare it with the
#' observed one.
#'
#' @param object A QAP fit.
#' @param newdata Optional named list of matrices or arrays, shaped like the
#'   data the model was fitted to, holding new values for the predictors.
#'   \code{NULL} (default) predicts for the fitted data.
#' @param type Character; \code{"response"} (default) or \code{"link"}, as
#'   for \code{\link[stats]{predict.glm}}, or \code{"matrix"} to reshape the
#'   response-scale predictions back into the network's dimensions.
#' @param ... Passed to the underlying predict method.
#'
#' @return A numeric vector, or -- for \code{type = "matrix"} -- a matrix or
#'   array, or a list of them when the model was fitted to several networks.
#' @export
#'
#' @examples
#' set.seed(1)
#' n  <- 10
#' x1 <- matrix(rnorm(n^2), n, n)
#' d  <- list(y = 0.5 * x1 + matrix(rnorm(n^2), n, n), x1 = x1)
#' fit <- QAP(y ~ x1, data = d, reps = 20)
#'
#' head(predict(fit))
#' dim(predict(fit, type = "matrix"))
predict.QAP <- function(object, newdata = NULL,
                        type = c("response", "link", "matrix"), ...) {
  type <- match.arg(type)

  if (qap_is_comparison(object))
    stop("predict() is not defined for a model fitted with comparison=: ",
         "each comparison is a separate model. Use ",
         "predict(fit$simple_fits[[i]]) on the one you want.")

  m <- object$simple_fit
  if (is.null(m))
    stop("No fitted model was retained; refit with less_mem = FALSE.")

  pv <- if (is.null(newdata)) {
    if (type == "link") stats::predict(m, type = "link", ...)
    else                stats::fitted(m, ...)
  } else {
    stats::predict(m, newdata = qap_predict_frame(object, newdata),
                   type = if (type == "matrix") "response" else type, ...)
  }
  pv <- as.vector(pv)

  if (type != "matrix") return(pv)

  # Put the vector back where it came from. The mask is rebuilt from the
  # shape recorded at fit time, so this needs neither the original data nor
  # the model frame.
  dims <- if (is.null(newdata)) object$dim else qap_newdata_dim(newdata)
  qap_vector_to_shape(pv, dims, object)
}


#' Vectorise newdata the way the fit vectorised its own data
#'
#' @param object A QAP fit.
#' @param newdata Named list of matrices or arrays.
#'
#' @return A data frame with the predictor columns the model expects.
#' @keywords internal

qap_predict_frame <- function(object, newdata) {
  if (!is.list(newdata) || is.null(names(newdata)))
    stop("newdata must be a named list of matrices or arrays.")

  wanted <- all.vars(stats::formula(object$simple_fit))
  vars   <- intersect(names(newdata), wanted)
  if (!length(vars))
    stop("newdata contains none of the model's predictors (",
         paste(setdiff(wanted, wanted[1]), collapse = ", "), ").")

  probe <- newdata[[vars[1]]]
  large <- is.list(probe)

  one <- function(v, net) {
    obj  <- if (large) newdata[[v]][[net]] else newdata[[v]]
    keep <- qap_valid_mask(dim(obj), css = isTRUE(object$css),
                           mode = object$mode, diag = isTRUE(object$diag),
                           multi_mode = isTRUE(object$multi_mode))
    av <- if (isTRUE(object$css)) {
      array_to_vector(obj, mode. = object$mode, diag. = isTRUE(object$diag))
    } else {
      as.vector(obj)
    }
    av[keep]
  }

  nets <- if (large) seq_along(probe) else 1L
  out  <- lapply(vars, function(v)
    unlist(lapply(nets, function(i) one(v, i)), use.names = FALSE))
  names(out) <- vars
  as.data.frame(out, stringsAsFactors = FALSE)
}


#' Dimensions implied by a newdata list
#'
#' @param newdata Named list of matrices or arrays.
#' @return A dim vector, or a list of them.
#' @keywords internal

qap_newdata_dim <- function(newdata) {
  probe <- newdata[[1]]
  if (is.list(probe)) lapply(probe, dim) else dim(probe)
}


#' Reshape a vector of predictions back into network shape
#'
#' @param v Numeric vector, one element per modelled cell.
#' @param dims A dim vector, or a list of them for several networks.
#' @param object A QAP fit, for the structural flags.
#'
#' @return A matrix/array, or a list of them.
#' @keywords internal

qap_vector_to_shape <- function(v, dims, object) {
  mask_of <- function(d)
    qap_valid_mask(d, css = isTRUE(object$css), mode = object$mode,
                   diag = isTRUE(object$diag),
                   multi_mode = isTRUE(object$multi_mode))

  fill <- function(d, vals) {
    keep <- mask_of(d)
    out  <- rep(NA_real_, length(keep))
    # Cells excluded for missingness leave fewer values than the structural
    # mask has slots; fill what there is and leave the rest NA.
    n <- min(sum(keep), length(vals))
    out[which(keep)[seq_len(n)]] <- vals[seq_len(n)]
    array(out, dim = d)
  }

  if (!is.list(dims)) return(fill(dims, v))

  # Several networks: the vector is their modelled cells concatenated.
  ns     <- vapply(dims, function(d) sum(mask_of(d)), integer(1))
  ends   <- cumsum(ns)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  lapply(seq_along(dims), function(i) fill(dims[[i]], v[starts[i]:ends[i]]))
}


#' Plot the permutation distribution of a QAP fit
#'
#' One panel per tested coefficient: a histogram of the permutation draws
#' with the observed statistic marked. This is the picture the p-value
#' summarises, and the quickest way to see whether an effect sits in the
#' tail of its null distribution or in the middle of it.
#'
#' @param x A QAP fit.
#' @param which Character or numeric; which coefficients to draw. Defaults
#'   to every tested one. The intercept is excluded under
#'   \code{nullhyp = "qapspp"}, which never permutes it.
#' @param statistic Character; plot the distribution of the \code{"t"}
#'   values (default) or of the coefficients, \code{"b"}.
#' @param breaks Passed to \code{\link[graphics]{hist}}.
#' @param ... Passed to \code{\link[graphics]{hist}}.
#'
#' @return Invisibly \code{x}.
#' @export
#'
#' @examples
#' set.seed(1)
#' n  <- 12
#' x1 <- matrix(rnorm(n^2), n, n)
#' x2 <- matrix(rnorm(n^2), n, n)
#' d  <- list(y = 0.6 * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
#' fit <- QAP(y ~ x1 + x2, data = d, reps = 100)
#'
#' op <- graphics::par(mfrow = c(1, 2))
#' plot(fit)
#' graphics::par(op)
plot.QAP <- function(x, which = NULL, statistic = c("t", "b"),
                     breaks = 30, ...) {
  statistic <- match.arg(statistic)

  if (qap_is_comparison(x))
    stop("plot() is not defined for a model fitted with comparison=. ",
         "Plot one comparison at a time from fit$null_dist[[name]].")
  if (identical(x$family, "multinom"))
    stop("plot() is not available for multinomial fits: the permutation ",
         "draws are a matrix per replication and are not retained.")

  nd <- x$null_dist
  if (is.null(nd) || is.null(nd[[statistic]]))
    stop("No permutation draws were retained, so there is nothing to plot. ",
         "Refit with less_mem = FALSE.")

  draws <- nd[[statistic]]
  obs   <- if (statistic == "t") x$t else x$coefficients

  # A column is plottable when it was actually permuted; the intercept
  # under qapspp is all NA, exactly as its p-values are.
  ok <- colnames(draws)[colSums(!is.na(draws)) > 0]
  if (is.null(which))     which <- ok
  if (is.numeric(which))  which <- colnames(draws)[which]

  bad <- setdiff(which, ok)
  if (length(bad))
    stop("No permutation draws for: ", paste(bad, collapse = ", "),
         if (identical(x$nullhyp, "qapspp") && "(Intercept)" %in% bad)
           ". The intercept is not tested under semi-partialling plus."
         else "")

  for (v in which) {
    d <- draws[, v]
    d <- d[!is.na(d)]
    o <- unname(obs[[v]])
    graphics::hist(
      d, breaks = breaks,
      xlim = range(c(d, o), finite = TRUE),
      main = v,
      xlab = paste0("permuted ", statistic, "   (p = ",
                    format_perm_p(x$abs[2L, v], x$reps), ")"),
      ...)
    graphics::abline(v = o, lwd = 2, col = "red")
    graphics::box()
  }
  invisible(x)
}
