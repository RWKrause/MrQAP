# ============================================================
# Fast path for the plain linear model.
#
# The general engine rebuilds the whole vectorised data frame and refits
# via the formula interface on every permutation, although only one column
# ever changes and the model is an OLS. This module precomputes everything
# invariant once and solves the normal equations directly.
#
# It must produce numerically identical results to the general path -- see
# tests/testthat/test-golden.R -- and consume randomness in the same order,
# which it does because it calls perm_networks() exactly once per rep.
# ============================================================

#' Is the plain-OLS fast path usable?
#'
#' Requires a linear model with no random effects, no absorbed fixed
#' effects, no comparisons, and complete data.  Completeness matters because
#' the fast path treats the set of valid cells as invariant: that holds when
#' the only exclusions are structural (the diagonal, or one triangle), since
#' a permutation maps the diagonal to the diagonal.  Any \code{NA} makes the
#' valid set depend on the permutation.
#'
#' @param family,estimator,has_random,use_fixest,comparison Model settings.
#' @param data Named list of matrices or arrays.
#' @param vars Character; the variables that will be vectorised.
#'
#' @return Logical.
#' @keywords internal

qap_ols_eligible <- function(family, estimator, has_random, use_fixest,
                             comparison, data, vars) {
  # Escape hatch used by the tests to check that the fast and general paths
  # agree, and available to users if the fast path is ever suspected.
  if (isTRUE(getOption("MrQAP.disable_fast_ols", FALSE))) return(FALSE)

  if (family != "gaussian" || estimator != "standard") return(FALSE)
  if (has_random || use_fixest || !is.null(comparison)) return(FALSE)

  for (v in vars) {
    obj <- data[[v]]
    if (is.list(obj)) {
      if (any(vapply(obj, anyNA, logical(1)))) return(FALSE)
    } else if (anyNA(obj)) {
      return(FALSE)
    }
  }
  TRUE
}


#' Precompute everything that does not change between permutations
#'
#' Builds the valid-cell index and an extractor that turns a (permuted)
#' matrix or array into the model's response/predictor column directly,
#' with no data frame construction.
#'
#' @param data Named list of matrices or arrays.
#' @param dep Character; dependent variable name.
#' @param main Character; predictor names, in model-matrix order.
#' @param diag,mode,css,large,multi_mode Structure flags.
#'
#' @return A list with \code{extract}, \code{X}, \code{y}, \code{n_obs} and
#'   \code{p}.
#' @keywords internal

qap_ols_template <- function(data, dep, main, diag, mode, css, large,
                             multi_mode) {

  # Vectorise one network exactly the way make_qap_data()/make_css_data() do.
  av <- if (css) {
    function(a) array_to_vector(a, mode. = mode, diag. = diag)
  } else {
    function(m) as.vector(m)
  }

  valid_of <- function(y) {
    d <- dim(y)
    if (css) {
      v <- array(TRUE, dim = d)
      if (!multi_mode) {
        if (!diag) for (i in seq_len(d[3])) diag(v[, , i]) <- FALSE
        if (mode == "undirected")
          for (i in seq_len(d[3])) v[, , i][lower.tri(v[, , i])] <- FALSE
      }
      av(v)
    } else {
      v <- matrix(TRUE, d[1], d[2])
      if (!multi_mode) {
        if (mode == "undirected") v[lower.tri(v)] <- FALSE
        if (!diag) diag(v) <- FALSE
      }
      as.vector(v)
    }
  }

  if (!large) {
    vv <- list(valid_of(data[[dep]]))
    extract <- function(obj) av(obj)[vv[[1]]]
  } else {
    vv <- lapply(data[[dep]], valid_of)
    extract <- function(obj) {
      unlist(lapply(seq_along(obj), function(i) av(obj[[i]])[vv[[i]]]),
             use.names = FALSE)
    }
  }

  y <- extract(data[[dep]])
  X <- cbind("(Intercept)" = 1,
             vapply(main, function(v) extract(data[[v]]), numeric(length(y))))
  colnames(X) <- c("(Intercept)", main)

  list(extract = extract, X = X, y = y, n_obs = length(y), p = ncol(X))
}


#' Solve an OLS fit and return coefficients and t-values
#'
#' @param X Model matrix (with intercept column).
#' @param y Response vector.
#' @param qrX Optional precomputed \code{qr(X)}; supplied when X is fixed
#'   across permutations.
#' @param xtxd Optional precomputed \code{diag(chol2inv(qr.R(qrX)))}.
#' @param robust Logical; use HC3 standard errors.  Must match what the
#'   baseline fit used, or the permutation t-values would be compared
#'   against a baseline computed on a different scale.
#'
#' @return A list with \code{coefficients} and \code{t}, named as the
#'   columns of X.
#' @keywords internal

qap_ols_solve <- function(X, y, qrX = NULL, xtxd = NULL, robust = FALSE) {
  if (is.null(qrX)) {
    qrX  <- qr(X)
    xtxd <- diag(chol2inv(qr.R(qrX)))
  }
  b   <- qr.coef(qrX, y)
  res <- y - as.vector(X %*% b)

  t <- if (robust) {
    # Same quantity fit_qap_model() produces for an lm via robust_se().
    b / HC3(X[, -1, drop = FALSE], res)
  } else {
    s2 <- sum(res^2) / (nrow(X) - ncol(X))
    b / sqrt(xtxd * s2)
  }

  names(b) <- names(t) <- colnames(X)
  list(coefficients = b, t = t)
}


#' Precompute the Frisch-Waugh-Lovell decomposition for one tested predictor
#'
#' Under \code{qapspp} only one column of the design changes per
#' permutation, and only that column's coefficient and t-value are ever
#' read -- \code{compare_perm_to_baseline(..., xi =)} discards the rest. A
#' full \code{qr()} of the whole design per replication therefore computes
#' \code{p - 1} coefficients to throw them away.
#'
#' Frisch-Waugh-Lovell gives the same number from a decomposition of the
#' \emph{fixed} columns, formed once here: with \eqn{Z} the untested columns
#' and \eqn{\tilde{y}}, \eqn{\tilde{x}} the residuals of \eqn{y} and the
#' tested column on \eqn{Z},
#' \deqn{b = \tilde{x}'\tilde{y} / \tilde{x}'\tilde{x}}
#' and the full model's residual sum of squares is
#' \eqn{RSS_y - b^2 \tilde{x}'\tilde{x}}. Both are exact, not approximations.
#'
#' @param X Model matrix including the intercept column.
#' @param y Response vector.
#' @param col Integer; the column of \code{X} that the permutation replaces.
#' @param robust Logical; precompute the leverages HC3 needs.
#'
#' @return A list consumed by \code{qap_fwl_solve()}.
#' @keywords internal

qap_fwl_template <- function(X, y, col, robust = FALSE) {
  Z   <- X[, -col, drop = FALSE]
  qrZ <- qr(Z)

  out <- list(qrZ = qrZ,
              y_t = qr.resid(qrZ, y),
              df  = nrow(X) - ncol(X),
              robust = robust)
  out$rss_y <- sum(out$y_t^2)

  if (robust) {
    # Leverages of the fixed columns. The full design's leverages are
    # h = h_Z + x_t^2 / sxx, so only this part is invariant.
    Q     <- qr.Q(qrZ)
    out$h_Z <- rowSums(Q^2)
  }
  out
}


#' Solve one permuted column via Frisch-Waugh-Lovell
#'
#' @param fw Output of \code{qap_fwl_template()}.
#' @param xcol Numeric vector; the permuted predictor column.
#'
#' @return A list with \code{b} and \code{t} for that column alone.
#' @keywords internal

qap_fwl_solve <- function(fw, xcol) {
  x_t <- qr.resid(fw$qrZ, xcol)
  sxx <- sum(x_t^2)
  b   <- sum(x_t * fw$y_t) / sxx

  if (!fw$robust) {
    # Residual SS of the FULL model, without forming its residuals.
    s2 <- (fw$rss_y - b^2 * sxx) / fw$df
    return(list(b = b, t = b / sqrt(s2 / sxx)))
  }

  # HC3 for a single coefficient is the sandwich specialised through the
  # same residualised regressor: sum(x_t^2 * omega) / sxx^2.
  e  <- fw$y_t - b * x_t          # full-model residuals
  h  <- fw$h_Z + x_t^2 / sxx      # full-model leverages
  v  <- sum(x_t^2 * e^2 / (1 - h)^2) / sxx^2
  list(b = b, t = b / sqrt(v))
}


#' Compare one predictor's permuted statistics against the baseline
#'
#' The single-column counterpart of \code{compare_perm_to_baseline()}, for
#' the FWL path, which never forms the other coefficients.
#'
#' @param b,t Numeric; the permuted coefficient and t-value.
#' @param base_b,base_t Numeric; the same from the baseline fit.
#'
#' @return A list shaped exactly like
#'   \code{compare_perm_to_baseline(..., xi =)}.
#' @keywords internal

qap_compare_one <- function(b, t, base_b, base_t) {
  pres <- stats::setNames(c(b, t), qap_stat_rows())
  bres <- stats::setNames(c(base_b, base_t), qap_stat_rows())
  list(lower  = pres <= bres,
       larger = pres >= bres,
       abs    = abs(pres) >= abs(bres),
       draw   = pres)
}


#' Run the permutation loop on the fast OLS path
#'
#' @return A list with \code{lower}, \code{larger}, \code{abs}.
#' @keywords internal
#' @noRd

qap_ols_perms <- function(tmpl, data, dep, main, groups, reps, base_fit,
                          nullhyp, css, multi_mode, pred, valid, valid_list,
                          large, rand_part, has_random, mode, diag, ncores,
                          use_robust_errors = FALSE) {

  if (!is.null(ncores)) {
    restore <- setup_future_plan(ncores)
    on.exit(restore(), add = TRUE)
  }

  extract <- tmpl$extract
  X       <- tmpl$X

  body <- function(p) {
    if (nullhyp == "qapy") {
      # X is fixed for the entire run, so its decomposition is computed once
      # and every permutation is a back-solve.
      qrX  <- qr(X)
      xtxd <- diag(chol2inv(qr.R(qrX)))
      yobj <- data[[dep]]

      res <- run_permutations(
        reps,
        function(i) {
          yp <- perm_networks(yobj, groups, CSS = css,
                              multi_mode = multi_mode)
          f  <- qap_ols_solve(X, extract(yp), qrX, xtxd,
                              robust = use_robust_errors)
          compare_perm_to_baseline(f$coefficients, f$t, base_fit, xi = NULL)
        },
        p = p
      )
      return(aggregate_perm_results(res, reps))
    }

    # --- qapspp: one column of X is replaced each rep ---
    out     <- qap_init_pmats(base_fit, NULL)
    per_xi  <- list()
    n_valid <- c()

    for (xi in main) {
      d <- qap_residualise(xi, data, pred, main, has_random, rand_part,
                           valid, valid_list, large, css, mode)
      if (is.null(d)) next

      xobj <- d[[xi]]
      col  <- match(xi, colnames(X))

      # Only this column's statistics are ever read, so decompose the fixed
      # columns once and solve for the one coefficient per replication.
      fw     <- qap_fwl_template(X, tmpl$y, col, robust = use_robust_errors)
      base_b <- base_fit$coefficients[[xi]]
      base_t <- base_fit$t[[xi]]

      res <- run_permutations(
        reps,
        function(i) {
          xp <- perm_networks(xobj, groups, CSS = css,
                              multi_mode = multi_mode)
          f  <- qap_fwl_solve(fw, extract(xp))
          qap_compare_one(f$b, f$t, base_b, base_t)
        },
        p = p
      )
      agg <- aggregate_perm_results(res, reps)
      n_valid[xi] <- agg$n_valid
      out$lower[,  xi] <- agg$lower
      out$larger[, xi] <- agg$larger
      out$abs[,    xi] <- agg$abs
      per_xi[[xi]] <- agg$draws
    }
    out$draws   <- qap_draws_by_predictor(per_xi, colnames(out$abs), reps)
    out$n_valid <- n_valid
    out
  }

  qap_with_progressor(
    body, steps = if (nullhyp == "qapy") reps else reps * length(main))
}
