# ============================================================
# Permutation drivers shared by network and CSS data.
#
# The former QAPglm() and QAPcss() entry points each carried a
# near-identical copy of everything below.
# ============================================================

#' Initialise the p-value matrices for qapspp
#'
#' Under qapspp each predictor is tested in its own set of permutations, so
#' the result matrices are filled one column at a time and start as NA.  The
#' intercept column is never filled: it has no meaning under semi-partialling.
#'
#' The shape is read off the baseline coefficients, so a multinomial fit --
#' whose coefficients are a matrix -- needs no special-casing by the caller.
#'
#' @param base_fit Baseline fit, or list of fits when comparisons are used.
#' @param comparison Comparison list or NULL.
#'
#' @return A list with \code{lower}, \code{larger}, \code{abs}.
#' @keywords internal

qap_init_pmats <- function(base_fit, comparison) {
  mk <- function(f) {
    cf <- f$coefficients
    # A multinomial fit's coefficients are a (ncat - 1) x k MATRIX, so the
    # column count is ncol(cf), not length(cf), and the names live in
    # colnames(cf), not names(cf).  Using the vector accessors built a
    # matrix (ncat - 1) times too wide with no column names, and the
    # per-predictor assignment out$lower[, xi] then failed outright.
    if (is.matrix(cf)) {
      cats <- rownames(cf)
      m <- matrix(NA_real_, nrow = 2L * length(cats), ncol = ncol(cf),
                  dimnames = list(c(paste0("b:", cats), paste0("t:", cats)),
                                  colnames(cf)))
      return(m)
    }
    m <- matrix(NA_real_, nrow = 2L, ncol = length(cf),
                dimnames = list(qap_stat_rows(), names(cf)))
    m
  }

  if (is.null(comparison)) {
    m <- mk(base_fit)
    return(list(lower = m, larger = m, abs = m))
  }

  per <- lapply(base_fit, mk)
  names(per) <- names(comparison)
  list(lower = per, larger = per, abs = per)
}


#' Aggregate permutation results into proportions
#'
#' Handles the plain case and the comparison case, where each permutation
#' returns one result per comparison.
#'
#' @param res List of permutation results (some possibly NULL).
#' @param reps Number of permutations attempted.
#' @param comparison Comparison list or NULL.
#'
#' @return A list with \code{lower}, \code{larger}, \code{abs}.
#' @keywords internal

qap_aggregate <- function(res, reps, comparison) {
  if (is.null(comparison)) return(aggregate_perm_results(res, reps))

  res_valid <- Filter(Negate(is.null), res)
  n_valid   <- length(res_valid)
  if (n_valid == 0L) stop("All permutations failed to converge.")
  if (n_valid < reps)
    warning(reps - n_valid, " of ", reps,
            " permutations failed and were excluded.")

  resL <- unlist(unlist(res_valid, recursive = FALSE), recursive = FALSE)
  pick <- function(cn, what)
    Reduce("+", resL[names(resL) == paste0(cn, ".", what)], 0) / n_valid

  cns <- names(comparison)
  out <- list(lower  = lapply(cns, pick, what = "lower"),
              larger = lapply(cns, pick, what = "larger"),
              abs    = lapply(cns, pick, what = "abs"))
  out$lower  <- stats::setNames(out$lower,  cns)
  out$larger <- stats::setNames(out$larger, cns)
  out$abs    <- stats::setNames(out$abs,    cns)

  out$draws <- stats::setNames(
    lapply(cns, function(cn)
      qap_stack_draws(resL[names(resL) == paste0(cn, ".draw")])),
    cns)
  out$n_valid <- n_valid
  out
}


#' Assemble the per-predictor null distributions of a qapspp run
#'
#' Under semi-partialling each predictor gets its own set of permutations,
#' so the draws arrive one column at a time and are slotted into a matrix
#' shaped like the coefficient vector.  Columns never tested -- the
#' intercept -- stay NA, exactly as their p-values do.
#'
#' @param per_xi Named list of \code{qap_stack_draws()} results, one per
#'   predictor.
#' @param coef_names Character vector of all coefficient names.
#' @param reps Integer; permutations requested per predictor.
#'
#' @return A list with \code{b} and \code{t}, each \code{reps x k}.
#' @keywords internal

qap_draws_by_predictor <- function(per_xi, coef_names, reps) {
  per_xi <- Filter(Negate(is.null), per_xi)
  if (!length(per_xi)) return(NULL)

  mk <- function(which) {
    m <- matrix(NA_real_, nrow = reps, ncol = length(coef_names),
                dimnames = list(NULL, coef_names))
    for (xi in names(per_xi)) {
      v <- per_xi[[xi]][[which]][, 1]
      m[seq_along(v), xi] <- v
    }
    m
  }
  list(b = mk("b"), t = mk("t"))
}


#' Residualise a predictor and write the residuals back into the data
#'
#' The qapspp step: regress \code{xi} on the other predictors, then put the
#' residuals back into the network structure so they can be permuted.
#'
#' @param xi Character; the predictor being tested.
#' @param data Named list of matrices or arrays.
#' @param pred Vectorised data frame.
#' @param main Character vector of predictor names.
#' @param has_random,rand_part Random-effect handling for the residualisation.
#' @param valid,valid_list Validity masks (CSS).
#' @param large,css,mode Structure flags.
#'
#' @return \code{data} with \code{xi} replaced by its residuals, or NULL when
#'   \code{xi} is not numeric and cannot be residualised.
#' @keywords internal

qap_residualise <- function(xi, data, pred, main, has_random, rand_part,
                            valid, valid_list, large, css, mode) {
  probe <- if (large) data[[xi]][[1]] else data[[xi]]
  if (!is.numeric(probe)) {
    warning("Cannot residualise non-numeric predictor '", xi,
            "'. Skipping qapspp for this variable.")
    return(NULL)
  }

  xR <- residualise_predictor(xi, pred, main, has_random = has_random,
                              rand_formula = rand_part)

  out <- data
  out[[xi]] <- if (css) {
    residuals_to_array(xR, data[[xi]], valid, pred, large, valid_list,
                       mode = mode)
  } else {
    residuals_to_matrix(xR, data[[xi]], pred, large,
                        mode = if (mode == "directed") "digraph" else "graph")
  }
  out
}


#' Run the CPU permutation loop
#'
#' @return A list with \code{lower}, \code{larger}, \code{abs}.
#' @keywords internal
#' @noRd

qap_cpu_perms <- function(data, parsed, mode, diag, groups, reps, base_fit,
                          nullhyp, main, data_vars, pred, valid, valid_list,
                          large, rand_part, mod, family, estimator,
                          use_fixest, fixest_se_cluster, use_robust_errors,
                          has_random, comparison, reference, ncores, css,
                          multi_mode = FALSE) {

  # Only touch the future plan when the caller asked for a core count;
  # otherwise leave whatever plan they set in place.
  if (!is.null(ncores)) {
    restore <- setup_future_plan(ncores)
    on.exit(restore(), add = TRUE)
  }

  dep  <- parsed$dependent
  mode_internal <- if (css) mode else {
    if (mode == "directed") "digraph" else "graph"
  }

  perm_args <- function(d, perm_var) {
    list(data. = d, perm_var. = perm_var, mode. = mode_internal,
         diag. = diag, mod. = mod, groups. = groups, fit. = base_fit,
         family. = family, estimator. = estimator, use_fixest. = use_fixest,
         fixest_se_cluster. = fixest_se_cluster,
         use_robust_errors. = use_robust_errors, has_random. = has_random,
         main_vars. = main, data_vars. = data_vars, parsed. = parsed,
         comp. = comparison, reference. = reference, css. = css,
         multi_mode. = multi_mode)
  }

  body <- function(p) {
    if (nullhyp == "qapy") {
      res <- do.call(run_permutations,
                     c(list(reps, QAPPermEst), perm_args(data, NULL),
                       list(p = p)))
      return(qap_aggregate(res, reps, comparison))
    }

    # --- qapspp: one set of permutations per predictor ---
    out <- qap_init_pmats(base_fit, comparison)
    per_xi  <- list()
    n_valid <- c()

    for (xi in main) {
      d <- qap_residualise(xi, data, pred, main, has_random, rand_part,
                           valid, valid_list, large, css, mode)
      if (is.null(d)) next

      res <- do.call(run_permutations,
                     c(list(reps, QAPPermEst), perm_args(d, xi),
                       list(p = p)))
      agg <- qap_aggregate(res, reps, comparison)
      # Each predictor is tested in its own set of permutations and can
      # lose a different number of them, so the divisor is per-predictor.
      n_valid[xi] <- agg$n_valid

      if (is.null(comparison)) {
        out$lower[,  xi] <- agg$lower
        out$larger[, xi] <- agg$larger
        out$abs[,    xi] <- agg$abs
        per_xi[[xi]] <- agg$draws
      } else {
        for (cn in names(comparison)) {
          out$lower[[cn]][,  xi] <- agg$lower[[cn]]
          out$larger[[cn]][, xi] <- agg$larger[[cn]]
          out$abs[[cn]][,    xi] <- agg$abs[[cn]]
        }
        per_xi[[xi]] <- agg$draws
      }
    }

    out$draws <- if (is.null(comparison)) {
      qap_draws_by_predictor(per_xi, colnames(out$abs), reps)
    } else {
      stats::setNames(lapply(names(comparison), function(cn)
        qap_draws_by_predictor(
          lapply(per_xi, `[[`, cn), colnames(out$abs[[cn]]), reps)),
        names(comparison))
    }
    out$n_valid <- n_valid
    out
  }

  if (requireNamespace("progressr", quietly = TRUE)) {
    total_reps <- if (nullhyp == "qapy") reps else reps * length(main)
    progressr::with_progress({
      p <- progressr::progressor(steps = total_reps)
      body(p)
    })
  } else {
    body(NULL)
  }
}


#' Run the GPU permutation loop
#'
#' @return A list with \code{lower}, \code{larger}, \code{abs}.
#' @keywords internal
#' @noRd

qap_gpu_perms <- function(tmpl, data, dep, main, groups, reps, base_fit,
                          nullhyp, css, multi_mode, pred, valid, valid_list,
                          large, rand_part, has_random, mode,
                          use_robust_errors = FALSE) {

  batch <- function(d, perm_var) {
    gpu_batch_ols(tmpl = tmpl, data = d, dep = dep, main = main,
                  groups = groups, reps = reps, baseline_fit = base_fit,
                  perm_var = perm_var, css = css, multi_mode = multi_mode,
                  use_robust_errors = use_robust_errors)
  }

  if (nullhyp == "qapy") return(batch(data, NULL))

  out     <- qap_init_pmats(base_fit, NULL)
  per_xi  <- list()
  n_valid <- c()
  for (xi in main) {
    d <- qap_residualise(xi, data, pred, main, has_random, rand_part,
                         valid, valid_list, large, css, mode)
    if (is.null(d)) next
    g <- batch(d, xi)
    n_valid[xi] <- g$n_valid
    out$lower[,  xi] <- g$lower[,  xi]
    out$larger[, xi] <- g$larger[, xi]
    out$abs[,    xi] <- g$abs[,    xi]
    # Keep only the column this run actually tested.
    per_xi[[xi]] <- list(b = g$draws$b[, xi, drop = FALSE],
                         t = g$draws$t[, xi, drop = FALSE])
  }
  out$draws   <- qap_draws_by_predictor(per_xi, colnames(out$abs), reps)
  out$n_valid <- n_valid
  out
}
