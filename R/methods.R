# ============================================================
# Standard extractor methods.
#
# Defined on the shared "QAP" parent class, so they work identically for
# QAPRegression, QAPGLM and QAPCSS fits. The specific classes keep their
# own print methods.
# ============================================================

#' Build the coefficient table for one fitted QAP model
#'
#' Single source of truth for what \code{summary()} shows and what
#' \code{as.data.frame()} returns, so the two can never disagree.
#'
#' @param base A fit list (\code{fit$base}, or one element of it).
#' @param lower,larger,abs The p-value matrices for that fit.
#' @param family Character; model family.
#' @param nullhyp Character; \code{"qapy"} or \code{"qapspp"}.
#'
#' @return A data frame with one row per coefficient.
#' @keywords internal

qap_coef_table <- function(base, lower, larger, abs, family, nullhyp) {
  est <- base$coefficients
  tv  <- base$t
  nm  <- names(est)

  # Multinomial coefficients are a matrix; flatten to one row per
  # category-by-term combination.  as.vector() is column-major, so the terms
  # vary slowest and the categories fastest -- the same order the p-value
  # blocks below unpack in.
  if (is.matrix(est)) {
    nm  <- as.vector(outer(rownames(est), colnames(est),
                           function(a, b) paste(a, b, sep = ":")))
    # One flag per flattened element, marking which came from the intercept
    # COLUMN.  Matching on the term string would not work: the flattened
    # names read "B:(Intercept)".
    is_icpt <- rep(colnames(est) == "(Intercept)", each = nrow(est))
    est <- as.vector(est)
    tv  <- as.vector(tv)
  } else {
    is_icpt <- names(est) == "(Intercept)"
  }

  # The p-value matrices stack the b-comparisons above the t-comparisons:
  # one row of each normally, and (ncat - 1) rows of each for a multinomial
  # fit.  Take the half, not the row -- indexing row 2 of a multinomial
  # matrix returns the second CATEGORY's coefficient comparisons.
  grab <- function(m, half) {
    if (is.null(m)) return(rep(NA_real_, length(est)))
    h   <- nrow(m) %/% 2L
    blk <- if (half == 1L) seq_len(h) else h + seq_len(h)
    v   <- as.vector(m[blk, , drop = FALSE])
    if (length(v) != length(est)) return(rep(NA_real_, length(est)))
    v
  }

  out <- data.frame(
    term      = nm,
    estimate  = as.numeric(est),
    statistic = as.numeric(tv),
    p_lower   = grab(lower,  2L),
    p_upper   = grab(larger, 2L),
    p_value   = grab(abs,    2L),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Under semi-partialling the intercept is never permuted, so it has no
  # p-value; make that explicit rather than showing a stale number.
  if (nullhyp == "qapspp") {
    out[is_icpt, c("p_lower", "p_upper", "p_value")] <- NA_real_
  }

  if (!family %in% c("gaussian")) {
    out$exp_estimate <- exp(out$estimate)
    out <- out[, c("term", "estimate", "exp_estimate", "statistic",
                   "p_lower", "p_upper", "p_value")]
  }
  out
}


#' Is this a comparison model?
#' @param object A QAP fit.
#' @return Logical.
#' @keywords internal
qap_is_comparison <- function(object) !is.null(object$comp)


#' Extract model coefficients from a QAP fit
#'
#' @param object An object of class \code{QAPRegression}, \code{QAPGLM} or
#'   \code{QAPCSS}.
#' @param ... Ignored.
#'
#' @return A named numeric vector, or -- for a model fitted with
#'   \code{comparison=} -- a named list with one vector per comparison.
#' @export
#'
#' @examples
#' set.seed(1)
#' n  <- 10
#' x1 <- matrix(rnorm(n^2), n, n)
#' d  <- list(y = 0.5 * x1 + matrix(rnorm(n^2), n, n), x1 = x1)
#' fit <- QAP(y ~ x1, data = d, reps = 20)
#' coef(fit)
coef.QAP <- function(object, ...) {
  if (qap_is_comparison(object)) return(lapply(object$base, `[[`, "coefficients"))
  object$coefficients
}


#' Number of observations used by a QAP fit
#'
#' The number of dyads (or triads, for CSS data) that entered the model,
#' after dropping the diagonal, the redundant triangle in undirected mode,
#' and any missing values.
#'
#' @param object A QAP fit.
#' @param ... Ignored.
#'
#' @return Integer, or \code{NA} if the fitted model was not retained
#'   (\code{less_mem = TRUE}).
#' @export
nobs.QAP <- function(object, ...) {
  m <- object$simple_fit
  if (is.null(m)) {
    nd <- object$null_dist
    return(NA_integer_)
  }
  tryCatch(stats::nobs(m), error = function(e) NA_integer_)
}


#' Fitted values and residuals from a QAP fit
#'
#' Taken from the baseline model, i.e. the fit to the unpermuted data.
#'
#' @param object A QAP fit.
#' @param ... Passed to the underlying method.
#'
#' @return A numeric vector.
#' @export
fitted.QAP <- function(object, ...) {
  m <- object$simple_fit
  if (is.null(m))
    stop("No fitted model was retained; refit with less_mem = FALSE.")
  stats::fitted(m, ...)
}

#' @rdname fitted.QAP
#' @export
residuals.QAP <- function(object, ...) {
  m <- object$simple_fit
  if (is.null(m))
    stop("No fitted model was retained; refit with less_mem = FALSE.")
  stats::residuals(m, ...)
}


#' Log-likelihood of the baseline model
#'
#' Enables \code{AIC()} and \code{BIC()}.  These describe the baseline fit
#' only; they say nothing about the permutation inference.
#'
#' @param object A QAP fit.
#' @param ... Passed to the underlying method.
#'
#' @return A \code{logLik} object.
#' @export
logLik.QAP <- function(object, ...) {
  m <- object$simple_fit
  if (is.null(m))
    stop("No fitted model was retained; refit with less_mem = FALSE.")
  stats::logLik(m, ...)
}


#' Model-based covariance matrix of the baseline fit
#'
#' @param object A QAP fit.
#' @param ... Passed to the underlying method.
#'
#' @return The covariance matrix of the baseline model.
#'
#' @note This is the covariance matrix the underlying model reports, which
#'   assumes independent observations -- the assumption MRQAP exists to
#'   avoid.  It does \strong{not} correspond to the permutation p-values.
#'   For uncertainty consistent with those, use \code{\link{confint.QAP}},
#'   which uses the spread of the permutation distribution.
#' @export
vcov.QAP <- function(object, ...) {
  m <- object$simple_fit
  if (is.null(m))
    stop("No fitted model was retained; refit with less_mem = FALSE.")
  stats::vcov(m, ...)
}


#' Confidence intervals from the permutation distribution
#'
#' Forms intervals as \eqn{\hat\beta \pm z_{1-\alpha/2} \cdot s}, where
#' \eqn{s} is the standard deviation of the permutation distribution of that
#' coefficient.  Because the permutations preserve the dependence structure
#' of the network, \eqn{s} reflects the dyadic dependence that a model-based
#' standard error ignores.
#'
#' @section Why not percentiles of the permutation distribution:
#' The permutation distribution is generated under the null of no effect, so
#' it is centred on zero.  Its percentiles form an interval around zero, not
#' around the estimate, and are not a confidence interval.  Only its
#' \emph{spread} is used here.
#'
#' This is a normal approximation and inherits its usual caveats: it assumes
#' the sampling distribution of the coefficient is roughly symmetric and
#' that its spread does not depend strongly on the true value.  For a
#' skewed statistic, or a coefficient near a boundary, treat it as
#' indicative rather than exact.
#'
#' Under \code{nullhyp = "qapspp"} the intercept is never permuted and its
#' interval is \code{NA}.
#'
#' @param object A QAP fit.
#' @param parm Character or numeric; which coefficients. Defaults to all.
#' @param level Confidence level (default 0.95).
#' @param ... Ignored.
#'
#' @return A matrix with lower and upper bounds, or a named list of them for
#'   a comparison model.
#' @export
#'
#' @examples
#' set.seed(1)
#' n  <- 12
#' x1 <- matrix(rnorm(n^2), n, n)
#' x2 <- matrix(rnorm(n^2), n, n)
#' d  <- list(y = 0.5 * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
#' fit <- QAP(y ~ x1 + x2, data = d, reps = 200)
#' confint(fit)
confint.QAP <- function(object, parm, level = 0.95, ...) {
  # A multinomial fit's draws are never stacked (the statistic is a matrix,
  # not a b/t pair), so null_dist is absent by construction. Say that,
  # rather than sending the user to less_mem, which cannot help here.
  if (identical(object$family, "multinom"))
    stop("Permutation confidence intervals are not available for ",
         "multinomial fits: the permutation draws are a matrix per ",
         "replication and are not retained. Use summary() for the ",
         "permutation p-values.")

  one <- function(est, nd) {
    if (is.null(nd) || is.null(nd$b))
      stop("No permutation draws were retained, so a permutation-based ",
           "interval cannot be formed. Refit with less_mem = FALSE.")

    s <- apply(nd$b, 2, stats::sd, na.rm = TRUE)
    s <- s[names(est)]
    z <- stats::qnorm(1 - (1 - level) / 2)

    ci <- cbind(est - z * s, est + z * s)
    colnames(ci) <- paste(format(100 * c((1 - level) / 2,
                                         1 - (1 - level) / 2),
                                 trim = TRUE, digits = 3), "%")
    rownames(ci) <- names(est)
    ci
  }

  if (qap_is_comparison(object)) {
    out <- lapply(names(object$comp), function(cn)
      one(object$base[[cn]]$coefficients, object$null_dist[[cn]]))
    names(out) <- names(object$comp)
    if (!missing(parm))
      out <- lapply(out, function(m) m[parm, , drop = FALSE])
    return(out)
  }

  ci <- one(object$coefficients, object$null_dist)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}


#' Summarise a QAP fit
#'
#' Returns the coefficient table and the settings the model was fitted
#' under.  Printing it gives the same report as printing the fit, with
#' significance codes added.
#'
#' @param object A QAP fit.
#' @param ... Ignored.
#'
#' @return An object of class \code{summary.QAP}: a list whose
#'   \code{coefficients} element is a data frame with one row per term
#'   (\code{estimate}, \code{statistic}, and the three p-values), plus the
#'   fitting settings.  For a comparison model, \code{coefficients} is a
#'   named list of such tables.
#' @export
#'
#' @examples
#' set.seed(1)
#' n  <- 12
#' x1 <- matrix(rnorm(n^2), n, n)
#' x2 <- matrix(rnorm(n^2), n, n)
#' d  <- list(y = 0.5 * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
#' fit <- QAP(y ~ x1 + x2, data = d, reps = 200)
#'
#' summary(fit)
#' summary(fit)$coefficients      # programmatic access
summary.QAP <- function(object, ...) {
  tab <- if (qap_is_comparison(object)) {
    stats::setNames(lapply(names(object$comp), function(cn)
      qap_coef_table(object$base[[cn]], object$lower[[cn]],
                     object$larger[[cn]], object$abs[[cn]],
                     object$family, object$nullhyp)),
      names(object$comp))
  } else {
    qap_coef_table(object, object$lower, object$larger, object$abs,
                   object$family, object$nullhyp)
  }

  out <- list(
    coefficients  = tab,
    family        = object$family,
    nullhyp       = object$nullhyp,
    reps          = object$reps,
    mode          = object$mode,
    diag          = object$diag,
    css           = isTRUE(object$css),
    multi_mode    = isTRUE(object$multi_mode),
    grouped       = !is.null(object$groups),
    robust_se     = isTRUE(object$robust_se),
    estimator     = object$estimator,
    random        = object$random,
    comp          = object$comp,
    reference     = object$reference,
    r.squared     = object$r.squared,
    adj.r.squared = object$adj.r.squared,
    theta         = object$theta,
    zi_coefficients   = object$zi_coefficients,
    random.intercepts = object$random.intercepts,
    confusion_matrix  = object$confusion_matrix,
    nobs          = suppressWarnings(nobs.QAP(object))
  )
  class(out) <- "summary.QAP"
  out
}


#' @param x A \code{summary.QAP} object.
#' @param digits Number of significant digits.
#' @param signif.stars Show significance codes?
#' @rdname summary.QAP
#' @export
print.summary.QAP <- function(x, digits = max(3L, getOption("digits") - 3L),
                              signif.stars = getOption("show.signif.stars",
                                                       TRUE),
                              ...) {
  kind <- if (x$css) "Cognitive Social Structure" else if (x$multi_mode)
    "Multi-Mode Network" else "Network"
  fam  <- if (x$family == "gaussian") "Linear" else "Generalized Linear"
  mixed <- if (!is.null(x$random.intercepts)) " Mixed" else ""

  cat("\n", fam, mixed, " ", kind, " Model (MRQAP)\n\n", sep = "")

  cat("Family:      ", x$family, "\n", sep = "")
  cat("Estimator:   ",
      if (identical(x$estimator, "gmm")) "generalized method-of-moments"
      else "standard", "\n", sep = "")
  cat("Outcome:     ", x$mode,
      if (x$diag) "; diagonal included" else "; diagonal ignored",
      "\n", sep = "")
  cat("Null:        ",
      if (identical(x$nullhyp, "qapspp"))
        "Dekker semi-partialling plus" else "outcome permutation",
      " (", format(x$reps), " permutations",
      if (identical(x$nullhyp, "qapspp")) " per predictor" else "",
      ")\n", sep = "")
  if (x$grouped)
    cat("Permutation: restricted to within groups\n")
  if (x$robust_se)
    cat("Std. errors: heteroskedasticity-consistent\n")
  if (!is.na(x$nobs))
    cat("Observations:", format(x$nobs), "\n")

  if (!is.null(x$r.squared))
    cat("\nR-squared: ", format(round(x$r.squared, 4)),
        "   Adjusted: ", format(round(x$adj.r.squared, 4)), "\n", sep = "")
  if (!is.null(x$theta))
    cat("\nNegative binomial dispersion (theta): ",
        format(round(x$theta, 4)), "\n", sep = "")
  if (!is.null(x$zi_coefficients))
    cat("\nZero-inflation coefficients: ",
        paste(names(x$zi_coefficients),
              format(round(x$zi_coefficients, 4)),
              sep = " = ", collapse = ", "), "\n", sep = "")

  if (is.list(x$coefficients) && !is.data.frame(x$coefficients)) {
    for (cn in names(x$coefficients)) {
      cat("\n--- Comparison: ", cn, " (", x$comp[[cn]][1], " vs ",
          x$comp[[cn]][2], ") ---\n", sep = "")
      qap_print_coefs(x$coefficients[[cn]], digits, signif.stars, x$reps)
    }
  } else {
    cat("\nCoefficients:\n")
    qap_print_coefs(x$coefficients, digits, signif.stars, x$reps)
  }

  if (identical(x$nullhyp, "qapspp"))
    cat("\nThe intercept is not tested under semi-partialling plus.\n")
  cat("p-values are proportions of ", format(x$reps),
      " permutations; the smallest distinguishable value is ",
      formatC(1 / x$reps, format = "f",
              digits = max(2, ceiling(log10(x$reps)))), ".\n", sep = "")

  if (!is.null(x$confusion_matrix)) {
    cat("\n")
    print(x$confusion_matrix)
  }
  cat("\n")
  invisible(x)
}


#' Format a permutation p-value honestly
#'
#' A permutation p-value is a proportion out of \code{reps} draws, so its
#' resolution is 1/reps.  An observed 0 means "no permutation was as
#' extreme", i.e. below that resolution -- not zero, and certainly not the
#' \code{<2e-16} that the usual coefficient formatter would print.
#'
#' @param p Numeric vector of p-values.
#' @param reps Integer; permutations per test.
#'
#' @return Character vector.
#' @keywords internal

format_perm_p <- function(p, reps) {
  res <- 1 / reps
  out <- formatC(p, format = "f", digits = max(2, ceiling(log10(reps))))
  out[!is.na(p) & p == 0] <- paste0("<", formatC(res, format = "f",
                                                 digits = max(2, ceiling(log10(reps)))))
  out[is.na(p)] <- "-"
  out
}


#' Render a QAP coefficient table
#'
#' Deliberately does not use \code{stats::printCoefmat()}: that formats the
#' p-value column as if it came from a continuous reference distribution and
#' would print an observed proportion of 0 as \code{<2e-16}.
#'
#' @param tab Output of \code{qap_coef_table()}.
#' @param digits Significant digits.
#' @param signif.stars Show significance codes?
#' @param reps Integer; permutations per test, used to bound the p-values.
#'
#' @return Invisibly NULL; prints.
#' @keywords internal

qap_print_coefs <- function(tab, digits = 4, signif.stars = TRUE,
                            reps = NULL) {
  # format() (not formatC) so a column shares one number of decimals and
  # the values line up on the decimal point.
  num <- function(v) format(round(v, digits), nsmall = 0, trim = FALSE)

  cols <- list(Estimate = num(tab$estimate))
  if (!is.null(tab$exp_estimate))
    cols[["exp(Est)"]] <- num(tab$exp_estimate)
  cols[["t value"]]   <- num(tab$statistic)
  if (is.null(reps)) reps <- 1000
  cols[["Pr(<=t)"]]   <- format_perm_p(tab$p_lower, reps)
  cols[["Pr(>=t)"]]   <- format_perm_p(tab$p_upper, reps)
  cols[["Pr(>=|t|)"]] <- format_perm_p(tab$p_value, reps)

  if (signif.stars) {
    p <- tab$p_value
    st <- ifelse(is.na(p), "",
          ifelse(p < 0.001, "***",
          ifelse(p < 0.01,  "**",
          ifelse(p < 0.05,  "*",
          ifelse(p < 0.1,   ".", " ")))))
    cols[[""]] <- st
  }

  m <- do.call(cbind, cols)
  rownames(m) <- tab$term
  print(noquote(m))

  if (signif.stars && any(!is.na(tab$p_value)))
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(NULL)
}


#' Coerce a QAP fit to a data frame
#'
#' The tidy form of the coefficient table: one row per term, with the
#' estimate, test statistic and permutation p-values.  Convenient for
#' assembling results across models.
#'
#' @param x A QAP fit.
#' @param row.names,optional Ignored; present for method consistency.
#' @param ... Ignored.
#'
#' @return A data frame.  For a comparison model the tables are stacked with
#'   a \code{comparison} column.
#' @export
#'
#' @examples
#' set.seed(1)
#' n  <- 10
#' x1 <- matrix(rnorm(n^2), n, n)
#' x2 <- matrix(rnorm(n^2), n, n)
#' d  <- list(y = 0.5 * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
#' as.data.frame(QAP(y ~ x1 + x2, data = d, reps = 50))
as.data.frame.QAP <- function(x, row.names = NULL, optional = FALSE, ...) {
  s <- summary(x)
  if (is.data.frame(s$coefficients)) return(s$coefficients)

  out <- do.call(rbind, lapply(names(s$coefficients), function(cn) {
    d <- s$coefficients[[cn]]
    cbind(comparison = cn, d, stringsAsFactors = FALSE)
  }))
  row.names(out) <- NULL
  out
}
