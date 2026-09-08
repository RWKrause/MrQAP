# ============================================================
# broom methods, registered only if broom is present.
#
# Registering at load time via registerS3method() keeps broom out of
# Imports: users who have it get tidy()/glance(); users who do not are
# unaffected and the package still installs without it.
# ============================================================

#' Tidy a QAP fit
#'
#' The coefficient table in \pkg{broom}'s column convention.  The
#' \code{p.value} column is the two-sided permutation p-value, the same
#' number \code{summary()} reports as \code{p_value}.
#'
#' @param x A QAP fit.
#' @param conf.int Logical; add \code{conf.low}/\code{conf.high} from
#'   \code{\link{confint.QAP}}, which uses the spread of the permutation
#'   distribution.
#' @param conf.level Confidence level for those intervals.
#' @param ... Ignored.
#'
#' @return A data frame with one row per term.
#' @keywords internal

tidy_QAP <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  tab <- as.data.frame(x)

  # broom's contract: term, estimate, std.error, statistic, p.value.
  # A QAP fit has no single standard error -- the t-value is what the
  # permutation test compares -- so std.error is deliberately absent
  # rather than filled with a model-based number the p-value ignores.
  names(tab)[names(tab) == "p_value"] <- "p.value"

  if (conf.int) {
    ci <- tryCatch(confint(x, level = conf.level), error = function(e) NULL)
    if (!is.null(ci)) {
      tab$conf.low  <- ci[match(tab$term, rownames(ci)), 1]
      tab$conf.high <- ci[match(tab$term, rownames(ci)), 2]
    }
  }
  tab
}


#' One-row summary of a QAP fit
#'
#' @param x A QAP fit.
#' @param ... Ignored.
#'
#' @return A one-row data frame.
#' @keywords internal

glance_QAP <- function(x, ...) {
  num <- function(v) if (is.null(v)) NA_real_ else as.numeric(v)
  quiet <- function(e) NA_real_

  data.frame(
    family        = if (is.null(x$family)) NA_character_ else x$family,
    nullhyp       = if (is.null(x$nullhyp)) NA_character_ else x$nullhyp,
    reps          = num(x$reps),
    nobs          = suppressWarnings(nobs(x)),
    r.squared     = num(x$r.squared),
    adj.r.squared = num(x$adj.r.squared),
    AIC           = tryCatch(as.numeric(AIC(x)), error = quiet),
    BIC           = tryCatch(as.numeric(BIC(x)), error = quiet),
    stringsAsFactors = FALSE
  )
}


# Register with broom if -- and only if -- it is installed.
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("broom", quietly = TRUE)) {
    registerS3method("tidy",   "QAP", tidy_QAP,   envir = asNamespace("broom"))
    registerS3method("glance", "QAP", glance_QAP, envir = asNamespace("broom"))
  }
  invisible()
}
