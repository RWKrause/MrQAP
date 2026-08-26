#' Perform HC3 correction for heteroskedasticity
#'
#' @param X matrix; Matrix of predictors
#' @param e numeric; vector of residuals
#'
#' @returns Numeric vector of HC3 standard errors, one per column of
#'   \code{cbind(1, X)}.
#'
#' @details Valid for linear models only. The leverages are computed from the
#'   unweighted \code{X}, so this must not be applied to a GLM fit; use
#'   \code{sandwich::vcovHC()} there instead.
#' @keywords internal

HC3 <- function(X, e) {
  XO     <- cbind(matrix(1, dim(X)[1], 1), X)
  XTXINV <- solve(t(XO) %*% XO)

  # Leverages = diag(XO %*% XTXINV %*% t(XO)), without forming the n x n hat
  # matrix or looping over rows.
  h  <- rowSums((XO %*% XTXINV) * XO)
  om <- e^2 / (1 - h)^2
  sqrt(diag(t(t(XTXINV %*% t(XO)) * om) %*% XO %*% XTXINV))
}
