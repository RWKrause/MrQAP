#' Internal auxiliary function to make QAPglm() data
#'
#' Converts network matrices into a long-format data frame suitable for
#' regression.  Also handles permutation when \code{perm = TRUE}.
#'
#' @param y Matrix; the dependent variable.
#' @param x Named list of matrices; the predictors.
#' @param g Vector; group memberships for constrained permutation.
#' @param diag Logical; include diagonal?
#' @param mode Character; "digraph" or "graph".
#' @param net Integer; network index (for multiple networks).
#' @param perm Logical; should a permutation be performed?
#' @param xi Character; name of variable to permute (qapspp) or NULL (qapy).
#'
#' @return A data frame with columns: location, yv, nv, sv, rv, and
#'   one column per predictor.
#' @keywords internal

make_qap_data <- function(y,
                          x,
                          g    = NULL,
                          diag = FALSE,
                          mode = "digraph",
                          net  = 1,
                          perm = FALSE,
                          xi   = NULL) {
  nx <- length(x)

  if (perm && is.null(xi)) {
    y <- RMPerm(y, g)
  } else if (perm && !is.null(xi)) {
    x[[xi]] <- RMPerm(x[[xi]], g)
  }

  n <- dim(y)[1]
  valid <- matrix(TRUE, n, n)
  if (!diag) diag(valid) <- FALSE

  for (var in seq_len(nx)) {
    valid[is.na(x[[var]])] <- FALSE
  }
  valid[is.na(y)] <- FALSE
  y[!valid] <- NA

  vv <- as.vector(valid)

  for (var in seq_len(nx)) {
    x[[var]][!valid] <- NA
  }

  pred <- data.frame(
    location = as.vector(matrix(seq_len(n^2), n, n))[vv],
    yv       = as.vector(y)[vv]
  )
  pred$nv <- as.factor(net)

  sv <- matrix(seq_len(n), n, n)
  sv[!valid] <- NA
  pred$sv <- as.vector(sv)[vv]

  rv <- t(matrix(seq_len(n), n, n))
  rv[!valid] <- NA
  pred$rv <- as.vector(rv)[vv]

  for (var in seq_len(nx)) {
    pred[[names(x)[var]]] <- as.vector(x[[var]])[vv]
  }
  return(pred)
}
