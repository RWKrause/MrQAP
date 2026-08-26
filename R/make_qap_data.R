#' Internal auxiliary function to make QAP() data
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
#' @param multi_mode Logical; do rows and columns index distinct node sets?
#'   If so the matrix may be rectangular and there is no diagonal to drop.
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
                          xi   = NULL,
                          multi_mode = FALSE) {
  nx <- length(x)
  mode <- match.arg(mode, c("digraph", "graph"))

  if (perm && is.null(xi)) {
    y <- RMPerm(y, g)
  } else if (perm && !is.null(xi)) {
    x[[xi]] <- RMPerm(x[[xi]], g)
  }

  nr <- dim(y)[1]
  nc <- dim(y)[2]
  valid <- matrix(TRUE, nr, nc)
  # Reciprocity and self-ties are only defined when rows and columns index
  # the same actors, so neither the triangle nor the diagonal applies to
  # multi-mode data.
  if (!multi_mode) {
    # "graph" (undirected): each dyad must contribute exactly one
    # observation, so only the upper triangle is retained.
    if (mode == "graph") valid[lower.tri(valid)] <- FALSE
    if (!diag) diag(valid) <- FALSE
  }

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
    location = seq_len(nr * nc)[vv],
    yv       = as.vector(y)[vv]
  )
  pred$nv <- as.factor(net)

  # Row and column indices.  For multi-mode data these index different node
  # sets, which is what makes (1|sv) and (1|rv) distinct grouping factors.
  sv <- matrix(seq_len(nr), nr, nc)
  sv[!valid] <- NA
  pred$sv <- as.vector(sv)[vv]

  rv <- matrix(seq_len(nc), nr, nc, byrow = TRUE)
  rv[!valid] <- NA
  pred$rv <- as.vector(rv)[vv]

  for (var in seq_len(nx)) {
    pred[[names(x)[var]]] <- as.vector(x[[var]])[vv]
  }
  return(pred)
}
