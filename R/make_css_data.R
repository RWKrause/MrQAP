#' Internal auxiliary function to shape data into internal QAP() format
#'
#' @param y array or list; see \code{QAP()}
#' @param x array or list; see \code{QAP()}
#' @param nets numeric; see \code{QAP()}
#' @param diag logical; see \code{QAP()}
#' @param mode character; see \code{QAP()}
#' @param multi_mode logical; does each dimension index a distinct node set?
#'   If so the array may be non-cubic, and neither the diagonal nor the
#'   triangle applies.
#'
#' @returns Returns vectorized CSS data
#' @keywords internal

make_css_data <- function(y,
                          x,
                          nets,
                          diag,
                          mode,
                          multi_mode = FALSE) {
  d  <- dim(y)
  nx <- length(x)
  valid <- array(TRUE, dim = d)

  # Self-ties and reciprocity are only defined when the sender and receiver
  # dimensions index the same actors.
  if (!multi_mode) {
    if (!diag) {
      for (i in seq_len(d[3])) {
        diag(y[, , i]) <- NA
        for (var in seq_len(nx)) diag(x[[var]][, , i]) <- NA
      }
    }
  }

  valid[is.na(y)] <- FALSE
  for (var in seq_len(nx)) valid[is.na(x[[var]])] <- FALSE

  if (!multi_mode && mode == 'undirected') {
    for (i in seq_len(d[3])) {
      y[, , i][lower.tri(y[, , i])] <- NA
      valid[, , i][lower.tri(valid[, , i])] <- FALSE
      for (var in seq_len(nx)) {
        x[[var]][, , i][lower.tri(x[[var]][, , i])] <- NA
      }
    }
  }

  y[!valid] <- NA
  for (var in seq_len(nx)) x[[var]][!valid] <- NA

  av <- function(a) array_to_vector(a, mode. = mode, diag. = diag)

  vv <- av(valid)
  pred <- data.frame(yv = av(y)[vv], nv = nets)

  # slice.index() gives the index along a dimension without a loop and
  # without assuming the array is cubic.
  pred$sv <- as.factor(av(slice.index(valid, 1L))[vv])
  pred$rv <- as.factor(av(slice.index(valid, 2L))[vv])
  pred$pv <- as.factor(av(slice.index(valid, 3L))[vv])

  for (var in seq_len(nx)) pred[[names(x)[var]]] <- av(x[[var]])[vv]

  list(pred = pred, valid = valid)
}
