#' Internal auxiliary function to transform an array into a vector
#'
#' Concatenates the perceiver slices of a 3D array into a single vector, in
#' the order the rest of the package expects.
#'
#' @param ar array; The \code{array} to be transformed.
#' @param mode. character; Is the \code{array} \code{mode. = 'directed'} or \code{mode. = 'undirected'}?
#' @param diag. logical; Should the diagonal of each perceiver slice of the array be used?
#'
#' @returns Returns a vector
#' @keywords internal

array_to_vector <- function(ar, mode., diag.) {
  # Directed mode concatenates every slice in order, which is exactly the
  # array's own column-major layout.  (The previous loop ran over dimension 1
  # while slicing dimension 3, so it silently truncated any non-cubic array.)
  if (mode. != 'undirected') return(as.vector(ar))

  # Undirected: one entry per dyad, taken from the upper triangle of each
  # perceiver slice.  The index is the same for every slice, so build it once
  # and preallocate rather than growing with c().
  n_p <- dim(ar)[3]
  idx <- upper.tri(matrix(0, dim(ar)[1], dim(ar)[2]), diag = diag.)
  parts <- vector("list", n_p)
  for (i in seq_len(n_p)) parts[[i]] <- ar[, , i][idx]
  unlist(parts, use.names = FALSE)
}
