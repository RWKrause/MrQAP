#' Permute a network matrix or CSS array
#'
#' Applies a random relabelling of nodes, the permutation step underlying QAP
#' inference.  Rows and columns (and for CSS, perceiver slices) are permuted
#' together, so the structure of the network is preserved and only the
#' correspondence between nodes and positions is randomised.
#'
#' For one-mode data every dimension indexes the same actors, so a single node
#' order is drawn and applied to all of them.  For multi-mode data each
#' dimension is its own node set and receives an independent order.
#'
#' @param m A \code{matrix}, \code{array}, or \code{list} of either.
#' @param groups Restricts permutation to occur only within groups.  Either a
#'   vector applying to every dimension, or -- for multi-mode data -- a named
#'   list with elements \code{sender}, \code{receiver} and (CSS)
#'   \code{perceiver} giving one grouping per dimension.  \code{NULL} permutes
#'   freely.
#' @param CSS Logical; is \code{m} a cognitive social structure, i.e. should
#'   one node order apply across all three dimensions jointly?  When
#'   \code{FALSE} a 3D array is treated as a stack of independent networks and
#'   each slice is permuted separately.
#' @param multi_mode Logical; does each dimension index a distinct node set
#'   (2- or 3-mode data)?  Defaults to \code{FALSE}, i.e. one-mode data, which
#'   requires a square matrix or a cubic array.
#'
#' @returns The permuted object, matching the structure of \code{m}.
#' @export
#'
#' @examples
#' m <- matrix(1:16, 4, 4)
#' RMPerm(m)
#'
#' # permute only within groups
#' RMPerm(m, groups = c(1, 1, 2, 2))
#'
#' # two-mode data: rows and columns are permuted independently
#' RMPerm(matrix(1:12, 3, 4), multi_mode = TRUE)

RMPerm <- function(m, groups = NULL, CSS = FALSE, multi_mode = FALSE) {

  if (is.list(m)) {
    return(lapply(m, RMPerm, groups = groups, CSS = CSS,
                  multi_mode = multi_mode))
  }

  d  <- dim(m)
  nd <- length(d)
  if (!nd %in% c(2L, 3L))
    stop("RMPerm() needs a matrix or a 3-dimensional array; got ",
         nd, " dimension(s).")

  if (!multi_mode && length(unique(d)) != 1L)
    stop("RMPerm() requires a square matrix or a cubic array; got ",
         paste(d, collapse = " x "),
         ". Use multi_mode = TRUE if each dimension is its own node set.")

  gs <- rmperm_dim_groups(groups, d, multi_mode)

  # sample() must never be handed a length-1 vector: sample(7) returns a
  # permutation of 1:7 rather than 7 itself, so a permutation group with a
  # single member would otherwise yield out-of-range indices.
  shuffle <- function(ix) ix[sample.int(length(ix))]

  draw <- function(k) {
    n <- d[k]
    g <- gs[[k]]
    # The unrestricted case is by far the most common and needs no
    # split()/unsplit().  It must still draw exactly one sample.int(n) so the
    # RNG stream matches the grouped path.
    if (is.null(g) || length(unique(g)) == 1L) return(sample.int(n))
    unsplit(lapply(split(seq_len(n), g), shuffle), g)
  }

  if (nd == 2L) {
    if (multi_mode) {
      o1 <- draw(1L); o2 <- draw(2L)
    } else {
      o1 <- o2 <- draw(1L)
    }
    return(matrix(m[o1, o2], nrow = d[1], ncol = d[2]))
  }

  p <- array(dim = d)

  if (CSS) {
    if (multi_mode) {
      o1 <- draw(1L); o2 <- draw(2L); o3 <- draw(3L)
    } else {
      o1 <- o2 <- o3 <- draw(1L)
    }
    p[, , ] <- array(m[o1, o2, o3])
    return(p)
  }

  # Not CSS: a stack of independent networks, one per first-dimension slice.
  for (i in seq_len(d[1])) {
    if (multi_mode) {
      o2 <- draw(2L); o3 <- draw(3L)
    } else {
      o2 <- o3 <- draw(2L)
    }
    p[i, , ] <- array(m[i, o2, o3])
  }
  p
}


#' Resolve a grouping specification to one vector per dimension
#'
#' @param groups Vector, named list (\code{sender}/\code{receiver}/
#'   \code{perceiver}), or NULL.
#' @param d Integer vector of dimensions.
#' @param multi_mode Logical.
#'
#' @return A list of length \code{length(d)}; each element is a grouping
#'   vector, or NULL for unrestricted permutation of that dimension.
#' @keywords internal

rmperm_dim_groups <- function(groups, d, multi_mode = FALSE) {
  nd <- length(d)
  if (is.null(groups)) return(vector("list", nd))

  dim_names <- c("sender", "receiver", "perceiver")[seq_len(nd)]

  if (is.list(groups)) {
    if (is.null(names(groups)) || any(!nzchar(names(groups))))
      stop("A list of groupings must be named: use ",
           "list(sender = ..., receiver = ...) to give one grouping per ",
           "dimension.")
    unknown <- setdiff(names(groups), dim_names)
    if (length(unknown))
      stop("Unknown grouping name(s): ", paste(unknown, collapse = ", "),
           ". Expected any of: ", paste(dim_names, collapse = ", "), ".")

    out <- vector("list", nd)
    for (k in seq_len(nd)) {
      g <- groups[[dim_names[k]]]
      if (!is.null(g)) {
        if (length(g) != d[k])
          stop("groups$", dim_names[k], " has length ", length(g),
               " but dimension ", k, " has ", d[k], " nodes.")
        out[[k]] <- as.character(g)
      }
    }
    return(out)
  }

  # A single vector applies to every dimension.
  g <- as.character(groups)
  if (any(d != length(g)))
    stop("groups has length ", length(g),
         " but the data have dimensions ", paste(d, collapse = " x "),
         ". Supply one grouping per dimension, e.g. ",
         "list(sender = ..., receiver = ...).")
  rep(list(g), nd)
}
