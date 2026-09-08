# ============================================================
# Input adapters.
#
# Everything downstream of QAP() works on matrices and arrays. This is the
# one place that turns the other shapes people actually have -- a dyadic
# data frame, igraph objects, network objects -- into that.
# ============================================================

#' Coerce QAP input into a named list of matrices or arrays
#'
#' Accepts, in addition to the native named list of matrices:
#' \itemize{
#'   \item a \strong{data frame} of dyadic (or triadic) observations, which
#'     is handed to \code{\link{df_to_mat}} -- the function that already
#'     does exactly this, including \code{multi_mode} and \code{split_by};
#'   \item list elements that are \pkg{igraph} or \pkg{network} objects,
#'     which are coerced to adjacency matrices.
#' }
#'
#' @param data A named list, or a data frame.
#' @param sender,receiver,perceiver Column names identifying the
#'   dimensions, required when \code{data} is a data frame.
#' @param split_by Optional column name to split a data frame by, giving
#'   one network per level.
#' @param mode Character; passed to \code{df_to_mat()}.
#' @param diag Logical; passed to \code{df_to_mat()} as \code{loops}.
#' @param multi_mode Logical or NULL; passed to \code{df_to_mat()}.
#' @param attr Character or NULL; for \pkg{igraph}/\pkg{network} input, the
#'   edge attribute to read.  \code{NULL} gives a binary adjacency matrix.
#'
#' @return A named list of matrices or arrays.
#' @keywords internal

qap_coerce_data <- function(data, sender = NULL, receiver = NULL,
                            perceiver = NULL, split_by = NULL,
                            mode = "directed", diag = FALSE,
                            multi_mode = NULL, attr = NULL) {

  # --- a dyadic data frame -----------------------------------------------
  if (is.data.frame(data)) {
    if (is.null(sender) || is.null(receiver))
      stop("When `data` is a data frame, `sender` and `receiver` must name ",
           "the columns identifying the dimensions, e.g. ",
           "QAP(y ~ x, data = df, sender = \"from\", receiver = \"to\").")

    return(df_to_mat(data, sender = sender, receiver = receiver,
                     perceiver = perceiver, mode = mode, loops = diag,
                     multi_mode = isTRUE(multi_mode), split_by = split_by))
  }

  if (!is.list(data))
    stop("`data` must be a named list of matrices or arrays, or a data ",
         "frame of dyadic observations.")

  # --- graph objects inside a list ---------------------------------------
  conv <- function(x) {
    if (inherits(x, "igraph"))  return(qap_from_igraph(x, attr))
    if (inherits(x, "network")) return(qap_from_network(x, attr))
    x
  }

  lapply(data, function(el) {
    # Convert BEFORE testing for a list of networks: igraph and network
    # objects are themselves lists, so descending into them first iterates
    # their internal structure instead of converting them.
    if (inherits(el, c("igraph", "network"))) return(conv(el))

    # A plain list of networks: convert each element. A matrix or array is
    # not a list, so it falls through untouched.
    if (is.list(el) && !is.data.frame(el) && is.null(dim(el)))
      return(lapply(el, conv))

    el
  })
}


#' Adjacency matrix from an igraph object
#'
#' @param g An \pkg{igraph} graph.
#' @param attr Character or NULL; edge attribute to use as the weight.
#'
#' @return A numeric matrix with vertex names as dimnames.
#' @keywords internal

qap_from_igraph <- function(g, attr = NULL) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required to use igraph objects as input ",
         "(install.packages('igraph')).")
  if (!is.null(attr) && !attr %in% igraph::edge_attr_names(g))
    stop("igraph object has no edge attribute '", attr, "'.")

  m <- as.matrix(igraph::as_adjacency_matrix(g, attr = attr, sparse = FALSE))
  storage.mode(m) <- "double"
  m
}


#' Adjacency matrix from a network object
#'
#' @param g A \pkg{network} object.
#' @param attr Character or NULL; edge attribute to use as the weight.
#'
#' @return A numeric matrix.
#' @keywords internal

qap_from_network <- function(g, attr = NULL) {
  if (!requireNamespace("network", quietly = TRUE))
    stop("Package 'network' is required to use network objects as input ",
         "(install.packages('network')).")

  m <- if (is.null(attr)) {
    network::as.matrix.network.adjacency(g)
  } else {
    network::as.matrix.network.adjacency(g, attrname = attr)
  }
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  m
}
