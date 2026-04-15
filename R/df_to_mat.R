#' Convert a dyadic or triaduc dataframe to a list of matrices or arrays
#'
#' Transforms a dataframe of dyadic or triadic observations into a named list of matrices
#' (or arrays if a perceiver is specified), one per variable. Optionally splits
#' the output by a grouping variable, returning a list of lists.
#'
#' @param df A dataframe containing dy/triadic data.
#' @param sender Character. Column name identifying the sender/row dimension.
#' @param receiver Character. Column name identifying the receiver/column dimension.
#' @param perceiver Character or NULL. Column name identifying the perceiver,
#'   creating a 3rd array dimension. Default \code{NULL} returns matrices.
#' @param mode Character. Either \code{"directed"} (default) for asymmetric
#'   matrices, or \code{"undirected"} to mirror values across the diagonal.
#' @param loops Logical. Whether to include self-ties (diagonal). Default \code{FALSE},
#'   which sets the diagonal to \code{NA}.
#' @param multi_mode Logical. If \code{TRUE}, treats each dimension as a distinct
#'   node set (2- or 3-mode data), allowing rectangular/non-cubic structures.
#'   If \code{FALSE} (default), all dimensions share the same nodeset (1-mode data).
#' @param split_by Character or NULL. Column name to split the data by (e.g., group or year),
#'   returning a list of lists where the outer list is by variable and the inner
#'   list is by split value. Each split may have differently sized matrices.
#'   Default \code{NULL} returns a single list of matrices.
#'
#' @returns A named list of matrices (one per variable), or a named list of lists
#'   (by variable, then by split level) if \code{split_by} is specified. If
#'   \code{perceiver} is provided, matrices are replaced by 3-dimensional arrays.
#' @export
#'
#' @importFrom purrr transpose
#'
#'
#' @examples
#' df <- data.frame(
#'   sender   = c("A", "A", "B"),
#'   receiver = c("B", "C", "C"),
#'   weight   = c(1, 2, 3)
#' )
#' df_to_mat(df, sender = "sender", receiver = "receiver")
#' df_to_mat(df, sender = "sender", receiver = "receiver", mode = "undirected")
#'
#'
df_to_mat <- function(df,
                      sender,
                      receiver,
                      perceiver  = NULL,
                      mode       = c("directed", "undirected"),
                      loops      = FALSE,
                      multi_mode = FALSE,
                      split_by   = NULL) {
  mode <- match.arg(mode)
  var_names <- setdiff(colnames(df), c(sender, receiver, perceiver, split_by))

  # --- Split recursion ---
  if (!is.null(split_by)) {
    result <- lapply(split(df, df[[split_by]]),
                     df_to_mat,
                     sender = sender,
                     receiver = receiver,
                     perceiver = perceiver,
                     mode = mode,
                     loops = loops,
                     multi_mode = multi_mode)
    return(purrr::transpose(result))
  }

  if (multi_mode) {
    nodes_s <- unique(df[[sender]])
    nodes_r <- unique(df[[receiver]])
    nodes_p <- if (!is.null(perceiver)) unique(df[[perceiver]]) else NULL
  } else {
    nodes_s <- nodes_r <- nodes_p <- unique(c(
      df[[sender]], df[[receiver]],
      if (!is.null(perceiver)) df[[perceiver]]
    ))
  }
  n_s <- length(nodes_s)
  n_r <- length(nodes_r)
  n_p <- if (!is.null(perceiver)) length(nodes_p) else NULL

  expected <- if (mode == "undirected") {
    if (loops) n_s * (n_s + 1) / 2 else n_s * (n_s - 1) / 2
  } else {
    if (loops) n_s * n_r else n_s * n_r - min(n_s, n_r)
  }
  if (anyNA(df[var_names]) || nrow(df) != expected) {
    warning("Incomplete dyadic data: some cells will be NA.",
            "\nCheck the data and consider coding matrices manually.")
  }

  make_structure <- function(var) {
    if (is.null(perceiver)) {
      mat <- matrix(NA_real_, nrow = n_s, ncol = n_r,
                    dimnames = list(nodes_s, nodes_r))
      mat[cbind(df[[sender]], df[[receiver]])] <- df[[var]]
      if (mode == "undirected")
        mat[cbind(df[[receiver]], df[[sender]])] <- df[[var]]
      if (!loops) diag(mat) <- NA
      mat
    } else {
      arr <- array(NA_real_, dim = c(n_s, n_r, n_p),
                   dimnames = list(nodes_s, nodes_r, nodes_p))
      arr[cbind(df[[sender]], df[[receiver]], df[[perceiver]])] <- df[[var]]
      if (mode == "undirected")
        arr[cbind(df[[receiver]], df[[sender]], df[[perceiver]])] <- df[[var]]
      if (!loops && !multi_mode)
        arr[cbind(nodes_s, nodes_s, rep(nodes_p, each = n_s))] <- NA
      arr
    }
  }

  setNames(lapply(var_names, make_structure), var_names)
}

