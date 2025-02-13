#' Internal auxiliary function to transform an array into a vector
#'
#' @param ar array; The \code{array} to be transformed.
#' @param mode. character; Is the \code{array} \code{mode. = 'directed'} or \code{mode. = 'undirected'}?
#' @param diag. logical; Should the diagonal of each perceiver slice of the array be used?
#'
#' @returns Returns a vector

array_to_vector <- function(ar, mode., diag.) {
  v <- c()
  for (i in 1:nrow(ar)) {
    if (mode. == 'undirected') {
      v <- c(v,
             as.vector(ar[,,i][upper.tri(
               ar[,,i], diag = diag.)]))
    } else {
      v <- c(v,as.vector(ar[,,i]))
    }
  }
  return(v)
}