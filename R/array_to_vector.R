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