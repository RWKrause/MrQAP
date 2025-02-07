#' Perform HC3 correction for heteroskedasticity
#'
#' @param X matrix; Matrix of predictors
#' @param e numeric; vector of residuals
#'
#' @returns Returns adjusted covariance matrix.

HC3 <- function(X,e) {
  XO <- cbind(matrix(1,dim(X)[1],1),X)
  XTXINV <- solve(t(XO) %*% XO)
  om <- c()

  hf <- function(z,XTXINV = XTXINV) {
   h <- z %*% XTXINV %*% z
   return(h)
  }
  h  <- apply(XO,1,hf,XTXINV = XTXINV)
  om <- e^2/(1 - h)^2
  x <- sqrt(diag( t(t(XTXINV %*% t(XO)) * om) %*% XO %*% XTXINV))
  gc()
  return(x)
}
