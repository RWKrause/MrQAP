#' Internal auxiliary function to make QAPglm() data
#'
#' @param y matrix; same as in \code{QAPglm()}
#' @param x matrix or list; same as in \code{QAPglm()}
#' @param g vector; created in same as in \code{QAPglm()} from groups
#' @param RIO matrix or list; same as \code{random_intercept_other} in \code{QAPglm()}
#' @param diag logical; same as in \code{QAPglm()}
#' @param mode character; same as in \code{QAPglm()}
#' @param net integer; internal parameter to create random intercept for each network in \code{y}
#' @param perm logical; should a permutation be performed
#' @param xi integer; either \code{NULL} or the number of the residualized variable in \code{x}
#'
#' @returns data in the format to be used in \code{QAPglm()} or \code{QAPglmPermEst}

make_qap_data <- function(y,
                          x,
                          g,
                          RIO,
                          diag,
                          mode,
                          net,
                          perm = FALSE,
                          xi = NULL) {
  nx <- length(x)

  if (perm && is.null(xi)) {
    y <- RMPerm(y, g)
  } else if (perm && !is.null(xi)) {
    x[[xi]] <- RMPerm(x[[xi]], g)
  }

  n <- dim(y)[1]
  valid <- matrix(TRUE,n,n)
  if (!diag) {
    diag(valid) <- FALSE
  }
  for (var in 1:nx) {
    valid[is.na(x[[var]])] <- FALSE
  }
  valid[is.na(y)] <- FALSE
  y[!valid] <- NA

  vv <- as.vector(valid)

  for (var in 1:nx) {
    x[[var]][!valid] <- NA
  }

  yv <- as.vector(y)[vv]

  sv <- matrix(1:n,n,n)
  sv[!valid] <- NA

  rv <- t(matrix(1:n,n,n))
  rv[!valid] <- NA

  pred <- data.frame(yv = yv,
                     nv = as.factor(net),
                     sv = as.vector(sv)[vv],
                     rv = as.vector(rv)[vv])
  pred$location <- as.vector(matrix(1:n**2,n,n))[vv]

  if (!is.null(RIO)) {
    if (!is.list(RIO)) {
      pred$ov <- as.factor(as.vector(RIO)[vv])
    } else {
      for (i in 1:length(RIO)) {
        pred[[paste0('ov',i)]] <- as.factor(as.vector(RIO[[i]])[vv])
      }
    }
  }

  for (var in c(1:nx)) {
    pred[[names(x)[var]]] <- as.vector(x[[var]])[vv]
  }
  return(pred)
}
