#' Internal auxiliary function to shape data into internal QAPcss() format
#'
#' @param y array or list; see \code{QAPcss()}
#' @param x array or list; see \code{QAPcss()}
#' @param g character or numeric; see \code{QAPcss()}
#' @param RIO logical; see \code{QAPcss()}
#' @param diag logical; see \code{QAPcss()}
#' @param mode character; see \code{QAPcss()}
#'
#' @returns Returns vectorized CSS data

make_css_data <- function(y,
                          x,
                          g,
                          RIO,
                          diag,
                          mode) {
  n <- dim(y)[1]
  nx <- length(x)
  valid <- array(TRUE, dim = c(n,n,n))

  if (!diag) {
    for (i in 1:n) {
      diag(y[,,i]) <- NA
      for (var in 1:nx) {
        diag(x[[var]][,,i]) <- NA
      }
    }
  }

  valid[is.na(y)] <- FALSE

  for (var in 1:nx) {
    valid[is.na(x[[var]])] <- FALSE
  }


  y[!valid] <- NA

  for (var in 1:nx) {
    x[[var]][!valid] <- NA
  }


  vv <- array_to_vector(valid,
                        mode. = mode,
                        diag. = diag)
  yv <- array_to_vector(y,
                        mode. = mode,
                        diag. = diag)[vv]
  gv <- array_to_vector(g,
                        mode. = mode,
                        diag. = diag)[vv]


  pred <- data.frame(yv = yv,
                     gv = gv)


  for (var in c(1:nx)) {
    pred[[names(x)[var]]] <- array_to_vector(x[[var]],
                                             mode. = mode,
                                             diag. = diag)[vv]
  }


  per <- sen <- rec <- array(NA, dim = c(n,n,n))

  for (i in 1:n) {
    sen[i,,] <- i
    rec[,i,] <- i
    per[,,i] <- i
  }

  pred$sv <- as.factor(array_to_vector(sen,
                                       mode. = mode,
                                       diag. = diag)[vv])
  pred$rv <- as.factor(array_to_vector(rec,
                                       mode. = mode,
                                       diag. = diag)[vv])
  pred$pv <- as.factor(array_to_vector(per,
                                       mode. = mode,
                                       diag. = diag)[vv])

  if (!is.null(RIO)) {
    if (!is.list(RIO)) {
      pred$ov <- as.factor(array_to_vector(RIO,
                                           mode. = mode,
                                           diag. = diag)[vv])
    } else {
      for (i in 1:length(RIO)) {
        pred[[paste0('ov',i)]] <- as.factor(array_to_vector(RIO[[i]],
                                                            mode. = mode,
                                                            diag. = diag)[vv])
      }
    }
  }

  return(list(pred = pred, valid = valid))
}
