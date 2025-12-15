#' Internal auxiliary function to estimate one CSS permutation
#'
#' @param i NULL; iterator of \code{parLapply()}
#' @param y. array or list; same as \code{y} in \code{QAPcss()}
#' @param x. array or list; same as \code{x} in \code{QAPcss()}
#' @param mode. character; same as \code{mode} in \code{QAPcss()}
#' @param diag. logical; same as \code{diag} in \code{QAPcss()}
#' @param rand. logical; TRUE if any random intercepts are requested
#' @param family. character; same as \code{family} in \code{QAPcss()}
#' @param groups. vector; same as \code{groups} in \code{QAPcss()}
#' @param fit. glm, lm, lmer, glmer, or multinom; Estimate on unpermuted data from \code{fit_base}
#' @param comp. list; same as \code{comparison} in \code{QAPcss()}
#' @param RIO. array or list; same as \code{random_intercept_other} in \code{QAPcss()}
#' @param use_robust_errors. logical; same as \code{use_robust_errors} in \code{QAPcss()}
#' @param xi. integer; either \code{NULL} or the number of the residualized variable in \code{x}
#' @param xRm. array; residualized \code{array} for predictor \code{xi}
#' @param reference. numeric or character; same as \code{reference} in \code{QAPcss()}
#' @param mod. formula; model to be estimated.
#'
#' @returns Returns permutation results for one CSS permutation.
#' @import lme4
#' @import nnet
#' @import dplyr
#'
QAPcssPermEst <- function(i,
                          y.,
                          x.,
                          mode.,
                          diag.,
                          rand. = rand,
                          family.,
                          groups.,
                          fit.,
                          comp.,
                          RIO.,
                          use_robust_errors.,
                          xi. = NULL,
                          xRm. = NULL,
                          reference. = NULL,
                          mod.) {

  nx <- length(x.)
  sufficient_data <- FALSE
  trial <- 0

  if (!is.null(RIO.)) {
    rio <- TRUE
  } else {
    rio <- FALSE
  }

  y_cat <- na.omit(unique(as.vector(unlist(y.))))

  while (!sufficient_data && trial < 10000) {
    trial <- trial + 1

    if (!is.list(y.)) {
      if (is.null(xRm.)) {
        y. <- RMPerm(y., groups., CSS = TRUE)
      } else {
        x.[[xi.]] <- RMPerm(xRm., groups., CSS = TRUE)
      }

      pred <- make_css_data(y = y.,
                            x = x.,
                            nets = 1,
                            RIO = RIO.,
                            rio = rio,
                            diag = diag.,
                            mode = mode.)$pred
    } else {
      if (is.null(xRm.)) {
        y. <- lapply(y., RMPerm, CSS = TRUE)
      } else {
        x.[[xi.]] <- lapply(x.[[xi.]], RMPerm,CSS = TRUE)
      }

      pred_list <- vector(mode = 'list', length = length(y.))

      for (gr in 1:length(y.)) {
        xgr <- list()
        for (var in 1:nx) {
          xgr[[names(x.)[var]]] <- x.[[var]][[gr]]
        }
        pred_list[[gr]] <- make_css_data(y = y.[[gr]],
                                         x = xgr,
                                         nets = gr,
                                         RIO = RIO.[[gr]],
                                         diag = diag.,
                                         rio = rio,
                                         mode = mode.)$pred
      }

      pred <- Reduce(f = 'rbind', pred_list)
    }



    if (family. != 'multinom' && is.null(comp.)) {
      y_ok <- ifelse(length(na.omit(unique(pred$yv))) > 1, TRUE, FALSE)
    } else {
      y2_cat <- pred$yv |>
        unique() |>
        na.omit()
      y_present <- ifelse(all(y_cat %in% y2_cat), TRUE, FALSE)
      y_mult <- all(table(pred$yv) > 2)
      y_ok <- y_mult * y_present
    }

    x_ok <- TRUE

    x_ok <- pred |>
      select(names(pred)[names(pred) %in% names(x.)]) |>
      select(where(is.numeric)) |>
      apply(MARGIN = 2, FUN = unique)
    if (!is.list(x_ok)) {
      x_ok <- ifelse(length(x_ok) > 1, TRUE, FALSE)
    } else {
      x_ok <- x_ok |>
        lapply(FUN = length) |>
        unlist() > 1
      x_ok <- all(x_ok)
    }




    if (nrow(pred) != 0 && x_ok && y_ok && !is.null(comp.)) {
      pred2 <- pred
      for (k in 1:length(comp.)) {
        pred2 <- pred
        pred2 <- pred2[pred2$yv %in% comp.[[k]],]
        pred2$yv <- ifelse(pred2$yv == comp.[[k]][1],0,1)

        tryCatch(cor(pred2[,c('yv',names(x.))],
                     use = 'complete.obs'),
                 error = function(e) {print(pred2); stop('ahh')}) -> cors
        if (any(is.na(cors))) {
          y_ok <- x_ok <- FALSE
        } else {
          diag(cors) <- 0
          cors <- abs(cors)
          if (any(cors > 0.9999)) {
            print(cors)
            y_ok <- x_ok <- FALSE
          }
        }
      }
    }



    sufficient_data <- as.logical(y_ok * x_ok)
  }

  if (trial == 10000) {
    stop('Too much missing data or wrong data entry.\n',
         'Cannot find permutation with valid data.\n')
  }




  if (is.null(comp.)) {
    xresL <- fit_perm(family. = family.,
                      predx = pred,
                      ref = reference.,
                      mod = mod.,
                      fitx = fit.,
                      xi. = xi.,
                      nx. = nx,
                      rand = rand.,
                      use_robust_errors = use_robust_errors.)
  } else {
    xresL <- vector(mode = 'list', length = length(comp.))
    for (k in 1:length(comp.)) {
      predK <- pred[pred$yv %in% comp.[[k]],]
      predK$yv <- ifelse(predK$yv == comp.[[k]][1],0,1)

      xresL[[k]] <- fit_perm(family. = family.,
                             predx = predK,
                             ref = reference.,
                             mod = mod.,
                             xi. = xi.,
                             fitx = fit.[[k]],
                             nx. = nx,
                             rand = rand.,
                             use_robust_errors = use_robust_errors.)
    }
  }

  return(xresL)
}
