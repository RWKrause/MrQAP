#' Internal auxiliary function to estimate one permutation
#'
#' @param i NULL; iterator of parLapply
#' @param y. matrix or list; same as \code{y} in \code{QAPglm()}
#' @param mode. character; same as \code{mode} in \code{QAPglm()}
#' @param diag. logical; same as \code{diag} in \code{QAPglm()}
#' @param mod. formula; model to be estimated.
#' @param family. character; same as \code{family} in \code{QAPglm()}
#' @param estimator character; same as \code{estimator} in \code{QAPglm()}
#' @param groups. vector; same as \code{groups} in \code{QAPglm()}
#' @param fit. glm, lm, lmer, or glmer; Estimate on unpermuted data.
#' @param RIO. matrix or list; same as \code{random_intercept_other} in \code{QAPglm()}
#' @param use_robust_errors. logical; same as \code{use_robust_errors} in \code{QAPglm()}
#' @param xi. integer; either \code{NULL} or the number of the residualized variable in \code{x}
#' @param xRm. matrix or list; same as \code{x} in \code{QAPglm()}
#' @param same_x_4_all_y. matrix or list; same as \code{y} in \code{QAPglm()}
#' @param rand. logical; TRUE if any random intercepts are requested
#'
#' @returns Results for one permutation.
#' @import lme4
#' @import gmm

QAPglmPermEst <- function(i,
                          y.,
                          mode.,
                          diag.,
                          mod.,
                          family.,
                          estimator.,
                          groups.,
                          fit.,
                          RIO.,
                          use_robust_errors.,
                          xi. = NULL,
                          xRm. = NULL,
                          same_x_4_all_y.,
                          rand.) {
  nx <- length(xRm.)

  if (!is.list(y.)) {
    pred <- make_qap_data(y = y.,
                          x = xRm.,
                          g = groups.,
                          RIO = RIO.,
                          diag = diag.,
                          mode = mode.,
                          net = 1,
                          perm = TRUE,
                          xi = xi.)

  } else {
    pred_list <- vector(mode = 'list', length = length(y.))
    for (net in 1:length(y.)) {
      if (!same_x_4_all_y.) {
        x2 <- vector(mode = 'list', length = nx)
        for (var in 1:nx) {
          x2[[var]] <- xRm.[[var]][[net]]
        }
        if (is.list(RIO.)) {
          RIO2 <- vector(mode = 'list', length = nx)

          for (ri in 1:length) {
            RIO2[[ri]] <- RIO.[[ri]][[net]]
          }
        } else {
          RIO2 <- RIO.
        }
      } else {
        x2 <- xRm.
      }
      names(x2) <- names(xRm.)

      pred_list[[net]] <- make_qap_data(y = y.[[net]],
                                        x = x2,
                                        g = groups.[[net]],
                                        RIO = RIO2,
                                        diag = diag.,
                                        mode = mode.,
                                        net = net,
                                        perm = TRUE,
                                        xi = xi.)
    }

    pred <- Reduce(f = 'rbind', pred_list)
  }

  xv. <- as.matrix(pred[,(ncol(pred) - nx + 1):ncol(pred)])

  if (!rand.) {
    if (estimator. == 'standard') {
      pm  <- glm(mod., data = pred, family = family.)
      if (use_robust_errors.) {
        pres <- rbind(pm$coefficients,
                      pm$coefficients / HC3(xv.,
                                            residuals(pm)))
      } else {
        pres <- rbind(pm$coefficients,
                      summary(pm)$coefficients[,3])
      }
    } else {
      if (family == 'binomial') {
        pm <- gmm(logit_moments,
                  x = list(y = pred$y,
                           x = pred[,names(x)]),
                  t0 = rnorm(nx + 1),
                  wmatrix = "optimal",
                  vcov = "MDS",
                  optfct = "nlminb",
                  control = list(eval.max = 10000))

        resid <- logit_resid(base_model)

      }
      if (family == 'poisson') {
        pm <- gmm(poisson_moments,
                          x = list(y = pred$y,
                                   x = pred[,names(x)]),
                          t0 = rnorm(nx + 1),
                          wmatrix = "optimal",
                          vcov = "MDS",
                          optfct = "nlminb",
                          control = list(eval.max = 10000))

        resid <- poisson_resid(base_model)
      }
      if (use_robust_errors.) {
        pres <- rbind(pm$coefficients,
                      pm$coefficients / HC3(xv.,resid))
      } else {
        pres <- rbind(pm$coefficients,
                      summary(pm)$coefficients[,3])
      }
    }

  } else {
    if (family. == 'gaussian') {
      pm <- lme4::lmer(mod., data = pred)

    } else {
      pm <- lme4::glmer(mod., data = pred, family = family.)

    }
    if (use_robust_errors.) {
      pres <- rbind(summary(pm)$coefficients[,1],
                    summary(pm)$coefficients[,1] /
                      HC3(xv.,residuals(pm)))
    } else {
      pres  <- rbind(summary(pm)$coefficients[,1],
                     summary(pm)$coefficients[,3])
    }
  }


  xresL <- list()
  if (is.null(xi.)) {
    xresL$lower  <- pres <= rbind(fit.$coefficients,
                                  fit.$t)
    xresL$larger <- pres >= rbind(fit.$coefficients,
                                  fit.$t)
    xresL$abs <- abs(pres) >= rbind(abs(fit.$coefficients),
                                    abs(fit.$t))
  } else {
    xresL$lower  <- (pres <= rbind(fit.$coefficients,
                                   fit.$t))[,xi.]
    xresL$larger <- (pres >= rbind(fit.$coefficients,
                                   fit.$t))[,xi.]
    xresL$abs <- (abs(pres) >= rbind(abs(fit.$coefficients),
                                     abs(fit.$t)))[,xi.]
  }
  return(xresL)
}
