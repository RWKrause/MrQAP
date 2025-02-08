#' Internal Auxiliary function to a permutation model for QAPcss()
#'
#' @param family. character; model family to be estimated.
#' @param predx data.frame; predictor matrix.
#' @param ref character; reference category if \code{family. == 'multinom'}.
#' @param xi. numeric; number of the predictor currently being permuted.
#' @param mod. formula; model to be estimated.
#' @param fitx list; output from \code{fit_base()}.
#' @param rand logical; are there any random intercepts requested?
#' @param use_robust_errors logical; should standard errors be corrected by HC3?
#' @param nx numeric; number of predictors.
#'
#' @returns Returns if a current permutation is more or less extreme than observed.
#' @import lme4
#' @import nnet


fit_perm <- function(family.,
                     predx,
                     ref,
                     xi.,
                     mod,
                     fitx,
                     rand,
                     use_robust_errors) {

  if (use_robust_errors) {
    xv <- as.matrix(pred[,c(3:(2 + nx))])
  }

  if (family. == 'multinom') {
    predx$yv <- as.factor(predx$yv)
    predx$yv <- relevel(predx$yv, ref = ref)

    pm   <- multinom(mod., data = predx, model = TRUE)
    pres <- rbind(coefficients(pm),
                  coefficients(pm) /
                    summary(pm)$standard.errors )


    presL <- list()



    presL$lower  <- pres <= rbind(fitx$coefficients,
                                  fitx$t)
    presL$larger <- pres >= rbind(fitx$coefficients,
                                  fitx$t)
    presL$abs <- abs(pres) >= rbind(abs(fitx$coefficients),
                                    abs(fitx$t))


  } else {
    if (!rand) {
      pm  <- glm(mod., data = predx, family = family.)

      if (use_robust_errors) {
        pres <- rbind(pm$coefficients,
                      pm$coefficients / HC3(xv,
                                            residuals(pm)))
      } else {
        pres <- rbind(pm$coefficients,
                      summary(pm)$coefficients[,3])
      }
    } else {
      if (family. == 'gaussian') {
        pm <- lmer(mod., data = predx)

      } else {
        pm <- glmer(mod., data = predx, family = family.)

      }
      if (use_robust_errors) {
        pres <- rbind(summary(pm)$coefficients[,1],
                      summary(pm)$coefficients[,1] /
                        HC3(xv,residuals(pm)))
      } else {
        pres  <- rbind(summary(pm)$coefficients[,1],
                       summary(pm)$coefficients[,3])
      }
    }
  }


  xresP <- list()
  if (is.null(xi.)) {

    xresP$lower  <- pres < rbind(fitx$coefficients,
                                 fitx$t)
    xresP$larger <- pres > rbind(fitx$coefficients,
                                 fitx$t)
    xresP$abs <- abs(pres) > rbind(abs(fitx$coefficients),
                                   abs(fitx$t))

  } else {
    xresP$lower  <- (pres < rbind(fitx$coefficients,
                                  fitx$t))[,xi.]
    xresP$larger <- (pres > rbind(fitx$coefficients,
                                  fitx$t))[,xi.]
    xresP$abs <- (abs(pres) > rbind(abs(fitx$coefficients),
                                    abs(fitx$t)))[,xi.]
  }

  return(xresP)
}
