#' Internal Auxiliary function to fit the base model in QAPcss()
#'
#' @param mod formula; model to be estimated.
#' @param rand logical; are there any random intercepts requested?
#' @param family character; model family to be estimated.
#' @param pred matrix; matrix of predictors.
#' @param y numeric or character; dependent data vector
#' @param use_robust_errors logical; should standard errors be corrected by HC3?
#' @param nx numeric; number of predictors.
#'
#' @returns Returns fit for the baseline model
#' @import lme4
#' @import nnet

fit_base <- function(mod,
                     rand,
                     family,
                     pred,
                     nx,
                     y,
                     use_robust_errors) {
  if (use_robust_errors) {
    xv <- as.matrix(pred[,c(3:(2 + nx))])
  }
  fit <- list()
  if (!rand) {
    if (family == 'multinom') {
      base_model       <- nnet::multinom(mod, data = pred)
      fit$coefficients <- coefficients(base_model)
      fit$t  <- coefficients(base_model)/summary(base_model)$standard.errors
    } else {
      if (family == 'gaussian') {
        base_model        <- lm(mod, data = pred)
        fit$r.squared     <- summary(base_model)$r.squared
        fit$adj.r.squared <- summary(base_model)$adj.r.squared
      } else {
        base_model      <- glm(mod, data = pred, family = family)
      }
      fit$coefficients  <- coefficients(base_model)


      if (use_robust_errors) {
        fit$t <- fit$coefficients / HC3(xv,base_model$residuals)
      } else {
        fit$t <- summary(base_model)$coefficients[,3]
      }
    }
  } else {
    if (family == 'gaussian') {
      base_model <- lme4::lmer(mod, data = pred)
    } else {
      base_model  <- lme4::glmer(mod, data = pred, family = family,
                           control = glmerControl(calc.derivs = FALSE,
                                                  optimizer = "bobyqa"),
                           nAGQ = 0)
      fit$log_lik <- summary(base_model)[[6]]
    }
    fit$coefficients <- summary(base_model)$coefficients[,1]


    if (use_robust_errors) {
      fit$t <- fit$coefficients / HC3(xv,residuals(base_model))
    } else {
      fit$t <- summary(base_model)$coefficients[,3]
    }
    fit$random.intercepts <- list()
    for (rV in names(coefficients(base_model))) {
      fit$random.intercepts[[rV]] <- coefficients(base_model)[[rV]][,1]
    }
  }
  fit$base_model <- base_model
  return(fit)
}
