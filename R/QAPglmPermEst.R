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
#' @import fixest
#' @import dplyr
#' @import stringr
#' @import reformulas
#' @import fixest

QAPglmPermEst <- function(i,
                          data.,
                          mode.,
                          diag.,
                          mod.,
                          family.,
                          estimator.,
                          groups.,
                          fit.,
                          use_fixest.,
                          fixest_se_cluster.,
                          use_robust_errors.,
                          xi. = NULL,
                          rand.) {

  make_qap_data <- function(y,
                            x,
                            g,
                            diag,
                            mode,
                            net,
                            perm = FALSE,
                            xi = NULL) {

    RMPerm <- function(m, groups = NULL, CSS = FALSE) {

      if (is.list(m)) {
        return(lapply(m, RMPerm, groups = groups))
      }

      if (is.null(groups)) {
        groups <- rep(1,dim(m)[2])
      } else {
        groups <- as.character(groups)
      }

      if (length(dim(m)) == 2) {
        o <- unsplit(lapply(split(1:dim(m)[1],groups), FUN = sample),groups)
        p <- matrix(data = m[o, o], nrow = dim(m)[1], ncol = dim(m)[2])
      } else if (CSS) {
        p <- array(dim = c(dim(m)[1], dim(m)[2], dim(m)[3]))
        o <- unsplit(lapply(split(1:dim(m)[2],groups), FUN = sample),groups)
        p[, , ] <- array(m[o, o, o])
      } else {
        p <- array(dim = c(dim(m)[1], dim(m)[2], dim(m)[3]))
        for (i in 1:dim(m)[1]) {
          o <- unsplit(lapply(split(1:dim(m)[2],groups), FUN = sample),groups)
          p[i, , ] <- array(m[i, o, o])
        }
      }
      return(p)
    }


    nx <- length(x)
    if (perm && is.null(xi)) {
      y <- RMPerm(y, groups = g)
    } else if (perm && !is.null(xi)) {
      x[[xi]] <- RMPerm(x[[xi]], groups = g)
    }

    print(dim(y))
    print(colnames(y))
    print(all(is.na(y)))


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
    pred <- data.frame(location = as.vector(matrix(1:n**2,n,n))[vv],
                       yv = as.vector(y)[vv])
    print(dim(pred))
    print(13)
    pred$nv <- as.factor(net)
    print(14)
    sv <- matrix(1:n,n,n)
    print(15)
    sv[!valid] <- NA
    print(16)
    pred$sv <- as.vector(sv)[vv]
    print(17)
    rv <- t(matrix(1:n,n,n))
    print(18)
    rv[!valid] <- NA
    print(19)
    pred$rv <- as.vector(rv)[vv]
      print(20)
    for (var in c(1:nx)) {
      pred[[names(x)[var]]] <- as.vector(x[[var]])[vv]
    }
      print(21)
    return(pred)
  }

  dependent <- all.vars(mod.)[1]
  predictors <- all.vars(mod.[-1])
  if (!is.null(fixest_se_cluster.)) {
    predictors <- c(predictors, fixest_se_cluster.)
  }
  nx <- length(predictors)
  if (any(stringr::str_detect(as.character(mod.), '\\|'))) {
    if (any(stringr::str_detect(as.character(mod.), '\\('))) {
      main  <- all.vars(lme4::reformulas(mod.))[-1]
      fe_re <- reformulas::findbars(mod.) |>
        lapply(all.vars) |>
        unlist() |>
        unique()
    } else {
      main <- all.vars(mod.[[3]][[2]])
      fe_re <- all.vars(mod.[[3]][[3]])
    }
  } else {
    main <- predictors
    fe_re <- NULL
  }

  nx <- length(predictors)

  if (!is.list(data.[[dependent]])) {
    pred <- make_qap_data(y = data.[[dependent]],
                          x = data.[predictors],
                          g = groups.,
                          diag = diag.,
                          mode = mode.,
                          net = 1,
                          perm = TRUE,
                          xi = NULL)
  } else {
    pred_list <- vector(mode = 'list', length = length(data.[[dependent]]))


    for (net in 1:length(data.[[dependent]])) {
      x2 <- vector(mode = 'list', length = nx)
      for (var in 1:nx) {
        x2[[var]] <- data.[predictors][[var]][[net]]
      }

      names(x2) <- predictors
      pred_list[[net]] <- make_qap_data(y = data.[[dependent]][[net]],
                                        x = x2,
                                        g = groups.[[net]],
                                        diag = diag.,
                                        mode = mode.,
                                        net = net,
                                        perm = TRUE,
                                        xi = NULL)
    }

    pred <- Reduce(f = 'rbind', pred_list)
  }
  pred <- pred |> rename(!!sym(dependent) := yv)

  print(2)

  if (!rand.) {
    if (estimator. == 'standard') {
      if (!use_fixest.) {
        pm  <- glm(mod., data = pred, family = family.)
        if (use_robust_errors.) {
          pres <- rbind(pm$coefficients,
                        pm$coefficients / HC3(pred[main], residuals(pm)))
        } else {
          pres <- rbind(pm$coefficients,
                        summary(pm)$coefficients[,3])
        }
      } else {
        pm <- fixest::feglm(mod.,
                            data = pred,
                            family = "gaussian",
                            cluster = fixest_se_cluster.)

        if (use_robust_errors.) {
          pres <- rbind(c(Intercept = NA, pm$coefficients),
                        c(Intercept = NA, pm$coefficients) /
                          HC3(as.matrix(pred[main]),residuals(pm)))
        } else {
          pres <- rbind(c(Intercept = NA, pm$coefficients),
                        c(Intercept = NA, fit$coefficients[main] /
                            sqrt(diag(vcov(pm)))))
        }
      }
    } else {
      if (family. == 'binomial') {
        pm <- gmm::gmm(logit_moments,
                       x = list(y = pred$yv,
                                x = cbind(1,as.matrix(pred[,names(xRm.)]))),
                       t0 = rnorm(nx + 1),
                       wmatrix = "optimal",
                       vcov = "MDS",
                       optfct = "nlminb",
                       control = list(eval.max = 10000))
        resid <- logit_resid(pm)

      } else if (family. == 'poisson') {
        pm <- gmm::gmm(poisson_moments,
                       x = list(y = pred$yv,
                                x = cbind(1,as.matrix(pred[,names(xRm.)]))),
                       t0 = rnorm(nx + 1),
                       wmatrix = "optimal",
                       vcov = "MDS",
                       optfct = "nlminb",
                       control = list(eval.max = 10000))

        resid <- poisson_resid(pm)

      }

      if (use_robust_errors.) {
        pres <- rbind(pm$coefficients,
                      pm$coefficients / HC3(xv., resid))
      } else {
        pres <- rbind(pm$coefficients,
                      summary(pm)$coefficients[,3])
      }
      colnames(pres) <- c('Intercept',names(xRm.))
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

  print(3)

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
