#' Parameter and p-value estimation with MRQAP
#'
#' @param y matrix or list; \code{y} needs to be a square \code{matrix}. Alternatively, \code{y} can be a \code{list} of matrices, if there are multiple networks that should be predicted at the same time.
#'
#' @param x matrix or list; needs to be a square \code{matrix} with the same dimensions as \code{y} and is the predictor for \code{y}. In most cases you have more than one predictor. Then \code{x} needs to be a \code{list} of matrices of the same dimensionality as \code{y}. If \code{y} is a \code{list}, then \code{x} should be a \code{list} of \code{list}s. Each predictor variable should be its own \code{list}, with each entrance being a \code{matrix} for each of the separate matrices in \code{y}. These \code{list}s are then combined into one \code{list} of \code{list}s (e.g., \code{x[[1]][[2]]} is the predictor array of the first predictor for the second group). It is highly recommended that \code{x} is named. The names will be carried to the output.
#'
#' @param family character; While there is controversy around using anything but a linear model in MRQAP (family = 'gaussian'), \code{QAPglm()} supports all natural \code{R} model families (see \code{?family}).
#'
#' @param mode character; indicates whether the \code{matrix} in \code{y} is 'directed' or 'undirected'. If it is 'undirected' only the upper triangle (\code{upper.tri()}) of each \code{matrix} will be used.
#'
#' @param diag logical; If \code{TRUE} diagonal values will also be included in the calculation. This is set to be \code{FALSE} by default and can potentially bias results when \code{TRUE}. (If \code{TRUE} is meaningful, it is recommended to run the model with and without \code{diag = TRUE} to see how relevant the diagonal is.)
#'
#' @param nullhyp character; Currently only two baseline models are available \code{nullhyp = 'qapy'} and \code{nullhyp = 'qapspp'}. In general, 'qapspp' is the recommended option (see Dekker, Krackhardt, & Snijders, 2007). However, it costs more time and if all \code{x} are uncorrelated with each other both 'qapy' and 'qapspp' will give the same results.
#'
#' @param reps integer; indicates how many permutations should be performed. Default is 1000 but larger numbers are highly recommended.
#'
#' @param seed integer; Given the random nature of the permutation, every call of \code{QAPglm()} will lead to different responses that will asymptotically converge with larger values for \code{reps}. To get consistent answers, one should specify a random number seed with the seed argument (e.g., \code{seed = 1402}).
#'
#' @param groups vector; It might be that a larger network is composed of qualitatively different groups. In that case, it might be desirable to only permute within groups. \code{groups} is a \code{vector} of \code{length = nrow(y)}, indicating the grouping of the nodes. If \code{y} is a \code{list} then \code{groups} needs to be a \code{list} of vectors.
#'
#' @param ncores integer; QAPglm() is parallelized (using the \code{parallel} package). If multiple cores are available, using them cuts the estimation in near linear relation by using multiple cores (e.g., \code{ncores = 10}). Be aware that depending on the \code{R} installation, parallelization can fail when too many cores are addressed at the same time. If you are using an HPC cluster, the recommendation is to submit many jobs, each only asking for 5-20 cores, and running only a few hundred \code{reps}. The resulting outputs can then be combined with the auxiliary function \code{combine_qap_estimates()}.
#'
#' @param same_x_4_all_y logical; If \code{y} is a \code{list} but all \code{x} should remain the same for all \code{y}, toggle this to TRUE. This might be relevant if similar but different networks between the same nodes are predicted. Should probably be combined with \code{random_intercept_groups = TRUE}.
#'
#' @param random_intercept_... logical; Multiple arguments exist to specify random intercepts. instead of relying on \code{glm()} estimation will use the \code{lmer()} or \code{glmer()} from the \code{lme4} package to obtain parameter estimates and t-values for the permutation assessment. \code{random_intercept_nets} includes a random intercept for each network when \code{y} is a \code{list}. \code{random_intercept_sender} and \code{random_intercept_receiver} add intercepts for each node on the corresponding dimension (row = sender, column = receiver). Additionally, there is \code{random_intercept_other}. This argument expects an input corresponding to every element in \code{y} (either a similarly sized \code{matrix} or a \code{list} of matrices corresponding to each element in \code{y} when is a \code{list}). This \code{matrix} (or \code{list} thereof) will be used to create cell-specifc random intercepts. This might be useful if some relationships are qualitatively different than others but all are expected to follow the same pattern (e.g., for some cells differences in (intercept and) residual variance are expected). You can submit \code{list}s of \code{list}s of matrices if you want multiple random intercept variables added
#'
#' @param use_robust_errors logical; indicates if internal standard errors should be adjusted for heteroskedasticity using the HC3 adjustment. This is by default \code{FALSE} but recommended when linear probability models are being used or heterogeneity is otherwise suspected. Results will overall be more conservative and thus significant findings more reliable.
#'
#' @param error_file character; Not meant for the end-user. Passed to \code{makeCluster(ncores, outfile = error_file)}. See \code{?parallel::makeCluster}. Helpful if you want to debug the code.
#'
#' @returns an object of \code{class} \code{QAPRegression} when \code{family = 'gaussian'} or \code{QAPGLM} otherwise. It contains the basic input parameters and estimated parameters and p-values.
#' import parallel
#' import lme4
#' @export
#' @import parallel
#' @import lme4


QAPglm <- function(y,
                    x,
                    family = 'gaussian',
                    mode = "directed",
                    diag = FALSE,
                    nullhyp = "qapspp",
                    reps = 1000,
                    seed = NULL,
                    groups = NULL,
                    ncores = NULL,
                    same_x_4_all_y = FALSE,
                    random_intercept_nets   = FALSE,
                    random_intercept_sender   = FALSE,
                    random_intercept_receiver = FALSE,
                    random_intercept_other = NULL,
                    use_robust_errors = FALSE,
                    error_file = NULL) {

  if (!is.list(y)) {
    large <- FALSE
  } else {
    large <- TRUE
  }


  rin <- random_intercept_nets
  ris <- random_intercept_sender
  rir <- random_intercept_receiver
  RIO <- random_intercept_other
  if (!is.null(RIO)) {
    rio <- TRUE
  } else {
    rio <- FALSE
  }


  if (mode == 'directed') {
    mode = 'digraph'
  }
  if (mode == 'undirected') {
    mode = 'graph'
  }


  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(ncores)) {
    ncores <- 1
  }
  clust <- parallel::makeCluster(ncores, outfile = error_file)


  if (is.list(x)) {
    nx <- length(x)
  } else {
    x <- list(x)
    nx <- 1
  }


  if (!is.list(y)) {
    pred <- make_qap_data(y = y,
                          x = x,
                          g = groups,
                          RIO = RIO,
                          diag = diag,
                          mode = mode,
                          net = 1,
                          perm = FALSE,
                          xi = NULL)
  } else {
    pred_list <- vector(mode = 'list', length = length(y))
    for (net in 1:length(y)) {
      if (!same_x_4_all_y) {
        x2 <- vector(mode = 'list', length = nx)
        for (var in 1:nx) {
          x2[[var]] <- x[[var]][[net]]
        }
        if (is.list(RIO)) {
          RIO2 <- vector(mode = 'list', length = nx)

          for (ri in 1:length) {
            RIO2[[ri]] <- RIO[[ri]][[net]]
          }
        } else {
          RIO2 <- RIO
        }
      } else {
        x2 <- x
      }
      names(x2) <- names(x)
      pred_list[[net]] <- make_qap_data(y = y[[net]],
                                        x = x2,
                                        g = groups[[net]],
                                        RIO = RIO2,
                                        diag = diag,
                                        mode = mode,
                                        net = net,
                                        perm = FALSE,
                                        xi = NULL)
    }

    pred <- Reduce(f = 'rbind', pred_list)
  }



  mod <- 'yv ~ 1'

  for (var in names(x)) {
    mod <- paste(mod, var, sep = ' + ')
  }

  rand_int <- ''

  if (rin) {
    rand_int <- paste(rand_int,'+ (1|nv)')
  }

  if (ris) {
    rand_int <- paste(rand_int,'+ (1|sv)')
  }

  if (rir) {
    rand_int <- paste(rand_int,'+ (1|rv)')
  }


  if (rio) {
    if (!is.list(RIO)) {
      rand_int <- paste(rand_int,'+ (1|ov)')
    } else {
      for (i in 1:length(RIO)) {
        rand_int <- paste0(rand_int,' + (1|ov',i,')')
      }
    }
  }

  mod <- paste(mod,rand_int)

  mod <- as.formula(mod)
  fit <- list()

  rand <- any(c(rin, ris, rir,  rio))


  # baseline estimate


  if (!rand) {
    if (family == 'gaussian') {
      base_model        <- lm(mod, data = pred)
      fit$r.squared     <- summary(base_model)$r.squared
      fit$adj.r.squared <- summary(base_model)$adj.r.squared
    } else {
      base_model      <- glm(mod, data = pred, family = family)
    }
    fit$coefficients  <- base_model$coefficients
    if (use_robust_errors) {
      fit$t <- fit$coefficients / HC3(xv,fit$residuals)
    } else {
      fit$t <- summary(base_model)$coefficients[,3]
    }


  } else {
    if (family == 'gaussian') {
      base_model <- lme4::lmer(mod, data = pred)


    } else {
      base_model  <- lme4::glmer(mod, data = pred, family = family)
      fit$log_lik <- summary(base_model)[[6]]
    }
    fit$coefficients  <- summary(base_model)$coefficients[,1]

    if (use_robust_errors) {
      fit$t <- fit$coefficients / HC3(xv,fit$residuals)
    } else {
      fit$t <- summary(base_model)$coefficients[,3]
    }

  }





  if ((nullhyp == "qapspp") && (nx == 1)) {
    nullhyp <- "qapy"
  }


  if (nullhyp == "qapy") {
    res <- parallel::parLapply(cl = clust, 1:reps,
                     fun = QAPglmPermEst,
                     y. = y,
                     xRm. = x,
                     xi. = NULL,
                     mode. = mode,
                     diag. = diag,
                     mod. = mod,
                     groups. = groups,
                     fit. = fit,
                     family. = family,
                     RIO. = RIO,
                     use_robust_errors. = use_robust_errors,
                     same_x_4_all_y. = same_x_4_all_y,
                     rand. = rand)

    resL <- unlist(res,recursive = FALSE)

    fit$lower  <- Reduce(f = '+', resL[names(resL) == 'lower'],0)/reps
    fit$larger <- Reduce(f = '+', resL[names(resL) == 'larger'],0)/reps
    fit$abs    <- Reduce(f = '+', resL[names(resL) == 'abs'],0)/reps



  } else if (nullhyp == "qapspp") {

    fit$lower  <- matrix(NA, nrow = 2, ncol = (nx + 1))
    fit$larger <- matrix(NA, nrow = 2, ncol = (nx + 1))
    fit$abs    <- matrix(NA, nrow = 2, ncol = (nx + 1))
    colnames(fit$lower) <- names(fit$coefficients)
    colnames(fit$larger) <- names(fit$coefficients)
    colnames(fit$abs) <- names(fit$coefficients)

    for (xi in names(x)) {
      modx <- paste(xi,'~ 1')
      for (varx in names(x)[names(x) != xi]) {
        modx <- paste(modx, varx, sep = ' + ')
      }
      if (!rand) {
        modx <- as.formula(modx)
        xm <- lm(modx, data = pred)
      } else {
        modx <- paste(modx, rand_int)
        modx <- as.formula(modx)
        xm <- lmer(modx, data = pred)
      }

      xR <- residuals(xm)
      xRm <- x

      if (!large) {
        xRm[[xi]][pred$location] <- xR
      } else {
        for (net in 1:length(y)) {
          xRm[[xi]][[net]][pred$location[
            as.character(pred$r_net) == net]] <- xR[
              as.character(pred$r_net) == net]
        }
      }


      res <- parallel::parLapply(cl = clust, 1:reps,
                       fun = QAPglmPermEst,
                       y. = y,
                       xi. = xi,
                       mode. = mode,
                       diag. = diag,
                       mod. = mod,
                       family. = family,
                       groups. = groups,
                       fit. = fit,
                       xRm. = xRm,
                       RIO. = RIO,
                       use_robust_errors. = use_robust_errors,
                       same_x_4_all_y. = same_x_4_all_y,
                       rand. = rand)

      resL <- unlist(res, recursive = FALSE)

      fit$lower[,xi]  <- Reduce(f = '+',
                                    resL[names(resL) == 'lower'],
                                    0)/reps
      fit$larger[,xi] <- Reduce(f = '+',
                                    resL[names(resL) == 'larger'],
                                    0)/reps
      fit$abs[,xi]    <- Reduce(f = '+',
                                    resL[names(resL) == 'abs'],
                                    0)/reps
    }



  }

  stopCluster(clust)
  if (!is.null(names(x))) {
    names(fit$coefficients) <- c('(Intercept)', names(x))
    names(fit$t)            <- c('(Intercept)', names(x))
    colnames(fit$lower)     <- c('(Intercept)', names(x))
    colnames(fit$larger)    <- c('(Intercept)', names(x))
    colnames(fit$abs)       <- c('(Intercept)', names(x))
  }



  fit$nullhyp <- nullhyp
  fit$diag <- diag
  fit$family <- family

  if (mode == 'digraph') {
    mode = 'directed'
  }
  if (mode == 'graph') {
    mode = 'undirected'
  }
  fit$mode <- mode
  fit$reps <- reps
  fit$groups <- unique(unlist(groups))
  fit$simple_fit <- base_model
  fit$robust_se <- use_robust_errors

  if (family == 'gaussian') {
    class(fit) <- 'QAPRegression'
  } else {
    class(fit) <- 'QAPGLM'
  }
  return(fit)
}
