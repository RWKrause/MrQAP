#' Parameter and p-value estimation with MRQAP for CSS data
#'
#' This function allows you to estimate \code{lm()}, \code{glm()}, \code{lmer()}, \code{glmer()}, and m\code{ultinom()} models on cubic CSS data. Data structure is expected to be \code{y[i,j,k]} with i = sender, j = receiver, and k = perceiver. If your data is differently structured, be careful when using the \code{random_intercept_...} arguments.
#'
#' Permutations are performed for all three dimensions simultaneously. If permutations should only be within perceiver slices, use \code{QAPglm()} instead.
#'
#' @param y array or list; \code{y} needs to be a cubic \code{array}. Alternatively, \code{y} can be a \code{list} of \code{array}s, if there are multiple CSS' that should be predicted at the same time.
#'
#' @param x array or list; needs to be a cubic \code{array} with the same dimensions as \code{y} and is the predictor for \code{y}. In most cases you have more than one predictor. Then \code{x} needs to be a \code{list} of \code{array}s of the same dimensionality as \code{y}. If \code{y} is a \code{list}, then \code{x} should be a \code{list} of \code{list}s. Each predictor variable should be its own \code{list}, with each entrance being an \code{array} for each of the separate \code{array}s in \code{y}. These \code{list}s are then combined into one \code{list} of \code{list}s (e.g., \code{x[[1]][[2]]} is the predictor \code{array} of the first predictor(\code{[[1]]}) for the second group (\code{[[2]]})). It is highly recommended that \code{x} is named. The names will be carried to the output.
#'
#' @param family character; While there is controversy around using anything but a linear model in MRQAP (family = 'gaussian'), \code{QAPglm()} supports all natural \code{R} model families (see \code{?family}). Use \code{family = 'multinom'} if a multinomial choice model is requested. Internally, estimation will be done with \code{multinom()} from the \code{nnet} package.
#'
#' @param mode character; indicates whether the \code{array} in \code{y} is 'directed' or 'undirected'. If it is 'undirected' only the upper triangle (\code{upper.tri()}) of each perceiver \code{matrix} will be used.
#'
#' @param diag logical; If \code{TRUE} diagonal values in each perceiver slice will also be included in the calculation. This is set to be \code{FALSE} by default and can potentially bias results when \code{TRUE}. (If \code{TRUE} is meaningful, it is recommended to run the model with and without \code{diag = TRUE} to see how relevant the diagonal is.)
#'
#' @param nullhyp character; Currently only two baseline models are available \code{nullhyp = 'qapy'} and \code{nullhyp = 'qapspp'}. In general, 'qapspp' is the recommended option (see Dekker, Krackhardt, & Snijders, 2007). However, it costs more time and if all \code{x} are uncorrelated with each other both 'qapy' and 'qapspp' will give the same results.
#'
#' @param reps integer; indicates how many permutations should be performed. Default is 1000 but larger numbers are highly recommended.
#'
#' @param seed integer; Given the random nature of the permutation, every call of \code{QAPglm()} will lead to different responses that will asymptotically converge with larger values for \code{reps}. To get consistent answers, one should specify a random number seed with the seed argument (e.g., \code{seed = 1402}).
#'
#' @param groups vector; It might be that a CSS is composed of qualitatively different groups. In that case, it might be desirable to only permute within groups. \code{groups} is a \code{vector} of \code{length = nrow(y)}, indicating the grouping of the nodes. (\code{groups} together with a \code{list} of \code{y} is not yet implemented but will be added soon.)
#'
#' @param ncores integer; QAPcss() is parallelized (using the \code{parallel} package). If multiple cores are available, using them cuts the estimation in near linear relation by using multiple cores (e.g., \code{ncores = 10}). Be aware that depending on the \code{R} installation, parallelization can fail when too many cores are addressed at the same time. If you are using an HPC cluster, the recommendation is to submit many jobs, each only asking for 5-20 cores, and running only a few hundred \code{reps}. The resulting outputs can then be combined with the auxiliary function \code{combine_qap_estimates()}.
#'
#' @param random_intercept_... logical; Multiple arguments exist to specify random intercepts. instead of relying on \code{glm()} estimation will use the \code{lmer()} or \code{glmer()} from the \code{lme4} package to obtain parameter estimates and t-values for the permutation assessment. \code{random_intercept_group} includes a random intercept for each group in groups. \code{random_intercept_sender}, \code{random_intercept_receiver}, and \code{random_intercept_perceiver} add intercepts for each node on the corresponding dimension. Additionally, there is \code{random_intercept_other}. This argument expects an input corresponding to every element in \code{y} (either a similarly sized \code{matrix} or a \code{list} of \code{array}s corresponding to each element in \code{y} when \code{y} is a \code{list}). This \code{array} (or \code{list} thereof) will be used to create cell-specifc random intercepts. This might be useful if some relationships are qualitatively different than others but all are expected to follow the same pattern (e.g., for some cells differences in (intercept and) residual variance are expected). You can submit \code{list}s of \code{list}s of matrices if you want multiple random intercept variables added
#'
#' @param use_robust_errors logical; indicates if internal standard errors should be adjusted for heteroskedasticity using the HC3 adjustment. This is by default \code{FALSE} but recommended when linear probability models are being used or heterogeneity is otherwise suspected. Results will overall be more conservative and thus significant findings more reliable.
#'
#' @param reference numeric or character; If \code{family = 'multinom'}, \code{reference} can be used to specify the reference group for the multinomial regression.
#'
#' @param comparison list; In case accuracy of perception is the outcome three different methods are supported. First, one can run a multinomial choice model comparing true positives (TP), false positives (FP), false negatives (FN), and true negatives (TN) - it is recommended to use TN as \code{reference}, because it is the most likely result in many graphs. The problem with this approach is that the comparisons TN vs TP and TN vs FN are not reasonably applicable, the alternative to a true negative is not at true positive, but if i and j have no relationship (Y_{ij} = 0) then there can only be a true negative (Y_{ijk} = 0) or a false positive (Y_{ijk} = 1). The same applies when Y_{ij} = 1, then only true positives and false negatives are meaningful comparisons. Thus, the second option is to run \code{QAPcss()} twice, once with only TN and FP, setting all TP and FN to \code{NA} in \code{y} and \code{x} (and another time with the reverse set to \code{NA}). If the corresponding values in \code{x} are not set to \code{NA} the result should be asymptotically the same as option three. The third option is to specify the \code{comparison}. Provide a (named) \code{list} of the two or more comparisons to be performed: \code{qapCSS(...,       comparison = list(comission = c('false_positive', 'true_negative'), omission  = c('false_negative', 'true_positive')))}. The difference to option two is that here permutations of the entire cube are performed before it is split into the two comparisons and the requested analyses (see \code{family}) are performed. The data are internally re-coded to 0 and 1, thus linear probability models (\code{family = 'gaussian'}) and logistic regression (\code{family = 'binomial'}) are possible. This means that a value in \code{x} that was linked to a the TN vs FP comparison can also predict values in the FN vs TP comparison. More research is necessary to identify when which method gives the most reliable results.
#'
#' @param error_file character; Not meant for the end-user. Passed to \code{makeCluster(ncores, outfile = error_file)}. See \code{?parallel::makeCluster}. Helpful if you want to debug the code.
#'
#' @returns an object of \code{class} \code{QAPRegression} when \code{family = 'gaussian'} or \code{QAPGLM} otherwise. It contains the basic input parameters and estimated parameters and p-values.
#' import parallel
#' import lme4
#' @export
#' @import parallel
#' @import lme4
#' @import nnet
#'

QAPcss <- function(y,
                   x,
                   mode = "directed",
                   diag = FALSE,
                   nullhyp = "qapy",
                   reps   = 1000,
                   seed   = NULL,
                   ncores = NULL,
                   family = 'gaussian',
                   groups = NULL,
                   reference = NULL,
                   comparison = NULL,
                   use_robust_errors = FALSE,
                   random_intercept_group = FALSE,
                   random_intercept_sender = FALSE,
                   random_intercept_receiver = FALSE,
                   random_intercept_perceiver = FALSE,
                   random_intercept_other = NULL,
                   error_file = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(ncores)) {
    ncores <- 1
  }

  if (!is.list(y)) {
    large <- FALSE
  } else {
    large <- TRUE
  }

  if (length(dim(y)) != 3 && !large) {
    stop("Wrong data entry!
         y must be a 3-dimensional array with the same length in each dimension!
         The expected format is:
         [sender, receiver, perceiver]")
  }

  if (large) {
    for (i in 1:length(y)) {
      if (length(dim(y[[i]])) != 3 && length(unique(dim(y[[i]]))) == 1) {
        stop("Wrong data entry in element ",i, " in y.",
             "\n each element in y ",
             "must be a 3-dimensional array with the same ",
             "length in each dimension!",
             "\nThe expected format is:
         [sender, receiver, perceiver]\n")
      }
    }
  }





  if (!is.list(x)) {
    x <- list(x)
  }


  nx <- length(x)

  rig <- random_intercept_group
  rip <- random_intercept_perceiver
  ris <- random_intercept_sender
  rir <- random_intercept_receiver
  RIO <- random_intercept_other
  if (!is.null(RIO)) {
    rio <- TRUE
  } else {
    rio <- FALSE
  }

  if (is.null(groups) && rig) {
    warning('Faulty argument!
            Random intercepts for groups requested but no groups given.
            Will proceed without random intercept estimation.')
    rig <- FALSE
  }

  rand <- any(c(rig, rip, ris, rir, rio))



  if (!large) {
    if (!is.array(y) || !as.logical(length(unique(dim(y))))) {
      stop("Wrong data entry!
         y must be an array with the same length in each dimension!
         The expected format is:
         [sender, receiver, perceiver]")
    }

    for (var in 1:nx) {
      if (any(dim(x[[var]]) != dim(y))) {
        stop('Wrong data entry!
           Not all arrays are the same size! Check variable ',
           var,
           ' in x.')
      }
    }

    n <- nrow(y)
    if (!is.null(groups)) {
      if (length(groups) != n) {
        stop('Wrong data entry!
           groups does not match N. Got ', length(groups), ' expected ',n,'.')
      }
      groups <- as.factor(groups)
    } else {
      groups <- as.factor(rep(1,n))
    }

    g <- array(NA, dim = c(n,n,n))
    for (i in 1:n) {
      g[,,i] <- groups[i]
    }

    sym <- c()
    for (i in 1:n) {
      sym <- c(isSymmetric(y[,,i]),sym)
    }
  } else {

    for (var in 1:nx) {
      for (gr in unique(groups)) {
        if (any(dim(x[[var]][[gr]]) != sum(groups == gr))) {
          stop('Wrong data entry in x ', var,' for group ', gr,
               "\n each element in x must be a list containing ",
               "a 3-dimensional array for each group.",
               "\n each of these arrays have the same dimensions as the ",
               "corresponding element in y.\n")
        }
      }
    }
    g <- NULL

    groups <- c()
    for (gr in 1:length(y)) {
      groups <- c(groups, rep(gr,dim(y[[gr]])[1]))
    }

    n <- length(groups)

    sym <- c()
    for (gr in unique(groups)) {
      for (i in 1:dim(y[[gr]])[1]) {
        sym <- c(isSymmetric(y[[gr]][,,i]),sym)
      }
    }
  }


  if (rand && family == 'multinom') {
    warning('Faulty argument!
            Random intercepts requested for multinomial choice.
            This is not implemented.
            Estimation will continue with standard nnet::multinom().')
    rig <- FALSE
    rip <- FALSE
    ris <- FALSE
    rir <- FALSE
    rio <- FALSE
    rand <- FALSE
  }


  if (!is.null(reference) &&
      !is.character(reference) &&
      family == 'multinom') {
    reference <- as.character(reference)
  }

  if (use_robust_errors && family == 'multinom') {
    warning('Faulty argument!
            Robust standard errors requested for multinomial choice.
            This is not implemented.
            Estimation will continue with standard nnet::multinom().')
    use_robust_errors <- FALSE
  }


  if ((nullhyp == "qapspp") && (nx == 1)) {
    nullhyp <- "qapy"
  }

  if (!(nullhyp %in% c("qapy","qapspp"))) {
    stop('Faulty argument!
         nullhyp does not match valid input.
         Got ', nullhyp, ' expected "qapy" or "qapspp".')
  }



  if (all(sym) && mode == 'directed') {
    warning('Mismatch between arguments and data.\n',
            ' mode = directed but y is symmetric for every perceiver.\n',
            ' Maybe y should be treated as undirected?\n')
  }

  if (!all(sym) && mode == 'undirected') {
    warning('Mismatch between arguments and data.\n',
            ' mode is undirected but y is not symmetric for every perceiver.\n',
            ' y will be treated as undirected.\n',
            'The upper triangle will be used - tri.upper().')
    mode <- 'undirected'
  }

  if (mode == 'undirected' && (ris || rir)) {
    warning('Mismatch between arguments/data.\n',
            'y is undirected and random intercepts for sender and/or receiver',
            ' requested.\n This is not possible.\n',
            'Random intercepts for sender and receiver will be set to FALSE')
    ris <- FALSE
    rir <- FALSE
  }

  if (diag) {
    warning('Results may not be valid when diagonal is used.\n')
  }


  char <- c()
  for (var in 1:nx) {
    char <- c(char,!is.numeric(x[[var]][[1]]))
  }


  if (all(char) && nullhyp == 'qapspp') {
    nullhyp <- 'qapy'
    cat('All predictors are characters/factors.\n',
          '"qapspp" is not implemented for this case.\n',
          'Using "qapy" instead.\n',
          'Maybe create separate dummy variables.\n\n')
  }


  if (is.null(names(x))) {
    warning('x is not named. Consider naming it...')
    names(x) <- paste0('x',1:nx)
  }



  if (!large) {
    cssd <- make_css_data(y = y,
                          x = x ,
                          g = g ,
                          RIO = RIO,
                          diag = diag,
                          mode = mode)
    pred <- cssd$pred
    valid <- cssd$valid
  } else {
    pred_list <- vector(mode = 'list', length = length(y))
    valid_list <- vector(mode = 'list', length = length(y))

    for (gr in 1:length(y)) {
      xgr <- list()
      for (var in 1:nx) {
        xgr[[names(x)[var]]] <- x[[var]][[gr]]
      }
      pred_list[[gr]] <- make_css_data(y = y[[gr]],
                                       x = xgr,
                                       g = array(gr, dim = dim(y[[gr]])),
                                       RIO = RIO[[gr]],
                                       diag = diag,
                                       mode = mode)$pred

      valid_list[[gr]] <- make_css_data(y = y[[gr]],
                                       x = xgr,
                                       g = array(gr, dim = dim(y[[gr]])),
                                       RIO = RIO[[gr]],
                                       diag = diag,
                                       mode = mode)$valid
    }

    pred <- Reduce(f = 'rbind', pred_list)
  }


  mod <- 'yv ~ 1'

  for (var in names(x)) {
    mod <- paste(mod, var, sep = ' + ')
  }

  rand_int <- ''

  if (rig) {
    rand_int <- paste(rand_int,'+ (1|gv)')
  }

  if (rip) {
    rand_int <- paste(rand_int,'+ (1|pv)')
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

  if (family == 'multinom') {
    pred$yv <- as.factor(pred$yv)

    if (is.null(reference)) {
      warning('No reference group provided for multinomial model.\n',
              'Reference group will be set to first:',
              levels(yv)[1])
    } else {
      pred$yv <- relevel(pred$yv, ref = reference)
    }
  }


  mod <- as.formula(mod)

  # baseline estimate

  rand <- any(c(rig,ris,rir,rip,rio))

  if (is.null(comparison)) {
    fit <- fit_base(mod = mod,
                    rand = rand,
                    family = family,
                    pred = pred,
                    nx = nx,
                    y = y,
                    use_robust_errors = use_robust_errors)
  } else {
    fit <- vector(mode = 'list', length = length(comparison))
    names(fit) <- names(comparison)
    for (k in 1:length(comparison)) {
      predK <- pred[pred$yv %in% comparison[[k]],]
      predK$yv <- ifelse(predK$yv == comparison[[k]][1],0,1)

      fit[[k]] <- fit_base(mod = mod,
                           rand = rand,
                           family = family,
                           pred = predK,
                           nx = nx,
                           y = y,
                           use_robust_errors = use_robust_errors)
    }
  }



  clust <- makeCluster(ncores, outfile = error_file)


  if (nullhyp == "qapy") {
    res <- parLapply(cl = clust, 1:reps,
                     fun = QAPcssPermEst,
                     y. = y,
                     x. = x,
                     g. = g,
                     mode. = mode,
                     diag. = diag,
                     rand. = rand,
                     family. = family,
                     groups. = groups,
                     fit. = fit,
                     comp. = comparison,
                     RIO. = RIO,
                     use_robust_errors. = use_robust_errors,
                     xi. = NULL,
                     xRm. = NULL,
                     reference. = reference,
                     mod. = mod)


    if (is.null(comparison)) {
      resL <- unlist(res, recursive = FALSE)

      fit$lower  <- Reduce(f = '+', resL[names(resL) == 'lower'],  0)/reps
      fit$larger <- Reduce(f = '+', resL[names(resL) == 'larger'], 0)/reps
      fit$abs    <- Reduce(f = '+', resL[names(resL) == 'abs'],    0)/reps
    } else {
      for (k in 1:length(comparison)) {
        resL <- unlist(res[[k]], recursive = FALSE)

        fit[[k]]$lower  <- Reduce(f = '+', resL[names(resL) == 'lower'], 0)/reps
        fit[[k]]$larger <- Reduce(f = '+', resL[names(resL) == 'larger'],0)/reps
        fit[[k]]$abs    <- Reduce(f = '+', resL[names(resL) == 'abs'],   0)/reps
      }
    }




  } else if (nullhyp == "qapspp") {

    if (family != 'multinom' && is.null(comparison)) {
      fit$lower  <- matrix(NA, nrow = 2,
                           ncol = length(fit$coefficients))
      colnames(fit$lower) <- names(fit$coefficients)

      fit$larger <- fit$abs <- fit$lower
    } else if (!is.null(comparison)) {
      for (k in 1:length(comparison)) {
        fit[[k]]$lower  <- matrix(NA, nrow = 2,
                                  ncol = length(fit[[k]]$coefficients))
        colnames(fit[[k]]$lower) <- names(fit[[k]]$coefficients)

        fit[[k]]$larger <- fit[[k]]$abs <- fit[[k]]$lower
      }

    } else {
      if (large) {
        ncat <- length(na.omit(unique(as.vector(unlist(y)))))
      } else {
        ncat <- length(na.omit(unique(as.vector(y))))
      }
      fit$lower  <- matrix(NA,
                           nrow = 2 * (ncat - 1),
                           ncol = ncol(coefficients(fit$base_model)))
      colnames(fit$lower) <- names(fit$coefficients)
      fit$larger <- fit$abs <- fit$lower
    }



    for (xi in names(x)[!char]) {
      modx <- paste(xi,'~ 1')
      for (varx in names(x)[names(x) != xi]) {
        modx <- paste(modx,varx, sep = ' + ')
      }
      if (!rand) {
        modx <- as.formula(modx)
        xm <- lm(modx, data = pred)
        xR <- xm$residuals
      } else {
        modx <- paste(modx, rand_int)
        modx <- as.formula(modx)

        xm <- lmer(modx, data = pred)
      }

      if (!large) {
        xR <- residuals(xm)

        xRm <- array(NA, dim = c(n,n,n))
        xRm[valid] <- xR
      } else {
        xRm <- x
        for (gr in unique(groups)) {
          xRm[[xi]][[gr]] <- array(NA, dim = c(sum(groups == gr),
                                               sum(groups == gr),
                                               sum(groups == gr)))
          xRm[[xi]][[gr]][valid_list[[gr]]] <- residuals(xm)[pred$gv == gr]
        }
      }

      res <- parLapply(cl = clust, 1:reps,
                       fun = QAPcssPermEst,
                       y. = y,
                       x. = x,
                       xi. = xi,
                       xRm. = xRm,
                       reference. = reference,
                       mode. = mode,
                       diag. = diag,
                       rand. = rand,
                       groups. = groups,
                       g. = g,
                       fit. = fit,
                       comp. = comparison,
                       family. = family,
                       RIO. = RIO,
                       use_robust_errors. = use_robust_errors,
                       mod. = mod)




      if (is.null(comparison)) {
        resL <- unlist(res, recursive = FALSE)

        fit$lower  <- Reduce(f = '+', resL[names(resL) == 'lower'],  0)/reps
        fit$larger <- Reduce(f = '+', resL[names(resL) == 'larger'], 0)/reps
        fit$abs    <- Reduce(f = '+', resL[names(resL) == 'abs'],    0)/reps
      } else {
        for (k in 1:length(comparison)) {

          resLL <- lapply(res, function(x){x[[k]]})

          resL <- unlist(resLL, recursive = FALSE)

          fit[[k]]$lower[,xi]  <- (Reduce(f = '+',
                                          resL[names(resL) == 'lower'],
                                          0)/reps)[,xi]
          fit[[k]]$larger[,xi] <- (Reduce(f = '+',
                                          resL[names(resL) == 'larger'],
                                          0)/reps)[,xi]
          fit[[k]]$abs[,xi]    <- (Reduce(f = '+',
                                          resL[names(resL) == 'abs'],
                                          0)/reps)[,xi]
        }
      }
    }

    if (any(char)) {
      res <- parLapply(cl = clust, 1:reps,
                       fun = QAPcssPermEst,
                       y. = y,
                       x. = x,
                       g. = g,
                       mode. = mode,
                       diag. = diag,
                       rand. = rand,
                       family. = family,
                       groups. = groups,
                       fit. = fit,
                       comp. = comparison,
                       RIO. = RIO,
                       use_robust_errors. = use_robust_errors,
                       xi. = NULL,
                       xRm. = NULL,
                       reference. = reference,
                       mod. = mod)



      if (is.null(comparison)) {
        resL <- unlist(res, recursive = FALSE)

        charV <- is.na( fit$lower[1,])
        charV[1] <- FALSE

        fit$lower[,charV]  <- (Reduce(f = '+', resL[names(resL) == 'lower'],
                                      0)/reps)[,charV]
        fit$larger[,charV] <- (Reduce(f = '+', resL[names(resL) == 'larger'],
                                      0)/reps)[,charV]
        fit$abs[,charV]    <- (Reduce(f = '+', resL[names(resL) == 'abs'],
                                      0)/reps)[,charV]
      } else {
        for (k in 1:length(comparison)) {

          resL <- unlist(res[[k]], recursive = FALSE)

          charV <- is.na( fit[[k]]$lower[1,])
          charV[1] <- FALSE

          fit$lower[[k]][,charV]  <- (Reduce(f = '+',
                                             resL[names(resL) == 'lower'],
                                             0)/reps)[,charV]
          fit$larger[[k]][,charV] <- (Reduce(f = '+',
                                             resL[names(resL) == 'larger'],
                                             0)/reps)[,charV]
          fit$abs[[k]][,charV]    <- (Reduce(f = '+',
                                             resL[names(resL) == 'abs'],
                                             0)/reps)[,charV]

        }
      }
    }
  }

  stopCluster(clust)

  fit$nullhyp   <- nullhyp
  fit$family    <- family
  fit$groups    <- unique(groups)
  fit$nullhyp   <- nullhyp
  fit$diag      <- diag
  fit$mode      <- mode
  fit$reps      <- reps
  fit$reference <- reference
  fit$comp      <- comparison
  fit$groups    <- groups
  fit$random    <- c(sender = ris,
                     receiver = rir,
                     perceiver = rip,
                     group = rig,
                     other = rio)
  fit$robust_se <- use_robust_errors
  if (family == 'multinom') {
    names(fit)[which(names(fit) == 't')] <- 'z'
  }


  class(fit) <- "QAPCSS"
  return(fit)
}
