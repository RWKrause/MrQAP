#' Parameter and p-value estimation with MRQAP for CSS data
#'
#' Estimates regression models on cubic Cognitive Social Structure (CSS) data
#' (3D arrays of format \code{[sender, receiver, perceiver]}) with
#' permutation-based inference.
#'
#' @param formula A formula describing the model.
#'   The left-hand side is the dependent variable name in \code{data};
#'   the right-hand side lists predictor variable names.
#'   Use \code{|} (without parentheses) for fixest fixed effects.
#'   Use \code{(1|var)} for lme4-style random effects.
#'
#' @param data Named list; each element is a cubic \code{array} (or a
#'   \code{list} of arrays if there are multiple CSS').  Names must
#'   correspond to the variables in \code{formula}.
#'
#' @param family Character; \code{"gaussian"} (default), \code{"binomial"},
#'   \code{"poisson"}, \code{"negbin"} (negative binomial),
#'   \code{"zip"} (zero-inflated Poisson), or \code{"multinom"}.
#'
#' @param mode Character; \code{"directed"} (default) or \code{"undirected"}.
#' @param diag Logical; include diagonal values (default \code{FALSE}).
#'
#' @param nullhyp Character; \code{"qapy"} (default) or \code{"qapspp"}.
#'
#' @param estimator Character; \code{"standard"} (default) or \code{"gmm"}
#'   (for binomial, poisson, negbin, and zip families).
#'
#' @param reps Integer; number of permutations (default 1000).
#' @param seed Integer; optional random seed.
#' @param ncores Integer; cores for parallel processing via the future
#'   framework.
#' @param groups Vector or list; permutation grouping.
#' @param fixest_se_cluster Character; cluster variable for fixest.
#'
#' @param reference Character or numeric; reference group for multinomial.
#' @param comparison Named list of length-2 character vectors for accuracy
#'   comparisons (e.g., \code{list(commission = c("false_positive",
#'   "true_negative"))}).
#'
#' @param random_intercept_nets,random_intercept_sender,random_intercept_receiver,random_intercept_perceiver
#'   Logical; add random intercepts.
#'
#' @param use_robust_errors Logical; use HC3 correction.
#' @param use_gpu Logical; use GPU-accelerated batch OLS (gaussian only,
#'   requires \code{torch} package).
#'
#' @return An object of class \code{QAPCSS}.
#'
#' @export

QAPcss <- function(formula,
                   data,
                   mode      = "directed",
                   diag      = FALSE,
                   nullhyp   = "qapy",
                   reps      = 1000,
                   seed      = NULL,
                   ncores    = NULL,
                   family    = "gaussian",
                   estimator = "standard",
                   groups    = NULL,
                   fixest_se_cluster = NULL,
                   reference  = NULL,
                   comparison = NULL,
                   use_robust_errors = FALSE,
                   random_intercept_nets      = FALSE,
                   random_intercept_sender    = FALSE,
                   random_intercept_receiver  = FALSE,
                   random_intercept_perceiver = FALSE,
                   use_gpu    = FALSE) {

  if (!is.null(seed)) set.seed(seed)

  # --- parse formula ---
  parsed <- parse_qap_formula(formula, fixest_se_cluster)
  dep       <- parsed$dependent
  main      <- parsed$main
  # Filter out structural vars (sv, rv, nv, pv) that are auto-generated
  data_vars <- intersect(parsed$all_data_vars, names(data))
  nx        <- length(main)

  # --- validate ---
  validate_qap_input(data, parsed, css = TRUE)
  large <- is.list(data[[dep]])

  # --- dimension checks ---
  if (!large) {
    y <- data[[dep]]
    if (length(dim(y)) != 3)
      stop("data[['", dep, "']] must be a 3-dimensional array ",
           "[sender, receiver, perceiver].")
  } else {
    for (i in seq_along(data[[dep]])) {
      if (length(dim(data[[dep]][[i]])) != 3)
        stop("data[['", dep, "']][[", i, "]] must be a 3D array.")
    }
  }

  # --- random intercept flags ---
  rin <- random_intercept_nets
  rip <- random_intercept_perceiver
  ris <- random_intercept_sender
  rir <- random_intercept_receiver

  # --- build internal formula ---
  mod <- build_internal_formula(formula,
                                rin = rin, ris = ris, rir = rir,
                                rip = rip)
  mod_str <- paste(deparse(mod, width.cutoff = 500), collapse = " ")
  has_random <- grepl("\\(", mod_str) || parsed$has_random
  use_fixest <- parsed$use_fixest
  if (has_random && use_fixest) {
    warning("Cannot combine fixest FE and lme4 random effects. ",
            "Using lme4 only.")
    use_fixest <- FALSE
  }
  mod <- as.formula(mod_str)

  # --- warning checks ---
  if (has_random && family == "multinom") {
    warning("Random intercepts not implemented for multinomial. ",
            "Using standard nnet::multinom().")
    has_random <- FALSE
  }
  if (!is.null(reference) && !is.character(reference) && family == "multinom")
    reference <- as.character(reference)
  if (use_robust_errors && family == "multinom") {
    warning("Robust SEs not implemented for multinomial.")
    use_robust_errors <- FALSE
  }
  if ((nullhyp == "qapspp") && (nx == 1)) nullhyp <- "qapy"
  if (mode == "undirected" && (ris || rir)) {
    warning("Undirected mode: sender/receiver random intercepts set to FALSE.")
    ris <- rir <- FALSE
  }
  if (diag) warning("Results may not be valid when diagonal is used.")

  # Build random-intercept string for residualisation
  rand_part <- ""
  if (rin) rand_part <- paste(rand_part, "+ (1|nv)")
  if (rip) rand_part <- paste(rand_part, "+ (1|pv)")
  if (ris) rand_part <- paste(rand_part, "+ (1|sv)")
  if (rir) rand_part <- paste(rand_part, "+ (1|rv)")

  # --- groups ---
  if (!large) {
    n <- dim(data[[dep]])[1]
    if (!is.null(groups)) {
      if (length(groups) != n)
        stop("groups length (", length(groups), ") != N (", n, ").")
      groups <- as.factor(groups)
    } else {
      groups <- as.factor(rep(1, n))
    }
  }

  # --- build data frame ---
  valid <- NULL; valid_list <- NULL
  if (!large) {
    x_list <- lapply(data_vars, function(v) data[[v]])
    names(x_list) <- data_vars
    cssd  <- make_css_data(y = data[[dep]], x = x_list,
                           nets = 1,
                           diag = diag, mode = mode)
    pred  <- cssd$pred
    valid <- cssd$valid
  } else {
    pred_list  <- vector("list", length(data[[dep]]))
    valid_list <- vector("list", length(data[[dep]]))
    for (gr in seq_along(data[[dep]])) {
      xgr <- lapply(data_vars, function(v) data[[v]][[gr]])
      names(xgr) <- data_vars
      cssd <- make_css_data(y = data[[dep]][[gr]], x = xgr,
                            nets = gr,
                            diag = diag, mode = mode)
      pred_list[[gr]]  <- cssd$pred
      valid_list[[gr]] <- cssd$valid
    }
    pred <- do.call(rbind, pred_list)
  }

  # rename yv -> dependent name
  names(pred)[names(pred) == "yv"] <- dep

  # --- baseline fit ---
  fit <- list()

  if (is.null(comparison)) {
    fit$base <- fit_qap_model(mod          = mod,
                              pred         = pred,
                              family       = family,
                              estimator    = estimator,
                              use_fixest   = use_fixest,
                              fixest_se_cluster = fixest_se_cluster,
                              use_robust_errors = use_robust_errors,
                              main_vars    = main,
                              has_random   = has_random,
                              reference    = reference)
  } else {
    fit$base <- vector("list", length(comparison))
    names(fit$base) <- names(comparison)
    for (k in seq_along(comparison)) {
      predK <- pred[pred[[dep]] %in% comparison[[k]], ]
      predK[[dep]] <- ifelse(predK[[dep]] == comparison[[k]][1], 0, 1)
      fit$base[[k]] <- fit_qap_model(mod          = mod,
                                     pred         = predK,
                                     family       = family,
                                     estimator    = estimator,
                                     use_fixest   = use_fixest,
                                     fixest_se_cluster = fixest_se_cluster,
                                     use_robust_errors = use_robust_errors,
                                     main_vars    = main,
                                     has_random   = has_random,
                                     reference    = reference)
    }
  }

  # --- GPU fast path (gaussian, no random, no fixest, no comparison) ---
  if (use_gpu && family == "gaussian" && !has_random && !use_fixest &&
      is.null(comparison) && !large) {

    if (nullhyp == "qapy") {
      gpu_res <- gpu_batch_ols_css(data         = data,
                                   parsed       = parsed,
                                   mode         = mode,
                                   diag         = diag,
                                   groups       = groups,
                                   reps         = reps,
                                   baseline_fit = fit$base,
                                   perm_var     = NULL)
      fit$lower  <- gpu_res$lower
      fit$larger <- gpu_res$larger
      fit$abs    <- gpu_res$abs

    } else if (nullhyp == "qapspp") {
      n_coefs <- length(fit$base$coefficients)
      fit$lower  <- matrix(NA, nrow = 2, ncol = n_coefs)
      fit$larger <- fit$abs <- fit$lower
      colnames(fit$lower) <- colnames(fit$larger) <-
        colnames(fit$abs)  <- names(fit$base$coefficients)

      for (xi in main) {
        test_val <- data[[xi]]
        if (!is.numeric(test_val)) {
          warning("Cannot residualise non-numeric predictor '", xi,
                  "'. Skipping qapspp for this variable.")
          next
        }
        xR <- residualise_predictor(xi, pred, main,
                                    has_random   = has_random,
                                    rand_formula = rand_part)
        data_resid <- data
        data_resid[[xi]] <- residuals_to_array(xR, data[[xi]], valid, pred,
                                                large, valid_list)

        gpu_res <- gpu_batch_ols_css(data         = data_resid,
                                     parsed       = parsed,
                                     mode         = mode,
                                     diag         = diag,
                                     groups       = groups,
                                     reps         = reps,
                                     baseline_fit = fit$base,
                                     perm_var     = xi)
        fit$lower[, xi]  <- gpu_res$lower[, xi]
        fit$larger[, xi] <- gpu_res$larger[, xi]
        fit$abs[, xi]    <- gpu_res$abs[, xi]
      }
    }

  } else {

  # --- set up future plan ---
  old_plan <- setup_future_plan(ncores)
  on.exit({
    future::plan(old_plan)
    options(future.globals.maxSize = attr(old_plan, "old_maxSize"))
  }, add = TRUE)

  # --- progress bar ---
  total_reps <- if (nullhyp == "qapy") reps else reps * length(main)
  has_progressr <- requireNamespace("progressr", quietly = TRUE)

  .run_cpu_perms <- function(p) {
  # --- permutation testing ---
  if (nullhyp == "qapy") {
    res <- run_permutations(
      reps, QAPcssPermEst,
      data.     = data,
      perm_var. = NULL,
      mode.     = mode,
      diag.     = diag,
      mod.      = mod,
      groups.   = groups,
      fit.      = if (is.null(comparison)) fit$base else fit$base,
      family.   = family,
      estimator. = estimator,
      use_fixest. = use_fixest,
      fixest_se_cluster. = fixest_se_cluster,
      use_robust_errors. = use_robust_errors,
      has_random. = has_random,
      main_vars. = main,
      data_vars. = data_vars,
      parsed.   = parsed,
      comp.     = comparison,
      reference. = reference,
      p         = p
    )

    if (is.null(comparison)) {
      agg <- aggregate_perm_results(res, reps)
      fit$lower  <<- agg$lower
      fit$larger <<- agg$larger
      fit$abs    <<- agg$abs
    } else {
      res_valid <- Filter(Negate(is.null), res)
      n_valid   <- length(res_valid)
      fit$lower <<- fit$larger <<- fit$abs <<-
        vector("list", length(comparison))
      names(fit$lower) <<- names(fit$larger) <<-
        names(fit$abs) <<- names(comparison)
      resL <- unlist(unlist(res_valid, recursive = FALSE), recursive = FALSE)
      for (k in seq_along(comparison)) {
        cn <- names(comparison)[k]
        fit$lower[[k]]  <<- Reduce("+", resL[names(resL) == paste0(cn, ".lower")], 0) / n_valid
        fit$larger[[k]] <<- Reduce("+", resL[names(resL) == paste0(cn, ".larger")], 0) / n_valid
        fit$abs[[k]]    <<- Reduce("+", resL[names(resL) == paste0(cn, ".abs")], 0) / n_valid
      }
    }

  } else if (nullhyp == "qapspp") {

    # Initialise p-value matrices
    if (is.null(comparison)) {
      if (family != "multinom") {
        n_coefs <- length(fit$base$coefficients)
        fit$lower  <<- matrix(NA, nrow = 2, ncol = n_coefs)
        fit$larger <<- fit$abs <<- fit$lower
        colnames(fit$lower) <<- colnames(fit$larger) <<-
          colnames(fit$abs)  <<- names(fit$base$coefficients)
      } else {
        ncat <- if (large) {
          length(na.omit(unique(as.vector(unlist(data[[dep]])))))
        } else {
          length(na.omit(unique(as.vector(data[[dep]]))))
        }
        n_coefs <- length(fit$base$coefficients)
        fit$lower  <<- matrix(NA, nrow = 2 * (ncat - 1), ncol = n_coefs)
        fit$larger <<- fit$abs <<- fit$lower
        colnames(fit$lower) <<- colnames(fit$larger) <<-
          colnames(fit$abs) <<- names(fit$base$coefficients)
      }
    } else {
      fit$lower <<- fit$larger <<- fit$abs <<-
        vector("list", length(comparison))
      names(fit$lower) <<- names(fit$larger) <<-
        names(fit$abs) <<- names(comparison)
      for (k in seq_along(comparison)) {
        n_coefs <- length(fit$base[[k]]$coefficients)
        fit$lower[[k]] <<- matrix(NA, nrow = 2, ncol = n_coefs)
        fit$larger[[k]] <<- fit$abs[[k]] <<- fit$lower[[k]]
        colnames(fit$lower[[k]]) <<- colnames(fit$larger[[k]]) <<-
          colnames(fit$abs[[k]]) <<- names(fit$base[[k]]$coefficients)
      }
    }

    for (xi in main) {
      # Check if predictor is character/factor (can't residualise)
      test_val <- if (!large) data[[xi]] else data[[xi]][[1]]
      if (!is.numeric(test_val)) {
        warning("Cannot residualise non-numeric predictor '", xi,
                "'. Skipping qapspp for this variable.")
        next
      }

      xR <- residualise_predictor(xi, pred, main,
                                  has_random   = has_random,
                                  rand_formula = rand_part)

      data_resid <- data
      data_resid[[xi]] <- residuals_to_array(xR, data[[xi]], valid, pred,
                                              large, valid_list)

      res <- run_permutations(
        reps, QAPcssPermEst,
        data.     = data_resid,
        perm_var. = xi,
        mode.     = mode,
        diag.     = diag,
        mod.      = mod,
        groups.   = groups,
        fit.      = if (is.null(comparison)) fit$base else fit$base,
        family.   = family,
        estimator. = estimator,
        use_fixest. = use_fixest,
        fixest_se_cluster. = fixest_se_cluster,
        use_robust_errors. = use_robust_errors,
        has_random. = has_random,
        main_vars. = main,
        data_vars. = data_vars,
        parsed.   = parsed,
        comp.     = comparison,
        reference. = reference,
        p         = p
      )

      if (is.null(comparison)) {
        agg <- aggregate_perm_results(res, reps)
        fit$lower[, xi]  <<- agg$lower
        fit$larger[, xi] <<- agg$larger
        fit$abs[, xi]    <<- agg$abs
      } else {
        res_valid <- Filter(Negate(is.null), res)
        n_valid   <- length(res_valid)
        resL <- unlist(unlist(res_valid, recursive = FALSE), recursive = FALSE)
        for (k in seq_along(comparison)) {
          cn <- names(comparison)[k]
          fit$lower[[k]][, xi]  <<- Reduce("+", resL[names(resL) == paste0(cn, ".lower")], 0) / n_valid
          fit$larger[[k]][, xi] <<- Reduce("+", resL[names(resL) == paste0(cn, ".larger")], 0) / n_valid
          fit$abs[[k]][, xi]    <<- Reduce("+", resL[names(resL) == paste0(cn, ".abs")], 0) / n_valid
        }
      }
    }
  }
  }  # end .run_cpu_perms

  if (has_progressr) {
    progressr::with_progress({
      p <- progressr::progressor(steps = total_reps)
      .run_cpu_perms(p)
    })
  } else {
    .run_cpu_perms(NULL)
  }
  } # end else (CPU path)

  # --- confusion matrix for binomial ---
  if (family == "binomial" && is.null(comparison)) {
    bm <- fit$base$base_model
    if (!inherits(bm, "gmm")) {
      predicted <- fitted(bm)
      actual    <- pred[[dep]]
      fit$confusion_matrix <- probabilistic_confusion_matrix(
        actual = actual, predicted_prob = predicted,
        n_draws = 1000, seed = seed
      )
    }
  }

  # --- package results ---
  fit$nullhyp   <- nullhyp
  fit$family    <- family
  fit$groups    <- unique(unlist(groups))
  fit$diag      <- diag
  fit$mode      <- mode
  fit$reps      <- reps
  fit$reference <- reference
  fit$comp      <- comparison
  fit$random    <- c(sender    = ris,
                     receiver  = rir,
                     perceiver = rip,
                     nets      = rin)
  fit$robust_se <- use_robust_errors
  fit$estimator <- estimator

  # Propagate family-specific parameters
  if (is.null(comparison) && !is.null(fit$base$theta))
    fit$theta <- fit$base$theta
  if (is.null(comparison) && !is.null(fit$base$zi_coefficients))
    fit$zi_coefficients <- fit$base$zi_coefficients

  if (family == "multinom")
    names(fit)[names(fit) == "t"] <- "z"

  class(fit) <- "QAPCSS"
  return(fit)
}
