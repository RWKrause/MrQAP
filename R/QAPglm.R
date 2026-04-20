#' Parameter and p-value estimation with MRQAP
#'
#' Estimates regression models on network data (matrices) with permutation-based
#' inference.
#'
#' @param formula A formula describing the model.
#'   The left-hand side is the dependent variable name in \code{data};
#'   the right-hand side lists predictor variable names.
#'   Use \code{|} (without parentheses) for fixest fixed effects:
#'   \code{y ~ x1 + x2 | fe1}.
#'   Use \code{(1|var)} for lme4-style random effects:
#'   \code{y ~ x1 + (1|sv)}.
#'
#' @param data Named list; each element is a square \code{matrix} (or a
#'   \code{list} of matrices if there are multiple networks).  Names must
#'   correspond to the variables in \code{formula}.
#'
#' @param family Character; model family (default \code{"gaussian"}).
#'   Supported: \code{"gaussian"}, \code{"binomial"}, \code{"poisson"},
#'   \code{"negbin"} (negative binomial), \code{"zip"} (zero-inflated
#'   Poisson), and \code{"multinom"} (multinomial).
#'
#' @param mode Character; \code{"directed"} (default) or \code{"undirected"}.
#'
#' @param diag Logical; include diagonal values (default \code{FALSE}).
#'
#' @param nullhyp Character; \code{"qapspp"} (default, recommended) or
#'   \code{"qapy"}.
#'
#' @param estimator Character; \code{"standard"} (default) or \code{"gmm"}
#'   (for binomial, poisson, negbin, and zip families).
#'
#' @param reps Integer; number of permutations (default 1000).
#' @param seed Integer; optional random seed.
#' @param groups Vector (or list of vectors); permutation grouping.
#' @param ncores Integer; number of cores for parallel processing via the
#'   \code{future} framework.
#' @param fixest_se_cluster Character; cluster variable for fixest.
#'
#' @param comparison Named list of length-2 character vectors specifying
#'   category comparisons (e.g., \code{list(commission =
#'   c("false_positive","true_negative"))}).
#' @param reference Character; reference category (for multinomial or
#'   comparison models).
#'
#' @param random_intercept_nets,random_intercept_sender,random_intercept_receiver
#'   Logical; add random intercepts for networks, senders, receivers.
#'
#' @param use_robust_errors Logical; use HC3 error correction.
#' @param less_mem Logical; do not store the fitted model object.
#' @param use_gpu Logical; use GPU-accelerated batch OLS (gaussian only,
#'   requires \code{torch} package).
#'
#' @return An object of class \code{QAPRegression} (gaussian) or
#'   \code{QAPGLM} (other families).
#'
#' @export

QAPglm <- function(formula,
                   data,
                   family    = "gaussian",
                   mode      = "directed",
                   diag      = FALSE,
                   nullhyp   = "qapspp",
                   estimator = "standard",
                   reps      = 1000,
                   seed      = NULL,
                   groups    = NULL,
                   ncores    = NULL,
                   fixest_se_cluster = NULL,
                   comparison = NULL,
                   reference  = NULL,
                   random_intercept_nets     = FALSE,
                   random_intercept_sender   = FALSE,
                   random_intercept_receiver = FALSE,
                   use_robust_errors = FALSE,
                   less_mem   = FALSE,
                   use_gpu    = FALSE) {

  # --- seed ---
  if (!is.null(seed)) set.seed(seed)

  # --- parse formula ---
  parsed <- parse_qap_formula(formula, fixest_se_cluster)
  dep        <- parsed$dependent
  main       <- parsed$main
  # Filter out structural vars (sv, rv, nv, pv) that are auto-generated
  data_vars  <- intersect(parsed$all_data_vars, names(data))

  # --- validate ---
  validate_qap_input(data, parsed, css = FALSE)
  large <- is.list(data[[dep]])

  # --- mode encoding ---
  mode_internal <- if (mode == "directed") "digraph" else "graph"

  rin <- random_intercept_nets
  ris <- random_intercept_sender
  rir <- random_intercept_receiver

  # --- build internal formula ---
  mod <- build_internal_formula(formula,
                                rin = rin, ris = ris, rir = rir)
  mod_str <- paste(deparse(mod, width.cutoff = 500), collapse = " ")
  has_random <- grepl("\\(", mod_str) || parsed$has_random
  use_fixest <- parsed$use_fixest
  if (has_random && use_fixest) {
    warning("Cannot combine fixest FE and lme4 random effects. ",
            "Using lme4 random effects only.")
    use_fixest <- FALSE
  }

  mod <- as.formula(mod_str)

  # --- build data frame ---
  if (!large) {
    pred <- make_qap_data(y    = data[[dep]],
                          x    = data[data_vars],
                          g    = groups,
                          diag = diag,
                          mode = mode_internal,
                          net  = 1,
                          perm = FALSE,
                          xi   = NULL)
  } else {
    pred_list <- vector("list", length(data[[dep]]))
    for (net in seq_along(data[[dep]])) {
      x2 <- lapply(data_vars, function(v) data[[v]][[net]])
      names(x2) <- data_vars
      g2 <- if (!is.null(groups)) groups[[net]] else NULL
      pred_list[[net]] <- make_qap_data(y    = data[[dep]][[net]],
                                        x    = x2,
                                        g    = g2,
                                        diag = diag,
                                        mode = mode_internal,
                                        net  = net,
                                        perm = FALSE,
                                        xi   = NULL)
    }
    pred <- do.call(rbind, pred_list)
  }

  # rename yv -> dependent name
  names(pred)[names(pred) == "yv"] <- dep

  # --- handle comparisons ---
  if (!is.null(comparison) && is.null(reference)) {
    reference <- NULL
  }

  # --- baseline fit ---
  fit <- list()

  # Build random-intercept string for residualisation later
  rand_part <- ""
  if (rin) rand_part <- paste(rand_part, "+ (1|nv)")
  if (ris) rand_part <- paste(rand_part, "+ (1|sv)")
  if (rir) rand_part <- paste(rand_part, "+ (1|rv)")

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

  # --- single-predictor override ---
  if ((nullhyp == "qapspp") && (length(main) == 1)) nullhyp <- "qapy"

  # --- GPU fast path (gaussian, no random, no fixest, no comparison) ---
  if (use_gpu && family == "gaussian" && !has_random && !use_fixest &&
      is.null(comparison) && !large) {

    if (nullhyp == "qapy") {
      gpu_res <- gpu_batch_ols(data         = data,
                               parsed       = parsed,
                               mode         = mode_internal,
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
        xR <- residualise_predictor(xi, pred, main,
                                    has_random   = has_random,
                                    rand_formula = rand_part)
        data_resid <- data
        data_resid[[xi]] <- residuals_to_matrix(xR, data[[xi]], pred, large)

        gpu_res <- gpu_batch_ols(data         = data_resid,
                                 parsed       = parsed,
                                 mode         = mode_internal,
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
    if (nullhyp == "qapy") {
      res <- run_permutations(
        reps, QAPglmPermEst,
        data.     = data,
        perm_var. = NULL,
        mode.     = mode_internal,
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
        names(fit$lower)  <<- names(comparison)
        names(fit$larger) <<- names(comparison)
        names(fit$abs)    <<- names(comparison)
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
        n_coefs <- length(fit$base$coefficients)
        fit$lower  <<- matrix(NA, nrow = 2, ncol = n_coefs)
        fit$larger <<- fit$abs <<- fit$lower
        colnames(fit$lower) <<- colnames(fit$larger) <<-
          colnames(fit$abs)  <<- names(fit$base$coefficients)
      } else {
        fit$lower <<- fit$larger <<- fit$abs <<-
          vector("list", length(comparison))
        names(fit$lower) <<- names(fit$larger) <<-
          names(fit$abs)  <<- names(comparison)
        for (k in seq_along(comparison)) {
          n_coefs <- length(fit$base[[k]]$coefficients)
          fit$lower[[k]] <<- matrix(NA, nrow = 2, ncol = n_coefs)
          fit$larger[[k]] <<- fit$abs[[k]] <<- fit$lower[[k]]
          colnames(fit$lower[[k]]) <<- colnames(fit$larger[[k]]) <<-
            colnames(fit$abs[[k]])  <<- names(fit$base[[k]]$coefficients)
        }
      }

      for (xi in main) {
        # Residualise xi on other predictors
        xR <- residualise_predictor(xi, pred, main,
                                    has_random   = has_random,
                                    rand_formula = rand_part)

        # Put residuals back into matrix form
        data_resid <- data
        data_resid[[xi]] <- residuals_to_matrix(xR, data[[xi]], pred, large)

        res <- run_permutations(
          reps, QAPglmPermEst,
          data.     = data_resid,
          perm_var. = xi,
          mode.     = mode_internal,
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
  }

  # --- package results ---
  if (is.null(comparison)) {
    fit$coefficients <- fit$base$coefficients
    fit$t            <- fit$base$t
    if (!is.null(fit$base$r.squared)) {
      fit$r.squared     <- fit$base$r.squared
      fit$adj.r.squared <- fit$base$adj.r.squared
    }
    if (!is.null(fit$base$random.intercepts))
      fit$random.intercepts <- fit$base$random.intercepts
    if (!is.null(fit$base$theta))
      fit$theta <- fit$base$theta
    if (!is.null(fit$base$zi_coefficients))
      fit$zi_coefficients <- fit$base$zi_coefficients
    if (!less_mem) fit$simple_fit <- fit$base$base_model
  } else {
    if (!less_mem) {
      fit$simple_fits <- lapply(fit$base, `[[`, "base_model")
    }
  }

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

  fit$nullhyp   <- nullhyp
  fit$diag      <- diag
  fit$family    <- family
  fit$mode      <- mode
  fit$reps      <- reps
  fit$groups    <- unique(unlist(groups))
  fit$robust_se <- use_robust_errors
  fit$estimator <- estimator
  fit$comp      <- comparison
  fit$reference <- reference

  if (family == "gaussian" && is.null(comparison)) {
    class(fit) <- "QAPRegression"
  } else {
    class(fit) <- "QAPGLM"
  }
  return(fit)
}
