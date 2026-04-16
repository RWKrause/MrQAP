#' Parse a QAP formula into its components
#'
#' Extracts dependent variable, main effects, fixed effects, and random effects
#' from a formula. Distinguishes between fixest-style FE (\code{|} without
#' parentheses) and lme4-style random effects (\code{(1|var)}).
#'
#' @param formula A formula object.
#' @param fixest_se_cluster Optional cluster variable name for fixest.
#'
#' @return A list with components: \code{dependent}, \code{main},
#'   \code{fixed_effects}, \code{has_random}, \code{use_fixest},
#'   \code{all_data_vars}.
#' @keywords internal

parse_qap_formula <- function(formula, fixest_se_cluster = NULL) {
  dependent <- all.vars(formula)[1]

  formula_str <- paste(deparse(formula, width.cutoff = 500), collapse = " ")
  has_pipe  <- grepl("\\|", formula_str)
  has_paren <- grepl("\\(", formula_str)

  if (has_pipe && has_paren) {
    # lme4-style random effects: y ~ x1 + x2 + (1|sv) + (1|rv)
    main <- all.vars(reformulas::nobars(formula))[-1]
    fixed_effects <- NULL
    use_fixest    <- FALSE
    has_random    <- TRUE
    all_data_vars <- main
  } else if (has_pipe && !has_paren) {
    # fixest-style fixed effects: y ~ x1 + x2 | fe1 + fe2
    main          <- all.vars(formula[[3]][[2]])
    fixed_effects <- all.vars(formula[[3]][[3]])
    has_random    <- FALSE
    use_fixest    <- TRUE
    all_data_vars <- c(main, fixed_effects)
  } else {
    main          <- all.vars(formula[-1])
    fixed_effects <- NULL
    has_random    <- FALSE
    use_fixest    <- !is.null(fixest_se_cluster)
    all_data_vars <- main
  }

  if (!is.null(fixest_se_cluster)) {
    use_fixest <- TRUE
    if (!(fixest_se_cluster %in% all_data_vars)) {
      all_data_vars <- c(all_data_vars, fixest_se_cluster)
    }
  }

  list(dependent     = dependent,
       main          = main,
       fixed_effects = fixed_effects,
       has_random    = has_random,
       use_fixest    = use_fixest,
       all_data_vars = all_data_vars)
}


#' Build internal formula with additional random intercepts
#'
#' Appends random-intercept terms for sender, receiver, perceiver, network,
#' and custom grouping to a user-supplied formula.
#'
#' @param formula A formula object.
#' @param rin Logical; random intercept for networks.
#' @param ris Logical; random intercept for senders.
#' @param rir Logical; random intercept for receivers.
#' @param rip Logical; random intercept for perceivers (CSS only).
#'
#' @return An updated formula.
#' @keywords internal

build_internal_formula <- function(formula,
                                   rin = FALSE,
                                   ris = FALSE,
                                   rir = FALSE,
                                   rip = FALSE) {
  mod_str <- paste(deparse(formula, width.cutoff = 500), collapse = " ")
  rand_int <- ""
  if (rin) rand_int <- paste(rand_int, "+ (1|nv)")
  if (ris) rand_int <- paste(rand_int, "+ (1|sv)")
  if (rir) rand_int <- paste(rand_int, "+ (1|rv)")
  if (rip) rand_int <- paste(rand_int, "+ (1|pv)")
  if (nchar(trimws(rand_int)) > 0) {
    as.formula(paste(mod_str, rand_int))
  } else {
    formula
  }
}


#' Fit a single model (shared by baseline and permutation fitting)
#'
#' Dispatches to lm, glm, MASS::glm.nb, pscl::zeroinfl, fixest::feglm,
#' lme4::lmer/glmer, glmmTMB::glmmTMB, gmm::gmm, or nnet::multinom
#' depending on arguments.
#'
#' @param mod Formula for the model.
#' @param pred Data frame of vectorised network data.
#' @param family Model family.
#' @param estimator \code{"standard"} or \code{"gmm"}.
#' @param use_fixest Logical; use fixest::feglm.
#' @param fixest_se_cluster Cluster variable name for fixest.
#' @param use_robust_errors Logical; use HC3 correction.
#' @param main_vars Character vector of main predictor names.
#' @param has_random Logical; any random effects?
#' @param reference Reference category for multinomial.
#'
#' @return A list with \code{coefficients}, \code{t}, and \code{base_model}.
#' @keywords internal

fit_qap_model <- function(mod, pred, family,
                          estimator        = "standard",
                          use_fixest       = FALSE,
                          fixest_se_cluster = NULL,
                          use_robust_errors = FALSE,
                          main_vars        = NULL,
                          has_random       = FALSE,
                          reference        = NULL) {
  fit <- list()
  dep_var <- all.vars(mod)[1]
  nx  <- length(main_vars)

  # --- multinomial ---
  if (family == "multinom") {
    pred[[dep_var]] <- as.factor(pred[[dep_var]])
    if (!is.null(reference)) {
      pred[[dep_var]] <- relevel(pred[[dep_var]], ref = reference)
    }
    base_model       <- nnet::multinom(mod, data = pred, trace = FALSE)
    fit$coefficients <- coefficients(base_model)
    fit$t            <- coefficients(base_model) /
                          summary(base_model)$standard.errors
    fit$base_model   <- base_model
    return(fit)
  }

  # --- GMM ---
  if (estimator == "gmm") {
    y_vec <- pred[[dep_var]]
    x_mat <- cbind(1, as.matrix(pred[, main_vars, drop = FALSE]))

    gmm_args <- list(
      x = list(y = y_vec, x = x_mat),
      t0 = rnorm(nx + 1),
      wmatrix = "optimal", vcov = "MDS",
      optfct = "nlminb",
      control = list(eval.max = 10000)
    )

    has_extra_param <- FALSE

    if (family == "binomial") {
      gmm_args$g <- logit_moments
      base_model <- do.call(gmm::gmm, gmm_args)
      resid <- logit_resid(base_model)
    } else if (family == "poisson") {
      gmm_args$g <- poisson_moments
      base_model <- do.call(gmm::gmm, gmm_args)
      resid <- poisson_resid(base_model)
    } else if (family == "negbin") {
      # Negative binomial GMM: estimate regression + log(alpha) jointly
      gmm_args$g  <- negbin_moments
      gmm_args$t0 <- rnorm(nx + 2)   # extra param for log(alpha)
      base_model   <- do.call(gmm::gmm, gmm_args)
      resid <- negbin_resid(base_model)
      has_extra_param <- TRUE
    } else if (family == "zip") {
      # Zero-inflated Poisson GMM: regression + logit(pi) for zero-inflation
      gmm_args$g  <- zip_moments
      gmm_args$t0 <- rnorm(nx + 2)   # extra param for logit(pi)
      base_model   <- do.call(gmm::gmm, gmm_args)
      resid <- zip_resid(base_model)
      has_extra_param <- TRUE
    } else {
      stop("GMM estimator is available for binomial, poisson, negbin, ",
           "and zip families.")
    }

    # Extract coefficients and t-values BEFORE any subsetting,
    # so summary.gmm sees consistent dimensions.
    all_coefs <- base_model$coefficients
    if (!use_robust_errors) {
      all_t <- summary(base_model)$coefficients[, 3]
    }

    # For negbin/zip, strip the extra nuisance parameter
    if (has_extra_param) {
      fit$coefficients <- all_coefs[1:(nx + 1)]
    } else {
      fit$coefficients <- all_coefs
    }
    names(fit$coefficients) <- c("(Intercept)", main_vars)

    if (use_robust_errors) {
      xv <- as.matrix(pred[, main_vars, drop = FALSE])
      fit$t <- fit$coefficients / HC3(xv, resid)
    } else {
      if (has_extra_param) {
        fit$t <- all_t[1:(nx + 1)]
      } else {
        fit$t <- all_t
      }
    }
    names(fit$t) <- names(fit$coefficients)
    fit$base_model <- base_model
    fit$estimator  <- "gmm"
    return(fit)
  }

  # --- zero-inflated Poisson (standard, no random, no fixest) ---
  if (family == "zip" && estimator == "standard") {
    if (has_random) {
      # glmmTMB for mixed ZIP
      if (!requireNamespace("glmmTMB", quietly = TRUE))
        stop("Package 'glmmTMB' is required for mixed ZIP models.")
      base_model <- glmmTMB::glmmTMB(mod, data = pred,
                                      family = poisson(),
                                      ziformula = ~1)
      fit$coefficients <- glmmTMB::fixef(base_model)$cond
      resid <- residuals(base_model, type = "response")
      fit$t <- summary(base_model)$coefficients$cond[, 3]
      names(fit$t) <- names(fit$coefficients)
      fit$zi_coefficients <- glmmTMB::fixef(base_model)$zi
      fit$random.intercepts <- list()
      re <- glmmTMB::ranef(base_model)$cond
      for (rV in names(re)) {
        fit$random.intercepts[[rV]] <- re[[rV]][, 1]
      }
    } else {
      if (!requireNamespace("pscl", quietly = TRUE))
        stop("Package 'pscl' is required for zero-inflated Poisson models.")
      base_model <- pscl::zeroinfl(mod, data = pred, dist = "poisson")
      fit$coefficients <- base_model$coefficients$count
      resid <- residuals(base_model, type = "response")
      if (use_robust_errors) {
        xv <- as.matrix(pred[, main_vars, drop = FALSE])
        fit$t <- fit$coefficients / HC3(xv, resid)
      } else {
        fit$t <- summary(base_model)$coefficients$count[, 3]
      }
      fit$zi_coefficients <- base_model$coefficients$zero
    }
    fit$base_model <- base_model
    return(fit)
  }

  # --- standard estimation ---
  if (!has_random) {
    if (use_fixest) {
      fe_family <- if (family == "negbin") "negbin" else family
      base_model <- fixest::feglm(mod, data = pred,
                                  family = fe_family,
                                  cluster = fixest_se_cluster)
      fit$coefficients <- c("(Intercept)" = NA, base_model$coefficients)
      resid <- residuals(base_model)

      if (use_robust_errors) {
        xv   <- as.matrix(pred[, main_vars, drop = FALSE])
        hc   <- HC3(xv, resid)
        fit$t <- fit$coefficients / c(NA, hc[-1])
      } else {
        fe_se <- sqrt(diag(vcov(base_model)))
        fit$t <- c("(Intercept)" = NA,
                    base_model$coefficients / fe_se)
      }
      names(fit$t) <- names(fit$coefficients)

      if (family == "gaussian") {
        r2s <- tryCatch(fixest::r2(base_model), error = function(e) NULL)
        if (!is.null(r2s)) {
          fit$r.squared     <- r2s[["r2"]]
          fit$adj.r.squared <- r2s[["ar2"]]
        }
      }
    } else {
      if (family == "gaussian") {
        base_model        <- lm(mod, data = pred)
        fit$r.squared     <- summary(base_model)$r.squared
        fit$adj.r.squared <- summary(base_model)$adj.r.squared
      } else if (family == "negbin") {
        if (!requireNamespace("MASS", quietly = TRUE))
          stop("Package 'MASS' is required for negative binomial models.")
        base_model <- MASS::glm.nb(mod, data = pred)
        fit$theta  <- base_model$theta
      } else {
        base_model <- glm(mod, data = pred, family = family)
      }
      fit$coefficients <- base_model$coefficients
      resid <- residuals(base_model)

      if (use_robust_errors) {
        xv <- as.matrix(pred[, main_vars, drop = FALSE])
        fit$t <- fit$coefficients / HC3(xv, resid)
      } else {
        fit$t <- summary(base_model)$coefficients[, 3]
      }
    }
  } else {
    # --- random effects ---
    if (family == "gaussian") {
      base_model <- lme4::lmer(mod, data = pred)
    } else if (family == "negbin") {
      # glmmTMB for mixed negative binomial
      if (!requireNamespace("glmmTMB", quietly = TRUE))
        stop("Package 'glmmTMB' is required for mixed negative binomial models.")
      base_model <- glmmTMB::glmmTMB(mod, data = pred,
                                      family = glmmTMB::nbinom2())
      fit$coefficients <- glmmTMB::fixef(base_model)$cond
      resid <- residuals(base_model, type = "response")
      if (use_robust_errors) {
        xv <- as.matrix(pred[, main_vars, drop = FALSE])
        fit$t <- fit$coefficients / HC3(xv, resid)
      } else {
        fit$t <- summary(base_model)$coefficients$cond[, 3]
        names(fit$t) <- names(fit$coefficients)
      }
      fit$theta <- glmmTMB::sigma(base_model)
      fit$random.intercepts <- list()
      re <- glmmTMB::ranef(base_model)$cond
      for (rV in names(re)) {
        fit$random.intercepts[[rV]] <- re[[rV]][, 1]
      }
      fit$base_model <- base_model
      return(fit)
    } else {
      base_model <- lme4::glmer(mod, data = pred, family = family,
                                control = lme4::glmerControl(
                                  calc.derivs = FALSE,
                                  optimizer   = "bobyqa"),
                                nAGQ = 0)
      fit$log_lik <- logLik(base_model)
    }

    fit$coefficients <- summary(base_model)$coefficients[, 1]
    resid <- residuals(base_model)

    if (use_robust_errors) {
      xv <- as.matrix(pred[, main_vars, drop = FALSE])
      fit$t <- fit$coefficients / HC3(xv, resid)
    } else {
      fit$t <- summary(base_model)$coefficients[, 3]
    }

    fit$random.intercepts <- list()
    for (rV in names(coefficients(base_model))) {
      fit$random.intercepts[[rV]] <- coefficients(base_model)[[rV]][, 1]
    }
  }

  fit$base_model <- base_model
  return(fit)
}


#' Compare permutation results to the baseline fit
#'
#' @param perm_coefs Named numeric vector of permutation coefficients.
#' @param perm_t Named numeric vector of permutation t-values.
#' @param base_fit List with \code{coefficients} and \code{t} from baseline.
#' @param xi Optional; name of the variable being tested (for qapspp).
#'
#' @return A list with logical matrices \code{lower}, \code{larger}, \code{abs}.
#' @keywords internal

compare_perm_to_baseline <- function(perm_coefs, perm_t, base_fit,
                                     xi = NULL) {
  pres <- rbind(perm_coefs, perm_t)
  bres <- rbind(base_fit$coefficients, base_fit$t)

  out <- list()
  if (is.null(xi)) {
    out$lower  <- pres <= bres
    out$larger <- pres >= bres
    out$abs    <- abs(pres) >= abs(bres)
  } else {
    out$lower  <- (pres <= bres)[, xi]
    out$larger <- (pres >= bres)[, xi]
    out$abs    <- (abs(pres) >= abs(bres))[, xi]
  }
  return(out)
}


#' Aggregate permutation results into proportions
#'
#' @param results List of permutation comparison results.
#' @param reps Number of permutations.
#'
#' @return A list with proportion matrices \code{lower}, \code{larger},
#'   \code{abs}.
#' @keywords internal

aggregate_perm_results <- function(results, reps) {
  # Drop NULLs (failed permutations)
  results <- Filter(Negate(is.null), results)
  n_valid <- length(results)
  if (n_valid == 0) stop("All permutations failed to converge.")
  if (n_valid < reps) {
    warning(reps - n_valid, " of ", reps,
            " permutations failed and were excluded.")
  }
  resL <- unlist(results, recursive = FALSE)
  list(
    lower  = Reduce("+", resL[names(resL) == "lower"],  0) / n_valid,
    larger = Reduce("+", resL[names(resL) == "larger"], 0) / n_valid,
    abs    = Reduce("+", resL[names(resL) == "abs"],    0) / n_valid
  )
}


#' Validate QAP input data
#'
#' Checks that the data list contains the expected variables as matrices or
#' lists of matrices.
#'
#' @param data Named list of matrices (QAPglm) or arrays (QAPcss).
#' @param parsed Output of \code{parse_qap_formula()}.
#' @param css Logical; TRUE for CSS (3D array) data.
#'
#' @return Invisible TRUE if valid; stops with an error otherwise.
#' @keywords internal

validate_qap_input <- function(data, parsed, css = FALSE) {
  dep <- parsed$dependent
  if (!(dep %in% names(data))) {
    stop("Dependent variable '", dep, "' not found in data.")
  }
  # sv, rv, nv, pv are auto-generated structural columns, not user data
  structural_vars <- c("sv", "rv", "nv", "pv")
  for (v in parsed$all_data_vars) {
    if (v %in% structural_vars) next
    if (!(v %in% names(data))) {
      stop("Predictor '", v, "' not found in data.")
    }
  }

  y <- data[[dep]]
  large <- is.list(y)

  if (!css) {
    if (!large) {
      if (!is.matrix(y)) stop("data[['", dep, "']] must be a matrix.")
    } else {
      for (i in seq_along(y)) {
        if (!is.matrix(y[[i]]))
          stop("data[['", dep, "']][[", i, "]] must be a matrix.")
      }
    }
  } else {
    if (!large) {
      if (length(dim(y)) != 3)
        stop("data[['", dep, "']] must be a 3-dimensional array.")
    } else {
      for (i in seq_along(y)) {
        if (length(dim(y[[i]])) != 3)
          stop("data[['", dep, "']][[", i, "]] must be a 3D array.")
      }
    }
  }
  invisible(TRUE)
}


#' Set up parallel processing with the future framework
#'
#' @param ncores Number of cores. If NULL or 1, uses sequential processing.
#'
#' @return Invisibly returns the previous plan.
#' @keywords internal

setup_future_plan <- function(ncores = NULL) {
  old_plan <- future::plan()
  # Raise the globals size limit so large network datasets can be
  # shipped to workers without hitting the default 500 MB cap.
  old_maxSize <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = +Inf)
  if (is.null(ncores) || ncores <= 1) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = ncores)
  }
  # Store old value so callers can restore it if desired
  attr(old_plan, "old_maxSize") <- old_maxSize
  invisible(old_plan)
}


#' Run permutations in parallel using future
#'
#' @param reps Number of permutations.
#' @param FUN Function to run per permutation.
#' @param ... Additional arguments passed to FUN.
#' @param p A \code{progressr::progressor} object or NULL.
#'
#' @return List of permutation results.
#' @keywords internal

run_permutations <- function(reps, FUN, ..., p = NULL) {
  if (!is.null(p)) {
    results <- future.apply::future_lapply(
      1:reps,
      function(i, ...) {
        res <- FUN(i, ...)
        p()
        res
      },
      ...,
      future.seed = TRUE
    )
  } else {
    results <- future.apply::future_lapply(
      1:reps,
      FUN,
      ...,
      future.seed = TRUE
    )
  }
  return(results)
}


#' Residualise a predictor on all others (for qapspp)
#'
#' Fits a linear model of \code{xi} on all other predictors in the data frame
#' and returns the residuals.
#'
#' @param xi Name of the variable to residualise.
#' @param pred Data frame with all predictors.
#' @param main_vars Character vector of main predictor names.
#' @param has_random Logical; use lme4?
#' @param rand_formula Character; random-effect part of formula (e.g.,
#'   \code{"+ (1|sv)"}).
#'
#' @return Numeric vector of residuals.
#' @keywords internal

residualise_predictor <- function(xi, pred, main_vars,
                                  has_random = FALSE,
                                  rand_formula = "") {
  others <- setdiff(main_vars, xi)
  modx_str <- paste(xi, "~ 1")
  if (length(others) > 0) {
    modx_str <- paste(modx_str, paste(others, collapse = " + "), sep = " + ")
  }
  if (has_random && nchar(trimws(rand_formula)) > 0) {
    modx_str <- paste(modx_str, rand_formula)
  }
  modx <- as.formula(modx_str)

  if (!has_random) {
    xm <- lm(modx, data = pred)
  } else {
    xm <- lme4::lmer(modx, data = pred)
  }
  residuals(xm)
}


#' Put residuals back into matrix form
#'
#' After residualising on the data frame, maps residuals back to the
#' original matrix structure using location indices.
#'
#' @param xR Numeric vector of residuals.
#' @param original_matrix The original matrix (or list of matrices).
#' @param pred Data frame with \code{location} and optionally \code{nv}.
#' @param large Logical; TRUE if y is a list of matrices.
#'
#' @return A matrix (or list of matrices) with residuals.
#' @keywords internal

residuals_to_matrix <- function(xR, original_matrix, pred, large = FALSE) {
  out <- original_matrix
  if (!large) {
    out[pred$location] <- xR
  } else {
    for (net_id in unique(as.character(pred$nv))) {
      idx <- as.character(pred$nv) == net_id
      net_num <- as.integer(net_id)
      out[[net_num]][pred$location[idx]] <- xR[idx]
    }
  }
  return(out)
}


#' Put residuals back into 3D array form (CSS)
#'
#' @param xR Numeric vector of residuals.
#' @param original_array The original 3D array (or list of arrays).
#' @param valid Logical array of valid positions.
#' @param pred Data frame with \code{nv}.
#' @param large Logical; TRUE if y is a list of arrays.
#' @param valid_list List of valid arrays (when large = TRUE).
#'
#' @return An array (or list of arrays) with residuals.
#' @keywords internal

residuals_to_array <- function(xR, original_array, valid, pred,
                               large = FALSE, valid_list = NULL) {
  if (!large) {
    n <- dim(original_array)[1]
    out <- array(NA, dim = c(n, n, n))
    out[valid] <- xR
    return(out)
  } else {
    out <- original_array
    for (gr in seq_along(original_array)) {
      out[[gr]] <- array(NA, dim = dim(original_array[[gr]]))
      out[[gr]][valid_list[[gr]]] <- xR[pred$nv == gr]
    }
    return(out)
  }
}


#' Internal auxiliary function to fit the base model
#'
#' Thin wrapper around \code{fit_qap_model()} retained for backward
#' compatibility.  New code should call \code{fit_qap_model()} directly.
#'
#' @param mod Formula.
#' @param rand Logical; any random intercepts?
#' @param family Character; model family.
#' @param pred Data frame.
#' @param nx Integer; number of predictors.
#' @param y Dependent data (unused, kept for compatibility).
#' @param use_robust_errors Logical.
#' @param estimator Character; "standard" or "gmm".
#' @param use_fixest Logical.
#' @param fixest_se_cluster Character or NULL.
#' @param main_vars Character vector.
#' @param reference Character or NULL.
#'
#' @return A list from \code{fit_qap_model()}.
#' @keywords internal

fit_base <- function(mod,
                     rand,
                     family,
                     pred,
                     nx = NULL,
                     y  = NULL,
                     use_robust_errors = FALSE,
                     estimator    = "standard",
                     use_fixest   = FALSE,
                     fixest_se_cluster = NULL,
                     main_vars    = NULL,
                     reference    = NULL) {

  # Derive main_vars if not provided (backward compat)
  if (is.null(main_vars)) {
    all_v <- all.vars(mod)
    dep   <- all_v[1]
    main_vars <- setdiff(names(pred), c(dep, "location", "nv", "sv", "rv",
                                        "pv", "ov",
                                        paste0("ov", 1:20)))
    main_vars <- intersect(all_v[-1], main_vars)
  }

  fit_qap_model(mod          = mod,
                pred         = pred,
                family       = family,
                estimator    = estimator,
                use_fixest   = use_fixest,
                fixest_se_cluster = fixest_se_cluster,
                use_robust_errors = use_robust_errors,
                main_vars    = main_vars,
                has_random   = rand,
                reference    = reference)
}


#' Internal auxiliary function to fit a permutation model and compare
#'
#' Thin wrapper that fits one permuted model and compares coefficients to
#' the baseline.  New code should call \code{fit_qap_model()} +
#' \code{compare_perm_to_baseline()} directly.
#'
#' @param family. Character; model family.
#' @param predx Data frame of permuted data.
#' @param ref Character; reference category (multinom).
#' @param xi. Character; variable being tested (qapspp) or NULL.
#' @param mod Formula.
#' @param fitx Baseline fit list.
#' @param rand Logical; random intercepts?
#' @param nx. Integer; number of predictors.
#' @param use_robust_errors Logical.
#' @param estimator Character.
#' @param use_fixest Logical.
#' @param fixest_se_cluster Character or NULL.
#' @param main_vars Character vector.
#'
#' @return A list with lower/larger/abs.
#' @keywords internal

fit_perm <- function(family.,
                     predx,
                     ref       = NULL,
                     xi.       = NULL,
                     mod,
                     fitx,
                     rand      = FALSE,
                     nx.       = NULL,
                     use_robust_errors = FALSE,
                     estimator    = "standard",
                     use_fixest   = FALSE,
                     fixest_se_cluster = NULL,
                     main_vars    = NULL) {

  if (is.null(main_vars)) {
    all_v <- all.vars(mod)
    dep   <- all_v[1]
    main_vars <- setdiff(names(predx), c(dep, "location", "nv", "sv", "rv",
                                         "pv", "ov",
                                         paste0("ov", 1:20)))
    main_vars <- intersect(all_v[-1], main_vars)
  }

  perm_fit <- fit_qap_model(mod          = mod,
                            pred         = predx,
                            family       = family.,
                            estimator    = estimator,
                            use_fixest   = use_fixest,
                            fixest_se_cluster = fixest_se_cluster,
                            use_robust_errors = use_robust_errors,
                            main_vars    = main_vars,
                            has_random   = rand,
                            reference    = ref)

  compare_perm_to_baseline(perm_fit$coefficients, perm_fit$t,
                            fitx, xi = xi.)
}


#' Internal auxiliary function to estimate one QAPglm permutation
#'
#' Called by \code{future.apply::future_lapply()} from \code{QAPglm()}.
#'
#' @param i Integer; iteration index.
#' @param data. Named list of matrices (same as \code{data} in QAPglm).
#' @param perm_var. NULL for qapy, or the name of the variable to permute
#'   (for qapspp, the residualised variable).
#' @param mode. Internal mode string ("digraph"/"graph").
#' @param diag. Logical.
#' @param mod. Formula for the model.
#' @param groups. Permutation groups.
#' @param fit. Baseline fit (list or list of lists for comparisons).
#' @param family. Model family.
#' @param estimator. Estimator type.
#' @param use_fixest. Logical.
#' @param fixest_se_cluster. Cluster variable.
#' @param use_robust_errors. Logical.
#' @param has_random. Logical.
#' @param main_vars. Character vector of main predictor names.
#' @param data_vars. Character vector of all data variable names.
#' @param parsed. Output of parse_qap_formula.
#' @param comp. Comparison list or NULL.
#' @param reference. Reference category or NULL.
#'
#' @return A list with lower/larger/abs comparison results.
#' @keywords internal

QAPglmPermEst <- function(i,
                          data.,
                          perm_var.,
                          mode.,
                          diag.,
                          mod.,
                          groups.,
                          fit.,
                          family.,
                          estimator.,
                          use_fixest.,
                          fixest_se_cluster.,
                          use_robust_errors.,
                          has_random.,
                          main_vars.,
                          data_vars.,
                          parsed.,
                          comp.,
                          reference.) {

  dep   <- parsed.$dependent
  large <- is.list(data.[[dep]])

  # --- permute ---
  d <- data.
  if (is.null(perm_var.)) {
    # qapy: permute dependent
    if (!large) {
      d[[dep]] <- RMPerm(d[[dep]], groups.)
    } else {
      d[[dep]] <- lapply(d[[dep]], RMPerm, groups = groups.)
    }
  } else {
    # qapspp: permute the (already residualised) variable
    if (!large) {
      d[[perm_var.]] <- RMPerm(d[[perm_var.]], groups.)
    } else {
      d[[perm_var.]] <- lapply(d[[perm_var.]], RMPerm, groups = groups.)
    }
  }

  # --- build data frame ---
  if (!large) {
    pred <- make_qap_data(y    = d[[dep]],
                          x    = d[data_vars.],
                          g    = groups.,
                          diag = diag.,
                          mode = mode.,
                          net  = 1,
                          perm = FALSE,
                          xi   = NULL)
  } else {
    pred_list <- vector("list", length(d[[dep]]))
    for (net in seq_along(d[[dep]])) {
      x2 <- lapply(data_vars., function(v) d[[v]][[net]])
      names(x2) <- data_vars.
      g2 <- if (!is.null(groups.)) groups.[[net]] else NULL
      pred_list[[net]] <- make_qap_data(y    = d[[dep]][[net]],
                                        x    = x2,
                                        g    = g2,
                                        diag = diag.,
                                        mode = mode.,
                                        net  = net,
                                        perm = FALSE,
                                        xi   = NULL)
    }
    pred <- do.call(rbind, pred_list)
  }

  names(pred)[names(pred) == "yv"] <- dep

  xi_arg <- if (!is.null(perm_var.)) perm_var. else NULL

  # --- no comparisons ---
  if (is.null(comp.)) {
    perm_fit <- tryCatch(
      fit_qap_model(mod          = mod.,
                    pred         = pred,
                    family       = family.,
                    estimator    = estimator.,
                    use_fixest   = use_fixest.,
                    fixest_se_cluster = fixest_se_cluster.,
                    use_robust_errors = use_robust_errors.,
                    main_vars    = main_vars.,
                    has_random   = has_random.,
                    reference    = reference.),
      error = function(e) NULL
    )
    if (is.null(perm_fit)) return(NULL)

    return(compare_perm_to_baseline(perm_fit$coefficients, perm_fit$t,
                                    fit., xi = xi_arg))
  }

  # --- with comparisons ---
  xresL <- vector("list", length(comp.))
  names(xresL) <- names(comp.)

  for (k in seq_along(comp.)) {
    predK <- pred[pred[[dep]] %in% comp.[[k]], ]
    predK[[dep]] <- ifelse(predK[[dep]] == comp.[[k]][1], 0, 1)

    perm_fit <- tryCatch(
      fit_qap_model(mod          = mod.,
                    pred         = predK,
                    family       = family.,
                    estimator    = estimator.,
                    use_fixest   = use_fixest.,
                    fixest_se_cluster = fixest_se_cluster.,
                    use_robust_errors = use_robust_errors.,
                    main_vars    = main_vars.,
                    has_random   = has_random.,
                    reference    = reference.),
      error = function(e) NULL
    )
    if (is.null(perm_fit)) return(NULL)

    xresL[[k]] <- compare_perm_to_baseline(perm_fit$coefficients, perm_fit$t,
                                            fit.[[k]], xi = xi_arg)
  }

  return(xresL)
}


#' Internal auxiliary function to estimate one CSS permutation
#'
#' Called by \code{future.apply::future_lapply()} from \code{QAPcss()}.
#'
#' @param i Integer; iteration index.
#' @param data. Named list of 3D arrays.
#' @param perm_var. NULL for qapy, or variable name (qapspp).
#' @param mode. Character; "directed" or "undirected".
#' @param diag. Logical.
#' @param mod. Formula.
#' @param groups. Permutation groups.
#' @param fit. Baseline fit.
#' @param family. Model family.
#' @param estimator. Estimator type.
#' @param use_fixest. Logical.
#' @param fixest_se_cluster. Cluster variable.
#' @param use_robust_errors. Logical.
#' @param has_random. Logical.
#' @param main_vars. Character vector of main predictor names.
#' @param data_vars. Character vector of all data variable names.
#' @param parsed. Output of parse_qap_formula.
#' @param comp. Comparison list or NULL.
#' @param reference. Reference category or NULL.
#'
#' @return A list with lower/larger/abs comparison results.
#' @keywords internal

QAPcssPermEst <- function(i,
                          data.,
                          perm_var.,
                          mode.,
                          diag.,
                          mod.,
                          groups.,
                          fit.,
                          family.,
                          estimator.,
                          use_fixest.,
                          fixest_se_cluster.,
                          use_robust_errors.,
                          has_random.,
                          main_vars.,
                          data_vars.,
                          parsed.,
                          comp.,
                          reference.) {

  dep   <- parsed.$dependent
  large <- is.list(data.[[dep]])
  nx    <- length(main_vars.)

  # Track categories for sufficient-data check
  y_cat <- na.omit(unique(as.vector(unlist(data.[[dep]]))))

  sufficient_data <- FALSE
  trial <- 0
  max_trials <- 10000

  while (!sufficient_data && trial < max_trials) {
    trial <- trial + 1

    # --- permute ---
    d <- data.
    if (is.null(perm_var.)) {
      # qapy: permute dependent
      if (!large) {
        d[[dep]] <- RMPerm(d[[dep]], groups., CSS = TRUE)
      } else {
        d[[dep]] <- lapply(d[[dep]], RMPerm, groups = groups., CSS = TRUE)
      }
    } else {
      # qapspp: permute residualised variable
      if (!large) {
        d[[perm_var.]] <- RMPerm(d[[perm_var.]], groups., CSS = TRUE)
      } else {
        d[[perm_var.]] <- lapply(d[[perm_var.]], RMPerm,
                                  groups = groups., CSS = TRUE)
      }
    }

    # --- build data frame ---
    if (!large) {
      x_list <- lapply(data_vars., function(v) d[[v]])
      names(x_list) <- data_vars.
      pred <- make_css_data(y = d[[dep]], x = x_list,
                            nets = 1,
                            diag = diag., mode = mode.)$pred
    } else {
      pred_list <- vector("list", length(d[[dep]]))
      for (gr in seq_along(d[[dep]])) {
        xgr <- lapply(data_vars., function(v) d[[v]][[gr]])
        names(xgr) <- data_vars.
        pred_list[[gr]] <- make_css_data(y = d[[dep]][[gr]], x = xgr,
                                         nets = gr,
                                         diag = diag., mode = mode.)$pred
      }
      pred <- do.call(rbind, pred_list)
    }

    names(pred)[names(pred) == "yv"] <- dep

    # --- sufficient data check ---
    if (family. != "multinom" && is.null(comp.)) {
      y_ok <- length(na.omit(unique(pred[[dep]]))) > 1
    } else {
      y2_cat <- na.omit(unique(pred[[dep]]))
      y_present <- all(y_cat %in% y2_cat)
      y_mult <- all(table(pred[[dep]]) > 2)
      y_ok <- y_present && y_mult
    }

    x_ok <- TRUE
    num_preds <- pred[, data_vars.[data_vars. %in% names(pred)], drop = FALSE]
    num_preds <- num_preds[, sapply(num_preds, is.numeric), drop = FALSE]
    if (ncol(num_preds) > 0) {
      x_ok <- all(sapply(num_preds, function(col) length(unique(col)) > 1))
    }

    # Correlation check for comparisons
    if (nrow(pred) != 0 && x_ok && y_ok && !is.null(comp.)) {
      for (k in seq_along(comp.)) {
        pred2 <- pred[pred[[dep]] %in% comp.[[k]], ]
        pred2[[dep]] <- ifelse(pred2[[dep]] == comp.[[k]][1], 0, 1)
        check_cols <- c(dep, intersect(main_vars., names(pred2)))
        if (length(check_cols) > 1) {
          cors <- tryCatch(
            cor(pred2[, check_cols, drop = FALSE], use = "complete.obs"),
            error = function(e) NULL
          )
          if (is.null(cors) || any(is.na(cors))) {
            y_ok <- x_ok <- FALSE
          } else {
            diag(cors) <- 0
            if (any(abs(cors) > 0.9999)) y_ok <- x_ok <- FALSE
          }
        }
      }
    }

    sufficient_data <- y_ok && x_ok
  }

  if (trial >= max_trials) {
    stop("Cannot find valid permutation after ", max_trials, " trials.")
  }

  xi_arg <- if (!is.null(perm_var.)) perm_var. else NULL

  # --- no comparisons ---
  if (is.null(comp.)) {
    perm_fit <- tryCatch(
      fit_qap_model(mod          = mod.,
                    pred         = pred,
                    family       = family.,
                    estimator    = estimator.,
                    use_fixest   = use_fixest.,
                    fixest_se_cluster = fixest_se_cluster.,
                    use_robust_errors = use_robust_errors.,
                    main_vars    = main_vars.,
                    has_random   = has_random.,
                    reference    = reference.),
      error = function(e) NULL
    )
    if (is.null(perm_fit)) return(NULL)

    return(compare_perm_to_baseline(perm_fit$coefficients, perm_fit$t,
                                    fit., xi = xi_arg))
  }

  # --- with comparisons ---
  xresL <- vector("list", length(comp.))
  names(xresL) <- names(comp.)

  for (k in seq_along(comp.)) {
    predK <- pred[pred[[dep]] %in% comp.[[k]], ]
    predK[[dep]] <- ifelse(predK[[dep]] == comp.[[k]][1], 0, 1)

    perm_fit <- tryCatch(
      fit_qap_model(mod          = mod.,
                    pred         = predK,
                    family       = family.,
                    estimator    = estimator.,
                    use_fixest   = use_fixest.,
                    fixest_se_cluster = fixest_se_cluster.,
                    use_robust_errors = use_robust_errors.,
                    main_vars    = main_vars.,
                    has_random   = has_random.,
                    reference    = reference.),
      error = function(e) NULL
    )
    if (is.null(perm_fit)) return(NULL)

    xresL[[k]] <- compare_perm_to_baseline(perm_fit$coefficients, perm_fit$t,
                                            fit.[[k]], xi = xi_arg)
  }

  return(xresL)
}
