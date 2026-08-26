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


#' Decide whether robust standard errors can actually be used
#'
#' Whether HC3 standard errors are available is fixed by the family,
#' estimator and random-effect structure, so it is resolved once up front
#' rather than inside \code{fit_qap_model()}, which runs for every
#' permutation and would emit one warning per rep.
#'
#' @param use_robust_errors Logical; what the user asked for.
#' @param family Character; model family.
#' @param estimator Character; \code{"standard"} or \code{"gmm"}.
#' @param has_random Logical; does the model have random effects?
#'
#' @return Logical; the effective value of \code{use_robust_errors}.
#' @keywords internal

resolve_robust_errors <- function(use_robust_errors, family, estimator,
                                  has_random) {
  if (!isTRUE(use_robust_errors)) return(FALSE)

  if (estimator == "gmm") {
    # gmm::gmm() is called with vcov = "MDS", so its standard errors are
    # already heteroskedasticity-robust: the flag is redundant rather than
    # unavailable.  It is switched off so the OLS-form HC3(), which is not a
    # valid sandwich for a nonlinear GMM fit, is not applied on top.
    warning("use_robust_errors is redundant for estimator = \"gmm\": ",
            "gmm::gmm() already reports heteroskedasticity-robust ",
            "(vcov = \"MDS\") standard errors.")
    return(FALSE)
  }

  if (has_random) {
    warning("Robust (HC3) standard errors are not defined for mixed models; ",
            "model-based standard errors are reported instead.")
    return(FALSE)
  }

  if (family == "zip") {
    # vcovHC() cannot compute leverages for a zeroinfl fit, so the
    # uncorrected (HC0) sandwich is used instead of HC3.
    message("family = \"zip\": robust standard errors use the uncorrected ",
            "(HC0) sandwich; the HC3 small-sample correction is not ",
            "available for zero-inflated fits.")
  }

  # Only the gaussian, no-random case is served by the closed-form HC3();
  # everything else needs the GLM sandwich.
  if (family != "gaussian" && !requireNamespace("sandwich", quietly = TRUE)) {
    warning("use_robust_errors = TRUE for family = \"", family,
            "\" requires the 'sandwich' package ",
            "(install.packages('sandwich')); ",
            "model-based standard errors are reported instead.")
    return(FALSE)
  }

  TRUE
}


#' Align a named vector to a set of coefficient names
#'
#' Returns a vector with one entry per name in \code{coef_names}, filled from
#' \code{v} where the names match and \code{NA} otherwise.  This keeps
#' \code{fit$t} the same length and order as \code{fit$coefficients} even when
#' \code{summary()} silently drops rows for aliased (perfectly collinear)
#' terms -- otherwise \code{compare_perm_to_baseline()} would recycle two
#' vectors of different lengths and produce meaningless p-values.
#'
#' @param v Named numeric vector.
#' @param coef_names Character vector of coefficient names.
#'
#' @return Named numeric vector of length \code{length(coef_names)}.
#' @keywords internal

align_to_coefs <- function(v, coef_names) {
  out <- rep(NA_real_, length(coef_names))
  names(out) <- coef_names
  if (!is.null(names(v))) {
    common <- intersect(names(v), coef_names)
    out[common] <- v[common]
  } else if (length(v) == length(coef_names)) {
    out[] <- v
  }
  out
}


#' Heteroskedasticity-consistent (HC3) standard errors for a fitted model
#'
#' For a linear model the closed-form \code{HC3()} is exact.  For any other
#' model class the leverages have to come from the IRLS-weighted design
#' matrix, so \pkg{sandwich} is used; applying \code{HC3()} to a GLM's
#' deviance residuals and unweighted \code{X} returns a number that is not an
#' HC3 standard error at all.  Model classes \pkg{sandwich} cannot handle
#' (\pkg{lme4}, \pkg{glmmTMB}) fall back to model-based standard errors with a
#' warning.
#'
#' @param base_model A fitted model object.
#' @param X Numeric matrix of predictors, excluding the intercept (lm only).
#' @param resid Numeric vector of residuals (lm only).
#'
#' @return Named numeric vector of standard errors, or \code{NULL} when robust
#'   standard errors are unavailable for this model class, in which case the
#'   caller falls back to model-based standard errors.
#' @keywords internal

robust_se <- function(base_model, X = NULL, resid = NULL) {
  if (inherits(base_model, "lm") && !inherits(base_model, "glm")) {
    se <- HC3(X, resid)
    names(se) <- c("(Intercept)", colnames(X))
    return(se)
  }

  if (inherits(base_model, c("merMod", "lmerMod", "glmerMod", "glmmTMB"))) {
    warning("Robust (HC3) standard errors are not defined for mixed models; ",
            "model-based standard errors are reported instead.")
    return(NULL)
  }

  if (!requireNamespace("sandwich", quietly = TRUE)) {
    warning("use_robust_errors = TRUE for a '", class(base_model)[1],
            "' model requires the 'sandwich' package ",
            "(install.packages('sandwich')); ",
            "model-based standard errors are reported instead.")
    return(NULL)
  }

  vc <- tryCatch(sandwich::vcovHC(base_model, type = "HC3"),
                 error = function(e) NULL)

  if (is.null(vc)) {
    # Some classes (pscl::zeroinfl among them) provide estfun() and bread()
    # but no hatvalues(), so vcovHC() cannot form the HC3 leverage
    # adjustment.  The uncorrected (HC0) sandwich is still available and is
    # still heteroskedasticity-consistent, just without the small-sample
    # correction.
    vc <- tryCatch(sandwich::sandwich(base_model), error = function(e) NULL)
  }

  if (is.null(vc)) {
    warning("Robust standard errors are unavailable for a '",
            class(base_model)[1],
            "' model; model-based standard errors are reported instead.")
    return(NULL)
  }
  sqrt(diag(vc))
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
#' @param baseline Logical; is this the baseline fit rather than a
#'   permutation?  Goodness-of-fit extras (R-squared, dispersion,
#'   zero-inflation coefficients) are only ever read off the baseline, so
#'   they are skipped for permutation fits, of which there are
#'   \code{reps * n_predictors}.
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
                          reference        = NULL,
                          baseline         = FALSE) {
  fit <- list()
  dep_var <- all.vars(mod)[1]
  nx  <- length(main_vars)

  # --- multinomial ---
  if (family == "multinom") {
    pred[[dep_var]] <- as.factor(pred[[dep_var]])
    if (!is.null(reference)) {
      pred[[dep_var]] <- relevel(pred[[dep_var]], ref = reference)
    }
    if (!requireNamespace("nnet", quietly = TRUE))
      stop("Package 'nnet' is required for multinomial models.")
    base_model       <- nnet::multinom(mod, data = pred, trace = FALSE)
    fit$coefficients <- coefficients(base_model)
    fit$t            <- coefficients(base_model) /
                          summary(base_model)$standard.errors
    fit$base_model   <- base_model
    return(fit)
  }

  # --- GMM ---
  if (estimator == "gmm") {
    if (!requireNamespace("gmm", quietly = TRUE))
      stop("Package 'gmm' is required for GMM estimation.")
    y_vec <- pred[[dep_var]]
    x_mat <- cbind(1, as.matrix(pred[, main_vars, drop = FALSE]))

    if (!family %in% c("binomial", "poisson", "negbin", "zip"))
      stop("GMM estimator is available for binomial, poisson, negbin, ",
           "and zip families.")

    # Deterministic, data-driven starting values. See gmm_start(): the
    # previous rnorm() start was redrawn on every permutation and overflowed
    # the exp() in the count moment conditions routinely.
    gmm_args <- list(
      x  = list(y = y_vec, x = x_mat),
      t0 = gmm_start(y_vec, x_mat, family),
      wmatrix = "optimal", vcov = "MDS",
      optfct = "nlminb",
      control = list(eval.max = 10000)
    )

    has_extra_param <- family %in% c("negbin", "zip")

    gmm_args$g <- switch(family,
                         binomial = logit_moments,
                         poisson  = poisson_moments,
                         negbin   = negbin_moments,   # coefs + log(alpha)
                         zip      = zip_moments)      # coefs + logit(pi)

    base_model <- do.call(gmm::gmm, gmm_args)

    # gmm() is already called with vcov = "MDS", which yields
    # heteroskedasticity-robust standard errors; resolve_robust_errors()
    # has already turned use_robust_errors off for this estimator.

    # Extract coefficients and t-values BEFORE any subsetting,
    # so summary.gmm sees consistent dimensions.
    all_coefs <- base_model$coefficients
    all_t     <- summary(base_model)$coefficients[, 3]

    # For negbin/zip, strip the extra nuisance parameter
    if (has_extra_param) {
      fit$coefficients <- all_coefs[1:(nx + 1)]
    } else {
      fit$coefficients <- all_coefs
    }
    names(fit$coefficients) <- c("(Intercept)", main_vars)

    fit$t <- if (has_extra_param) all_t[1:(nx + 1)] else all_t
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
      if (use_robust_errors) robust_se(base_model)  # warns; not available
      fit$t <- align_to_coefs(summary(base_model)$coefficients$cond[, 3],
                              names(fit$coefficients))
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

      se <- if (use_robust_errors) robust_se(base_model) else NULL
      if (!is.null(se)) {
        # vcovHC() on a zeroinfl fit covers count and zero components; the
        # count-component entries are prefixed "count_".
        names(se) <- sub("^count_", "", names(se))
        fit$t <- fit$coefficients /
                   align_to_coefs(se, names(fit$coefficients))
      } else {
        fit$t <- align_to_coefs(summary(base_model)$coefficients$count[, 3],
                                names(fit$coefficients))
      }
      fit$zi_coefficients <- base_model$coefficients$zero
    }
    fit$base_model <- base_model
    return(fit)
  }

  # --- standard estimation ---
  if (!has_random) {
    if (use_fixest) {
      if (!requireNamespace("fixest", quietly = TRUE))
        stop("Package 'fixest' is required when fixest_se_cluster or fixed effects with | ")
      fe_family <- if (family == "negbin") "negbin" else family
      base_model <- fixest::feglm(mod, data = pred,
                                  family = fe_family,
                                  cluster = fixest_se_cluster)
      fit$coefficients <- c("(Intercept)" = NA, base_model$coefficients)

      # fixest computes its own heteroskedasticity-robust vcov on the demeaned
      # design.  Recomputing HC3 from pred[, main_vars] would use leverages
      # from a design matrix that excludes the absorbed fixed effects.
      fe_vcov <- if (use_robust_errors) {
        vcov(base_model, vcov = "hetero")
      } else {
        vcov(base_model)
      }
      fe_se <- sqrt(diag(fe_vcov))
      fit$t <- align_to_coefs(base_model$coefficients / fe_se,
                              names(fit$coefficients))

      if (baseline && family == "gaussian") {
        r2s <- tryCatch(fixest::r2(base_model), error = function(e) NULL)
        if (!is.null(r2s)) {
          fit$r.squared     <- r2s[["r2"]]
          fit$adj.r.squared <- r2s[["ar2"]]
        }
      }
    } else {
      if (family == "gaussian") {
        base_model <- lm(mod, data = pred)
      } else if (family == "negbin") {
        if (!requireNamespace("MASS", quietly = TRUE))
          stop("Package 'MASS' is required for negative binomial models.")
        base_model <- MASS::glm.nb(mod, data = pred)
        if (baseline) fit$theta <- base_model$theta
      } else {
        base_model <- glm(mod, data = pred, family = family)
      }
      fit$coefficients <- base_model$coefficients

      se <- if (use_robust_errors) {
        robust_se(base_model,
                  X     = as.matrix(pred[, main_vars, drop = FALSE]),
                  resid = residuals(base_model, type = "response"))
      } else {
        NULL
      }

      # summary() is not cheap and was previously called up to three times
      # per fit; build it at most once and only when it is actually needed.
      sm <- if (is.null(se) || (baseline && family == "gaussian")) {
        summary(base_model)
      } else {
        NULL
      }

      if (baseline && family == "gaussian") {
        fit$r.squared     <- sm$r.squared
        fit$adj.r.squared <- sm$adj.r.squared
      }

      if (!is.null(se)) {
        fit$t <- fit$coefficients /
                   align_to_coefs(se, names(fit$coefficients))
      } else {
        fit$t <- align_to_coefs(sm$coefficients[, 3],
                                names(fit$coefficients))
      }
    }
  } else {
    # --- random effects ---
    if (family == "gaussian") {
      if (!requireNamespace("lme4", quietly = TRUE))
        stop("Package 'lme4' is required for  random effects.")
      base_model <- lme4::lmer(mod, data = pred)
    } else if (family == "negbin") {
      # glmmTMB for mixed negative binomial
      if (!requireNamespace("glmmTMB", quietly = TRUE))
        stop("Package 'glmmTMB' is required for mixed negative binomial models.")
      base_model <- glmmTMB::glmmTMB(mod, data = pred,
                                      family = glmmTMB::nbinom2())
      fit$coefficients <- glmmTMB::fixef(base_model)$cond
      if (use_robust_errors) robust_se(base_model)  # warns; not available
      fit$t <- align_to_coefs(summary(base_model)$coefficients$cond[, 3],
                              names(fit$coefficients))
      fit$theta <- glmmTMB::sigma(base_model)
      fit$random.intercepts <- list()
      re <- glmmTMB::ranef(base_model)$cond
      for (rV in names(re)) {
        fit$random.intercepts[[rV]] <- re[[rV]][, 1]
      }
      fit$base_model <- base_model
      return(fit)
    } else {
      if (!requireNamespace("lme4", quietly = TRUE))
        stop("Package 'lme4' is required for  random effects.")
      base_model <- lme4::glmer(mod, data = pred, family = family,
                                control = lme4::glmerControl(
                                  calc.derivs = FALSE,
                                  optimizer   = "bobyqa"),
                                nAGQ = 0)
      fit$log_lik <- logLik(base_model)
    }

    fit$coefficients <- summary(base_model)$coefficients[, 1]
    if (use_robust_errors) robust_se(base_model)  # warns; not available
    fit$t <- align_to_coefs(summary(base_model)$coefficients[, 3],
                            names(fit$coefficients))

    fit$random.intercepts <- list()
    for (rV in names(coefficients(base_model))) {
      fit$random.intercepts[[rV]] <- coefficients(base_model)[[rV]][, 1]
    }
  }

  fit$base_model <- base_model
  return(fit)
}


#' Warn once about coefficients that could not be estimated
#'
#' Called on the baseline fit only.  Permutation fits reuse the same design,
#' so warning there would emit one message per permutation.
#'
#' @param fit A fit list from \code{fit_qap_model()}, or a list of them.
#'
#' @return Invisibly NULL.
#' @keywords internal

warn_aliased_coefs <- function(fit) {
  coefs <- fit$coefficients
  if (is.null(coefs) || !is.numeric(coefs)) return(invisible(NULL))
  bad <- names(coefs)[is.na(coefs)]
  # A fixest fit deliberately carries an NA intercept (it is absorbed).
  bad <- setdiff(bad, "(Intercept)")
  if (length(bad))
    warning("Coefficient(s) ", paste(bad, collapse = ", "),
            " could not be estimated (perfect collinearity among ",
            "predictors); their permutation p-values will be NA.")
  invisible(NULL)
}


#' Row labels for the permutation p-value matrices
#'
#' \code{lower}, \code{larger} and \code{abs} have one row per test
#' statistic: the coefficient itself, and its t-value.  Defined in one place
#' so every path -- CPU, GPU, qapy and qapspp -- labels them identically.
#'
#' @return Character vector of length 2.
#' @keywords internal

qap_stat_rows <- function() c("b", "t")


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

  # rbind() would otherwise name the rows after this function's own
  # arguments ("perm_coefs", "perm_t"), which leaks into the user-facing
  # result. Row 1 compares the coefficient, row 2 the t-value.
  #
  # Multinomial fits stack a whole BLOCK of each: their coefficients are a
  # (ncat - 1) x k matrix, so rbind() gives (ncat - 1) b-rows above
  # (ncat - 1) t-rows. Label both blocks so this path and qap_init_pmats()
  # produce identically-named matrices.
  h <- nrow(pres) %/% 2L
  rownames(pres) <- rownames(bres) <- if (h == 1L) {
    qap_stat_rows()
  } else {
    cats <- rownames(base_fit$coefficients)
    c(paste0("b:", cats), paste0("t:", cats))
  }

  out <- list()
  if (is.null(xi)) {
    out$lower  <- pres <= bres
    out$larger <- pres >= bres
    out$abs    <- abs(pres) >= abs(bres)
    # The permuted statistics themselves, retained so confint() can form
    # percentile intervals from the same draws the p-values come from.
    out$draw   <- pres
  } else {
    out$lower  <- (pres <= bres)[, xi]
    out$larger <- (pres >= bres)[, xi]
    out$abs    <- (abs(pres) >= abs(bres))[, xi]
    out$draw   <- pres[, xi]
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
  draws <- resL[names(resL) == "draw"]

  list(
    lower  = Reduce("+", resL[names(resL) == "lower"],  0) / n_valid,
    larger = Reduce("+", resL[names(resL) == "larger"], 0) / n_valid,
    abs    = Reduce("+", resL[names(resL) == "abs"],    0) / n_valid,
    # One row per successful permutation. For a whole-model permutation each
    # element is a 2 x k matrix (b and t); when a single predictor is being
    # tested it is a length-2 vector.
    draws  = qap_stack_draws(draws)
  )
}


#' Stack per-permutation draws into b and t matrices
#'
#' @param draws List of per-permutation results: 2 x k matrices, or length-2
#'   vectors when a single predictor was tested.
#'
#' @return A list with \code{b} and \code{t}, each \code{n_perm x k}, or NULL
#'   if nothing was retained.
#' @keywords internal

qap_stack_draws <- function(draws) {
  if (!length(draws)) return(NULL)
  first <- draws[[1]]

  if (is.matrix(first)) {
    if (nrow(first) != 2L) return(NULL)   # multinomial: not a b/t pair
    nm <- colnames(first)
    b  <- do.call(rbind, lapply(draws, function(d) d[1, ]))
    tv <- do.call(rbind, lapply(draws, function(d) d[2, ]))
    colnames(b) <- colnames(tv) <- nm
  } else {
    if (length(first) != 2L) return(NULL)
    b  <- matrix(vapply(draws, function(d) d[[1]], numeric(1)), ncol = 1)
    tv <- matrix(vapply(draws, function(d) d[[2]], numeric(1)), ncol = 1)
  }
  list(b = b, t = tv)
}


#' Validate QAP input data
#'
#' Checks that the data list contains the expected variables as matrices or
#' lists of matrices.
#'
#' @param data Named list of matrices (QAP) or arrays (QAP).
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

  what  <- if (css) "a 3-dimensional array" else "a matrix"
  ndim  <- if (css) 3L else 2L
  ok_shape <- function(v) {
    if (css) length(dim(v)) == 3L else is.matrix(v)
  }

  if (!large) {
    if (!ok_shape(y)) stop("data[['", dep, "']] must be ", what, ".")
  } else {
    for (i in seq_along(y)) {
      if (!ok_shape(y[[i]]))
        stop("data[['", dep, "']][[", i, "]] must be ", what, ".")
    }
  }

  # Predictors must have exactly the same shape as the dependent variable.
  # Without this check a mismatched predictor makes the internal
  # valid[is.na(x)] <- FALSE recycle a wrong-length logical index and produce
  # silently incorrect data rather than an error.
  preds <- setdiff(parsed$all_data_vars, c(dep, structural_vars))
  preds <- intersect(preds, names(data))

  for (v in preds) {
    xv <- data[[v]]
    if (large) {
      if (!is.list(xv))
        stop("Predictor '", v, "' must be a list of ", ndim,
             "-dimensional objects, matching data[['", dep, "']].")
      if (length(xv) != length(y))
        stop("Predictor '", v, "' has ", length(xv), " network(s) but '",
             dep, "' has ", length(y), ".")
      for (i in seq_along(xv)) {
        if (!ok_shape(xv[[i]]))
          stop("data[['", v, "']][[", i, "]] must be ", what, ".")
        if (!identical(dim(xv[[i]]), dim(y[[i]])))
          stop("data[['", v, "']][[", i, "]] has dimensions ",
               paste(dim(xv[[i]]), collapse = " x "), " but data[['", dep,
               "']][[", i, "]] has ", paste(dim(y[[i]]), collapse = " x "),
               ".")
      }
    } else {
      if (is.list(xv))
        stop("Predictor '", v, "' is a list but '", dep,
             "' is a single network.")
      if (!ok_shape(xv))
        stop("data[['", v, "']] must be ", what, ".")
      if (!identical(dim(xv), dim(y)))
        stop("data[['", v, "']] has dimensions ",
             paste(dim(xv), collapse = " x "), " but data[['", dep,
             "']] has ", paste(dim(y), collapse = " x "), ".")
    }
  }

  invisible(TRUE)
}


#' Set up parallel processing with the future framework
#'
#' Returns the function that undoes the change, rather than leaving the
#' caller to remember which pieces of global state were touched.  Use it as
#' \code{restore <- setup_future_plan(ncores); on.exit(restore(), add = TRUE)}.
#'
#' @param ncores Number of cores. If NULL or 1, uses sequential processing.
#'
#' @return A zero-argument function that restores the previous plan and the
#'   previous \code{future.globals.maxSize}.
#' @keywords internal

setup_future_plan <- function(ncores = NULL) {
  old_plan    <- future::plan()
  old_maxSize <- getOption("future.globals.maxSize")

  # Raise the globals size limit so large network datasets can be shipped to
  # workers without hitting the default cap.
  options(future.globals.maxSize = +Inf)
  if (is.null(ncores) || ncores <= 1) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = ncores)
  }

  function() {
    future::plan(old_plan)
    options(future.globals.maxSize = old_maxSize)
    invisible(NULL)
  }
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
    if (!requireNamespace("lme4", quietly = TRUE))
      stop("Package 'lme4' is required for  random effects.")
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
#' @param mode Character; \code{"digraph"} or \code{"graph"}.  In
#'   \code{"graph"} mode only the upper triangle is vectorised, so the
#'   residuals are mirrored onto the lower triangle to keep the matrix
#'   symmetric under permutation.
#'
#' @return A matrix (or list of matrices) with residuals.
#' @keywords internal

residuals_to_matrix <- function(xR, original_matrix, pred, large = FALSE,
                                mode = "digraph") {
  # In graph mode make_qap_data() keeps only the upper triangle, so writing
  # the residuals back would leave the lower triangle holding the original
  # (unresidualised) values.  RMPerm() permutes rows and columns together and
  # moves cells across the diagonal, so the matrix must be symmetrised.
  mirror <- function(m) {
    if (mode != "graph") return(m)
    lt <- lower.tri(m)
    m[lt] <- t(m)[lt]
    m
  }

  out <- original_matrix
  if (!large) {
    out[pred$location] <- xR
    out <- mirror(out)
  } else {
    for (net_id in unique(as.character(pred$nv))) {
      idx <- as.character(pred$nv) == net_id
      net_num <- as.integer(net_id)
      out[[net_num]][pred$location[idx]] <- xR[idx]
      out[[net_num]] <- mirror(out[[net_num]])
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
#' @param mode Character; \code{"directed"} or \code{"undirected"}.  In
#'   \code{"undirected"} mode only the upper triangle of each perceiver slice
#'   is vectorised, so the residuals are mirrored onto the lower triangle.
#'
#' @return An array (or list of arrays) with residuals.
#' @keywords internal

residuals_to_array <- function(xR, original_array, valid, pred,
                               large = FALSE, valid_list = NULL,
                               mode = "directed") {
  # As in residuals_to_matrix(): make_css_data() vectorises only the upper
  # triangle of each slice in undirected mode, while RMPerm(CSS = TRUE)
  # permutes cells across the diagonal.  Mirror so the slices stay symmetric.
  mirror <- function(a) {
    if (mode != "undirected") return(a)
    for (i in seq_len(dim(a)[3])) {
      sl <- a[, , i]
      lt <- lower.tri(sl)
      sl[lt] <- t(sl)[lt]
      a[, , i] <- sl
    }
    a
  }

  if (!large) {
    out <- array(NA, dim = dim(original_array))
    out[valid] <- xR
    return(mirror(out))
  } else {
    out <- original_array
    for (gr in seq_along(original_array)) {
      out[[gr]] <- array(NA, dim = dim(original_array[[gr]]))
      out[[gr]][valid_list[[gr]]] <- xR[pred$nv == gr]
      out[[gr]] <- mirror(out[[gr]])
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


#' Permute one network, or a list of networks with per-network groups
#'
#' When \code{groups} is a list it holds one grouping vector per network, so
#' each network must be permuted with its own vector.  Passing the whole list
#' to \code{RMPerm()} would coerce it with \code{as.character()} and produce a
#' meaningless grouping.
#'
#' @param m Matrix/array, or list thereof.
#' @param groups Grouping vector; an unnamed list of grouping vectors (one
#'   per network); a named list of per-dimension groupings (multi-mode); or
#'   NULL.
#' @param CSS Logical; TRUE for 3D CSS arrays.
#' @param multi_mode Logical; each dimension its own node set?
#'
#' @return The permuted object, matching the structure of \code{m}.
#' @keywords internal

perm_networks <- function(m, groups = NULL, CSS = FALSE,
                          multi_mode = FALSE) {
  if (!is.list(m))
    return(RMPerm(m, groups, CSS = CSS, multi_mode = multi_mode))

  # A NAMED list is one grouping per dimension and applies to every network;
  # an UNNAMED list is one grouping per network.
  per_network <- is.list(groups) && is.null(names(groups))

  if (per_network) {
    if (length(groups) != length(m))
      stop("groups is a list of length ", length(groups),
           " but there are ", length(m), " networks.")
    return(Map(function(mi, gi) RMPerm(mi, gi, CSS = CSS,
                                       multi_mode = multi_mode), m, groups))
  }

  lapply(m, RMPerm, groups = groups, CSS = CSS, multi_mode = multi_mode)
}

