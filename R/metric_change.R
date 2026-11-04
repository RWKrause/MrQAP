# Session state: the LPM bias warning is given once per session, not once
# per fit.
.qap_state <- new.env(parent = emptyenv())


#' Seed for the probabilistic confusion-matrix draws
#'
#' The draws are always seeded, so a fit's classification metrics -- and the
#' differences in \code{metric_change}, which compare draws across models --
#' come out the same on every run.  The seed given to \code{QAP()} is used
#' when there is one.
#'
#' @param seed The \code{seed} argument of \code{QAP()}, or NULL.
#'
#' @return An integer seed.
#' @keywords internal

qap_confusion_seed <- function(seed) if (is.null(seed)) 350415L else seed


#' Does a vector take exactly two distinct values?
#'
#' Any type: a comparison model's raw outcome is often character.  Cheaper
#' than \code{length(unique(y)) == 2} on a long vector, since it never builds
#' the hash table.
#'
#' @param y A vector.
#'
#' @return Logical.
#' @keywords internal

qap_is_binary <- function(y) {
  y <- y[!is.na(y)]
  if (!length(y)) return(FALSE)
  other <- y != y[1]
  if (!any(other)) return(FALSE)
  all(!other | y == y[which.max(other)])
}


#' Observed outcome and predicted probabilities of a binary fit
#'
#' \describe{
#'   \item{glm}{Read off the model rather than the data frame it was fitted
#'     to, because \code{fixest::feglm()} drops observations (fixed-effect
#'     groups with no variation in the outcome) and the two would no longer
#'     line up.}
#'   \item{lpm}{A linear probability model.  Outcome and fitted values are
#'     mapped by \code{(v - lo) / (hi - lo)} with \code{lo}, \code{hi} the two
#'     observed outcome values, which is exactly the fit the 0/1 recoded
#'     outcome would have given.  Fitted values outside [0, 1] are counted,
#'     then truncated.}
#'   \item{gmm}{\code{gmm::gmm()} exposes no fitted values, so the logit
#'     probabilities are rebuilt from the coefficients on the same design
#'     \code{fit_qap_model()} hands it.}
#' }
#'
#' @param f A fit list from \code{fit_qap_model()}.
#' @param kind One of \code{"glm"}, \code{"lpm"}, \code{"gmm"}.
#' @param pred,dep,main_vars The data the model was fitted to (\code{gmm}).
#' @param levels The two outcome values, lower first (\code{lpm}).
#'
#' @return A list with \code{actual} (0/1), \code{prob}, and for an LPM
#'   \code{n_out}.
#' @keywords internal

qap_binary_outcome <- function(f, kind = "glm", pred = NULL, dep = NULL,
                               main_vars = NULL, levels = c(0, 1)) {
  bm <- f$base_model
  if (kind == "gmm") {
    x   <- cbind(1, as.matrix(pred[, main_vars, drop = FALSE]))
    eta <- drop(x %*% f$coefficients)
    if (!is.null(pred[[".qap_offset"]])) eta <- eta + pred[[".qap_offset"]]
    return(list(actual = pred[[dep]], prob = stats::plogis(eta)))
  }

  fv <- as.numeric(stats::fitted(bm))
  y  <- fv + as.numeric(stats::residuals(bm, type = "response"))
  if (kind == "glm") return(list(actual = round(y), prob = fv))

  span  <- levels[2] - levels[1]
  p     <- (fv - levels[1]) / span
  n_out <- sum(p < 0 | p > 1)
  list(actual = round((y - levels[1]) / span),
       prob   = pmin(pmax(p, 0), 1),
       n_out  = n_out)
}


#' Confusion matrix of a binary baseline fit
#'
#' @inheritParams qap_binary_outcome
#' @param seed Integer seed for the draws.
#'
#' @return A \code{QAPConfusionMatrix}, or NULL when there is no model to
#'   read.  An LPM's carries \code{n_out}, the number of fitted values that
#'   were outside [0, 1].
#' @keywords internal

qap_fit_confusion <- function(f, seed, kind = "glm", pred = NULL, dep = NULL,
                              main_vars = NULL, levels = c(0, 1)) {
  if (is.null(f) || (kind != "gmm" && is.null(f$base_model))) return(NULL)
  o  <- qap_binary_outcome(f, kind, pred, dep, main_vars, levels)
  cm <- probabilistic_confusion_matrix(actual = o$actual,
                                       predicted_prob = o$prob,
                                       n_draws = 1000, seed = seed)
  if (kind == "lpm") cm$n_out <- o$n_out
  cm
}


#' The eight classification metrics of one confusion matrix
#'
#' @param cm A \code{QAPConfusionMatrix}.
#'
#' @return Named numeric vector: each metric at the 0.5 cut-off
#'   (\code{_cut}) and averaged over the probabilistic draws (\code{_draw}).
#' @keywords internal

qap_class_metrics <- function(cm) {
  c(accuracy_cut     = cm$threshold_accuracy,
    accuracy_draw    = cm$accuracy,
    sensitivity_cut  = cm$threshold_sensitivity,
    sensitivity_draw = cm$sensitivity,
    specificity_cut  = cm$threshold_specificity,
    specificity_draw = cm$specificity,
    precision_cut    = cm$threshold_precision,
    precision_draw   = cm$precision)
}


#' Fixed-effect terms of a model formula, and the formula without one
#'
#' Random-effect bars are left alone; so is the fixed-effect part after a
#' fixest \code{|}, which \code{update()} cannot see into -- it reads
#' \code{x1 + x2 | fe} as a single term.
#'
#' @param mod The internal model formula.
#' @param use_fixest Logical; may the formula carry a fixest \code{|}?
#' @param term A term label to drop.
#'
#' @return \code{qap_fixed_terms()}: character vector of term labels.
#'   \code{qap_drop_term()}: a formula.  \code{qap_fixed_rhs()}: the
#'   fixed-effect right-hand side as a call.
#' @keywords internal

qap_fixed_terms <- function(mod, use_fixest) {
  rhs <- qap_fixed_rhs(mod, use_fixest)
  attr(stats::terms(stats::as.formula(call("~", rhs))), "term.labels")
}

#' @rdname qap_fixed_terms
#' @keywords internal
qap_drop_term <- function(mod, term, use_fixest) {
  if (use_fixest && is.call(mod[[3]]) && identical(mod[[3]][[1]], as.name("|"))) {
    keep <- setdiff(qap_fixed_terms(mod, use_fixest), term)
    rhs  <- if (length(keep)) paste(keep, collapse = " + ") else "1"
    fe   <- paste(deparse(mod[[3]][[3]], width.cutoff = 500), collapse = " ")
    return(stats::as.formula(paste(deparse(mod[[2]]), "~", rhs, "|", fe),
                             env = environment(mod)))
  }
  stats::update(mod, stats::as.formula(paste(". ~ . -", term)))
}

#' @rdname qap_fixed_terms
#' @keywords internal
qap_fixed_rhs <- function(mod, use_fixest) {
  rhs <- mod[[3]]
  if (use_fixest && is.call(rhs) && identical(rhs[[1]], as.name("|")))
    return(rhs[[2]])
  reformulas::nobars(mod)[[3]]
}


#' Change in classification metrics from including each term
#'
#' Refits the baseline model once per fixed-effect term with that term left
#' out -- no permutations and no residualisation -- and reports the full
#' model's metrics minus the reduced model's.  A positive value means the
#' term improves that metric.  The draws of every reduced model use the same
#' seed as the full model's, so the differences are not swamped by
#' Monte Carlo noise.
#'
#' Only terms involving a tested predictor get a row: \code{x1:x2} is kept
#' when \code{x1} or \code{x2} is in \code{tested}.  A GMM fit's design is
#' its predictors alone, so its terms are the tested predictors.
#'
#' @param full_cm The full model's \code{QAPConfusionMatrix}.
#' @param mod,pred,family,estimator,use_fixest,has_random,main_vars As for
#'   \code{fit_qap_model()}; \code{pred} is the frame the full model was
#'   fitted to, weights and offset columns included.
#' @param tested Character; the tested predictors.
#' @param seed Integer seed for the draws.
#' @inheritParams qap_binary_outcome
#'
#' @return A data frame with a \code{term} column, the eight differences,
#'   and \code{one_class}: TRUE when every fitted probability of the
#'   reduced model falls on the same side of 0.5, so that at the cut-off it
#'   predicts a single class.  NULL when the model has no terms to drop.
#' @keywords internal

qap_metric_change <- function(full_cm, mod, pred, dep, family, estimator,
                              use_fixest, has_random, main_vars, tested,
                              seed, kind = "glm", levels = c(0, 1)) {
  if (is.null(full_cm)) return(NULL)
  terms <- if (kind == "gmm") main_vars else qap_fixed_terms(mod, use_fixest)
  terms <- terms[vapply(strsplit(terms, ":", fixed = TRUE),
                        function(v) any(v %in% tested), logical(1))]
  if (!length(terms)) return(NULL)
  full <- qap_class_metrics(full_cm)

  failed <- character(0)
  rows <- lapply(terms, function(tm) {
    red_mod  <- qap_drop_term(mod, tm, use_fixest)
    red_vars <- intersect(main_vars,
                          all.vars(qap_fixed_rhs(red_mod, use_fixest)))
    # The refits repeat every warning the baseline fit already gave once
    # (separation, convergence), so they are kept quiet; a refit that fails
    # outright is reported below instead.
    cm <- tryCatch(suppressWarnings(suppressMessages({
      f <- fit_qap_model(mod = red_mod, pred = pred, family = family,
                         estimator = estimator, use_fixest = use_fixest,
                         main_vars = red_vars, has_random = has_random)
      qap_fit_confusion(f, seed, kind, pred, dep, red_vars, levels)
    })), error = function(e) NULL)
    if (is.null(cm)) {
      failed <<- c(failed, tm)
      return(c(rep(NA_real_, length(full)), NA))
    }
    c(full - qap_class_metrics(cm),
      any(colSums(cm$threshold_confusion_matrix) == 0))
  })

  if (length(failed))
    warning("The model without ", paste0("'", failed, "'", collapse = ", "),
            " could not be fitted; its metric changes are NA.", call. = FALSE)

  out <- data.frame(term = terms, do.call(rbind, rows),
                    stringsAsFactors = FALSE, row.names = NULL)
  names(out) <- c("term", names(full), "one_class")
  out$one_class <- as.logical(out$one_class)
  out
}


#' Literature note on LPM fitted values outside [0, 1]
#'
#' Shared by the once-per-session warning and the printed confusion matrix.
#'
#' @return A character string.
#' @keywords internal

qap_lpm_bias_note <- function() {
  paste("Linear probability models with fitted values outside [0, 1] can",
        "give biased and inconsistent coefficients (Horrace & Oaxaca, 2006,",
        "Economics Letters 90, 321-327); how many fall outside is, however,",
        "no reliable guide to how good the approximation is (Chen, Martin &",
        "Wooldridge, 2023, arXiv:2308.15338).")
}


#' Warn, once per session, about LPM fitted values outside [0, 1]
#'
#' @param cm A \code{QAPConfusionMatrix} or a named list of them.
#'
#' @return Invisibly NULL.
#' @keywords internal

qap_warn_lpm <- function(cm) {
  if (isTRUE(.qap_state$lpm_warned)) return(invisible(NULL))
  cms <- if (inherits(cm, "QAPConfusionMatrix")) list(cm) else cm
  n_out <- sum(vapply(cms, function(x)
    if (is.null(x$n_out)) 0 else as.numeric(x$n_out), numeric(1)))
  if (n_out > 0) {
    .qap_state$lpm_warned <- TRUE
    warning(format(n_out), " fitted value(s) of the linear probability ",
            "model lie outside [0, 1] and were truncated for the confusion ",
            "matrix. ", qap_lpm_bias_note(),
            " This warning is shown once per session.", call. = FALSE)
  }
  invisible(NULL)
}


#' Print a confusion matrix, or one per comparison
#'
#' @param cm A \code{QAPConfusionMatrix} or a named list of them.
#'
#' @return Invisibly NULL; prints.
#' @keywords internal

qap_print_confusion <- function(cm) {
  if (is.null(cm)) return(invisible(NULL))
  if (inherits(cm, "QAPConfusionMatrix")) {
    cat("\n")
    print(cm)
    return(invisible(NULL))
  }
  for (cn in names(cm)) {
    if (is.null(cm[[cn]])) next
    cat("\n--- Comparison: ", cn, " ---", sep = "")
    print(cm[[cn]])
  }
  invisible(NULL)
}


#' Print the change-in-metrics table, or one per comparison
#'
#' @param mc A data frame from \code{qap_metric_change()} or a named list
#'   of them.
#' @param digits Number of decimal places.
#'
#' @return Invisibly NULL; prints.
#' @keywords internal

qap_print_metric_change <- function(mc, digits = 4) {
  if (is.null(mc) || (is.list(mc) && !is.data.frame(mc) &&
                      all(vapply(mc, is.null, logical(1)))))
    return(invisible(NULL))
  starred <- FALSE
  one <- function(tab) {
    m   <- as.matrix(tab[, 2:9, drop = FALSE])
    out <- matrix(ifelse(is.na(m), "NA",
                         formatC(m, format = "f", digits = digits)),
                  nrow = nrow(m))
    # Tables from before the one_class column carry no stars.
    one_class <- if (is.null(tab$one_class)) rep(FALSE, nrow(m)) else
      tab$one_class %in% TRUE
    starred <<- starred || any(one_class)
    cut <- c(1, 3, 5, 7)
    out[, cut] <- paste0(out[, cut], ifelse(one_class, "*", " "))
    dimnames(out) <- list(tab$term,
                          paste0(rep(c("Acc", "Sens", "Spec", "Prec"),
                                     each = 2),
                                 c("(0.5)", "(draw)")))
    print(noquote(out), right = TRUE)
  }
  cat("\nChange in classification from including each term",
      "(full model minus model without it):\n")
  if (is.data.frame(mc)) {
    one(mc)
  } else {
    for (cn in names(mc)) {
      if (is.null(mc[[cn]])) next
      cat("\n--- Comparison: ", cn, " ---\n", sep = "")
      one(mc[[cn]])
    }
  }
  cat("(0.5) = 0.5 cut-off; (draw) = mean over 1000 Bernoulli draws from",
      "the fitted probabilities.\n")
  if (starred)
    cat(strwrap(paste(
      "* Every fitted probability of the model without this term lies on",
      "the same side of 0.5, so at the cut-off it predicts one class for",
      "all observations: sensitivity 1 and specificity 0 when that class is",
      "1, the reverse when it is 0.  Its accuracy is then the share of that",
      "class, and its precision the share of 1s (NA when it predicts only",
      "0s).  The (0.5) changes in sensitivity and specificity therefore",
      "trade off against each other; read them together, or use the (draw)",
      "columns, where a model without information scores at chance.")),
      sep = "\n")
  invisible(NULL)
}
