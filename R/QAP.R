#' Parameter and p-value estimation with MRQAP
#'
#' Estimates regression models on network data with permutation-based
#' inference.  Handles both ordinary network matrices and cubic Cognitive
#' Social Structure (CSS) arrays of format \code{[sender, receiver,
#' perceiver]}; which one you have is detected from the data.
#'
#' @param formula A formula describing the model.
#'   The left-hand side is the dependent variable name in \code{data};
#'   the right-hand side lists predictor variable names.
#'   Use \code{|} (without parentheses) for fixest fixed effects:
#'   \code{y ~ x1 + x2 | fe1}.
#'   Use \code{(1|var)} for lme4-style random effects: \code{y ~ x1 + (1|sv)}.
#'
#' @param data Named list.  Each element is either a square \code{matrix}
#'   (network data) or a cubic \code{array} (CSS data), or a \code{list} of
#'   those if there are several networks.  Names must correspond to the
#'   variables in \code{formula}.
#'
#' @param family Character; model family (default \code{"gaussian"}).
#'   Supported: \code{"gaussian"}, \code{"binomial"}, \code{"poisson"},
#'   \code{"negbin"} (negative binomial), \code{"zip"} (zero-inflated
#'   Poisson) and \code{"multinom"} (multinomial).
#'
#' @param mode Character; \code{"directed"} (default) or \code{"undirected"}.
#'   In undirected mode each dyad contributes a single observation.
#'
#' @param diag Logical; include diagonal values (self-ties, default
#'   \code{FALSE}).
#'
#' @param nullhyp Character; \code{"qapspp"} (default, recommended) or
#'   \code{"qapy"}.  \code{"qapspp"} is Dekker's semi-partialling plus
#'   procedure, which residualises each predictor on the others before
#'   permuting it; \code{"qapy"} permutes the outcome.  With a single
#'   predictor the two coincide and \code{"qapy"} is used.
#'
#' @param estimator Character; \code{"standard"} (default) or \code{"gmm"}
#'   (for binomial, poisson, negbin and zip families).
#'
#' @param reps Integer; number of permutations (default 1000).  Under
#'   \code{"qapspp"} this many permutations are run \emph{per predictor}.
#' @param seed Integer; optional random seed.  The RNG state is restored on
#'   exit.
#' @param groups Vector, or list of vectors (one per network); restricts
#'   permutation to occur only within groups.
#' @param ncores Integer; number of cores for parallel processing via the
#'   \code{future} framework.  When \code{NULL} the active \code{future}
#'   plan is used unchanged.
#' @param fixest_se_cluster Character; cluster variable for fixest.
#'
#' @param comparison Named list of length-2 character vectors specifying
#'   category comparisons, e.g. \code{list(commission =
#'   c("false_positive", "true_negative"))}.  Each comparison is estimated
#'   and reported separately.
#' @param reference Character; reference category (multinomial, or
#'   comparison models).
#'
#' @param random_intercept_nets,random_intercept_sender,random_intercept_receiver,random_intercept_perceiver
#'   Logical; add random intercepts for networks, senders, receivers and
#'   (CSS only) perceivers.
#'
#' @param use_robust_errors Logical; use heteroskedasticity-consistent
#'   standard errors.  HC3 for linear models, the GLM sandwich otherwise.
#'   Not available for mixed models.
#' @param less_mem Logical; do not store the fitted model object.
#' @param use_gpu Logical; use GPU-accelerated batch OLS (gaussian only,
#'   requires the \code{torch} package).
#'
#' @param css Logical or \code{NULL}; is this CSS (3D array) data?
#'   \code{NULL} (default) detects it from the dimensions of the dependent
#'   variable.  Supplying \code{TRUE} or \code{FALSE} asserts it, and errors
#'   if the data do not match.
#'
#' @param multi_mode Logical or \code{NULL}; does each dimension index its
#'   own node set (2- or 3-mode data, as produced by
#'   \code{\link{df_to_mat}(multi_mode = TRUE)})?  Each dimension is then
#'   permuted independently.  \code{NULL} (default) infers \code{TRUE} when
#'   the dimensions differ in size and \code{FALSE} otherwise -- a square
#'   matrix is genuinely ambiguous, so state it explicitly if your node sets
#'   happen to be the same size.  Self-ties and undirected mode are
#'   undefined for multi-mode data and are rejected.
#'
#' @return An object of class \code{QAPCSS} (CSS data), \code{QAPRegression}
#'   (gaussian, no comparisons) or \code{QAPGLM} (everything else).
#'
#' @references
#' Dekker, D., Krackhardt, D., & Snijders, T. A. B. (2007).
#' Sensitivity of MRQAP tests to collinearity and autocorrelation conditions.
#' \emph{Psychometrika}, 72(4), 563--581. \doi{10.1007/s11336-007-9016-1}
#'
#' Krackhardt, D. (1988). Predicting with networks: Nonparametric multiple
#' regression analysis of dyadic data. \emph{Social Networks}, 10(4),
#' 359--381. \doi{10.1016/0378-8733(88)90004-4}
#'
#' @seealso \code{\link{df_to_mat}} to build the input from a dyadic data
#'   frame, and \code{\link{combine_qap_estimates}} to pool several runs.
#'
#' @export
#'
#' @examples
#' # --- a one-mode network ------------------------------------------------
#' set.seed(1)
#' n  <- 12
#' x1 <- matrix(rnorm(n^2), n, n)
#' x2 <- matrix(rnorm(n^2), n, n)
#' y  <- 1 + 0.8 * x1 - 0.4 * x2 + matrix(rnorm(n^2), n, n)
#'
#' fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2), reps = 50)
#' fit
#'
#' # --- a cognitive social structure --------------------------------------
#' # A 3D array [sender, receiver, perceiver] is detected automatically.
#' a1 <- array(rnorm(6^3), dim = c(6, 6, 6))
#' ya <- 0.7 * a1 + array(rnorm(6^3), dim = c(6, 6, 6))
#' QAP(y ~ x1, data = list(y = ya, x1 = a1), reps = 20)
#'
#' # --- two-mode data: rows and columns are different node sets ------------
#' b1 <- matrix(rnorm(9 * 5), 9, 5)
#' yb <- 0.6 * b1 + matrix(rnorm(9 * 5), 9, 5)
#' QAP(y ~ x1, data = list(y = yb, x1 = b1), reps = 20)

QAP <- function(formula,
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
                random_intercept_nets       = FALSE,
                random_intercept_sender     = FALSE,
                random_intercept_receiver   = FALSE,
                random_intercept_perceiver  = FALSE,
                use_robust_errors = FALSE,
                less_mem   = FALSE,
                use_gpu    = FALSE,
                css        = NULL,
                multi_mode = NULL) {

  mode      <- match.arg(mode, c("directed", "undirected"))
  nullhyp   <- match.arg(nullhyp, c("qapspp", "qapy"))
  estimator <- match.arg(estimator, c("standard", "gmm"))
  family    <- match.arg(family, c("gaussian", "binomial", "poisson",
                                   "negbin", "zip", "multinom"))

  # --- seed: do not leave the caller's RNG state altered ---
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
              add = TRUE)
    }
    set.seed(seed)
  }

  # --- parse formula ---
  parsed    <- parse_qap_formula(formula, fixest_se_cluster)
  dep       <- parsed$dependent
  main      <- parsed$main
  # Structural vars (sv, rv, nv, pv) are auto-generated, not user data.
  data_vars <- intersect(parsed$all_data_vars, names(data))
  nx        <- length(main)

  if (!(dep %in% names(data)))
    stop("Dependent variable '", dep, "' not found in data.")

  # --- matrix or CSS? ---
  css <- detect_css(data[[dep]], css = css, dep = dep)
  multi_mode <- detect_multi_mode(data[[dep]], multi_mode = multi_mode,
                                  dep = dep)

  validate_qap_input(data, parsed, css = css)
  large <- is.list(data[[dep]])

  if (multi_mode) {
    # Neither self-ties nor reciprocity exist between disjoint node sets.
    if (diag)
      stop("diag = TRUE is not meaningful for multi-mode data: there are no ",
           "self-ties between distinct node sets.")
    if (mode == "undirected")
      stop("mode = \"undirected\" is not meaningful for multi-mode data: ",
           "there is no reciprocal dyad across distinct node sets.")
  }

  rin <- random_intercept_nets
  ris <- random_intercept_sender
  rir <- random_intercept_receiver
  rip <- random_intercept_perceiver

  if (rip && !css) {
    warning("random_intercept_perceiver applies to CSS data only; ignored.")
    rip <- FALSE
  }

  # --- validation that must precede formula construction -------------
  # These switch random-effect flags off.  build_internal_formula() has to
  # see the final values, or the formula keeps terms the warnings claim
  # were dropped.
  if (!is.null(reference) && !is.character(reference) && family == "multinom")
    reference <- as.character(reference)
  if (use_robust_errors && family == "multinom") {
    warning("Robust SEs not implemented for multinomial.")
    use_robust_errors <- FALSE
  }
  if (nullhyp == "qapspp" && nx == 1) nullhyp <- "qapy"
  if (mode == "undirected" && (ris || rir)) {
    warning("Undirected mode: sender/receiver random intercepts set to FALSE.")
    ris <- rir <- FALSE
  }
  if (family == "multinom" && (rin || ris || rir || rip || parsed$has_random)) {
    warning("Random intercepts not implemented for multinomial. ",
            "Using standard nnet::multinom().")
    rin <- ris <- rir <- rip <- FALSE
  }
  if (diag) warning("Results may not be valid when the diagonal is used.")

  # --- internal formula ---
  mod     <- build_internal_formula(formula, rin = rin, ris = ris,
                                    rir = rir, rip = rip)
  mod_str <- paste(deparse(mod, width.cutoff = 500), collapse = " ")
  has_random <- grepl("\\(", mod_str) || parsed$has_random
  if (family == "multinom") has_random <- FALSE
  use_fixest <- parsed$use_fixest
  if (has_random && use_fixest) {
    warning("Cannot combine fixest FE and lme4 random effects. ",
            "Using lme4 random effects only.")
    use_fixest <- FALSE
  }
  mod <- as.formula(mod_str)
  if (family == "multinom" && parsed$has_random) {
    # nnet::multinom() cannot parse (1|g) terms.
    mod <- reformulas::nobars(mod)
  }

  use_robust_errors <- resolve_robust_errors(use_robust_errors, family,
                                             estimator, has_random)

  # Random-effect suffix used when residualising for qapspp.
  rand_part <- ""
  if (rin) rand_part <- paste(rand_part, "+ (1|nv)")
  if (rip) rand_part <- paste(rand_part, "+ (1|pv)")
  if (ris) rand_part <- paste(rand_part, "+ (1|sv)")
  if (rir) rand_part <- paste(rand_part, "+ (1|rv)")

  # --- groups ---
  gp      <- normalise_groups(groups, data[[dep]], large = large)
  groups  <- gp$groups
  grouped <- gp$grouped

  # --- vectorise the data ---
  built <- qap_build_pred(data, dep, data_vars, groups, diag, mode,
                          css = css, large = large, multi_mode = multi_mode)
  pred       <- built$pred
  valid      <- built$valid
  valid_list <- built$valid_list

  # --- can the accelerated paths be used? ---
  # Both rest on the same assumptions, so they share one eligibility test.
  # Resolved before the baseline fit so an unusable use_gpu is reported
  # immediately rather than after the expensive part.
  fast_ok <- qap_ols_eligible(family, estimator, has_random, use_fixest,
                              comparison, data, c(dep, main))
  gpu_ok  <- use_gpu && fast_ok && !use_robust_errors

  if (use_gpu && !gpu_ok)
    warning("use_gpu = TRUE was requested but this model cannot use the GPU ",
            "path (it is limited to gaussian models with complete data, no ",
            "random or fixed effects, no comparisons and no robust standard ",
            "errors). Falling back to the CPU.")

  # --- baseline fit ---
  fit <- list()
  fit$base <- qap_fit_baseline(mod, pred, dep, family, estimator, use_fixest,
                               fixest_se_cluster, use_robust_errors, main,
                               has_random, reference, comparison)

  # --- permutations ---
  if (fast_ok) {
    tmpl <- qap_ols_template(data, dep, main, diag, mode, css, large,
                             multi_mode)
  }

  if (gpu_ok) {
    agg <- qap_gpu_perms(tmpl, data, dep, main, groups, reps, fit$base,
                         nullhyp, css, multi_mode, pred, valid, valid_list,
                         large, rand_part, has_random, mode,
                         use_robust_errors = use_robust_errors)
  } else if (fast_ok) {
    # Plain OLS with complete data: skip the formula interface and the
    # per-permutation data frame entirely.
    agg <- qap_ols_perms(tmpl, data, dep, main, groups, reps, fit$base,
                         nullhyp, css, multi_mode, pred, valid, valid_list,
                         large, rand_part, has_random, mode, diag, ncores,
                         use_robust_errors = use_robust_errors)
  } else {
    agg <- qap_cpu_perms(data, parsed, mode, diag, groups, reps, fit$base,
                         nullhyp, main, data_vars, pred, valid, valid_list,
                         large, rand_part, mod, family, estimator,
                         use_fixest, fixest_se_cluster, use_robust_errors,
                         has_random, comparison, reference, ncores, css,
                         multi_mode)
  }
  fit$lower  <- agg$lower
  fit$larger <- agg$larger
  fit$abs    <- agg$abs
  # How many permutations actually contributed -- the divisor of the three
  # proportions above. A scalar under qapy; one entry per predictor under
  # qapspp, which tests each in its own set of permutations. Used by
  # combine_qap_estimates() to pool by weight rather than by intent.
  fit$n_valid <- agg$n_valid
  # The permutation draws themselves, used by confint(). Cheap to keep
  # (reps x k doubles) but suppressed by less_mem like the model object.
  if (!less_mem) fit$null_dist <- agg$draws

  # --- package results ---
  # The object has the same shape whatever the input shape was.
  if (is.null(comparison)) {
    fit$coefficients <- fit$base$coefficients
    fit$t            <- fit$base$t
    for (nm in c("r.squared", "adj.r.squared", "random.intercepts",
                 "theta", "zi_coefficients")) {
      if (!is.null(fit$base[[nm]])) fit[[nm]] <- fit$base[[nm]]
    }
    if (!less_mem) fit$simple_fit <- fit$base$base_model
  } else if (!less_mem) {
    fit$simple_fits <- lapply(fit$base, `[[`, "base_model")
  }

  if (family == "binomial" && is.null(comparison)) {
    bm <- fit$base$base_model
    if (!inherits(bm, "gmm")) {
      fit$confusion_matrix <- probabilistic_confusion_matrix(
        actual = pred[[dep]], predicted_prob = fitted(bm),
        n_draws = 1000, seed = seed
      )
    }
  }

  fit$nullhyp   <- nullhyp
  fit$diag      <- diag
  fit$family    <- family
  fit$mode      <- mode
  fit$reps      <- reps
  fit$css       <- css
  fit$multi_mode <- multi_mode
  fit$groups    <- if (grouped) unique(unlist(groups)) else NULL
  fit$random    <- c(sender = ris, receiver = rir,
                     perceiver = rip, nets = rin)
  # GMM always reports heteroskedasticity-robust (MDS) standard errors,
  # whether or not use_robust_errors was set.
  fit$robust_se <- use_robust_errors || estimator == "gmm"
  fit$estimator <- estimator
  fit$comp      <- comparison
  fit$reference <- reference

  # The specific class drives the print method; the shared "QAP" parent
  # carries coef(), summary(), confint() and the rest, which are identical
  # across all three.
  class(fit) <- c(
    if (css) "QAPCSS"
    else if (family == "gaussian" && is.null(comparison)) "QAPRegression"
    else "QAPGLM",
    "QAP")
  fit
}


#' Decide whether the data are CSS
#'
#' CSS data are cubic arrays \code{[sender, receiver, perceiver]}; network
#' data are matrices.  The distinction is read off the number of dimensions.
#'
#' @param y The dependent variable from \code{data}, possibly a list.
#' @param css Logical or NULL; NULL auto-detects, TRUE/FALSE asserts.
#' @param dep Character; the dependent variable's name, for messages.
#'
#' @return Logical.
#' @keywords internal

detect_css <- function(y, css = NULL, dep = "y") {
  probe <- if (is.list(y)) y[[1]] else y
  nd    <- length(dim(probe))
  found <- nd == 3L

  if (is.null(css)) {
    if (!nd %in% c(2L, 3L))
      stop("data[['", dep, "']] must be a matrix or a 3-dimensional array; ",
           "it has ", nd, " dimension(s).")
    return(found)
  }

  if (!is.logical(css) || length(css) != 1L || is.na(css))
    stop("`css` must be TRUE, FALSE, or NULL.")
  if (css && !found)
    stop("css = TRUE, but data[['", dep, "']] has ", nd,
         " dimension(s); CSS data must be a 3-dimensional array ",
         "[sender, receiver, perceiver].")
  if (!css && found)
    stop("css = FALSE, but data[['", dep, "']] is a 3-dimensional array. ",
         "Drop the argument to detect this automatically.")
  css
}


#' Decide whether each dimension is its own node set
#'
#' Dimensions of different sizes can only be distinct node sets.  Equal sizes
#' are genuinely ambiguous -- a square matrix may be one-mode, or two-mode
#' with equally large node sets -- so the safe default is one-mode and the
#' caller can override.
#'
#' @param y The dependent variable from \code{data}, possibly a list.
#' @param multi_mode Logical or NULL; NULL infers, TRUE/FALSE asserts.
#' @param dep Character; the dependent variable's name, for messages.
#'
#' @return Logical.
#' @keywords internal

detect_multi_mode <- function(y, multi_mode = NULL, dep = "y") {
  probe  <- if (is.list(y)) y[[1]] else y
  d      <- dim(probe)
  ragged <- length(unique(d)) != 1L

  if (is.null(multi_mode)) return(ragged)

  if (!is.logical(multi_mode) || length(multi_mode) != 1L || is.na(multi_mode))
    stop("`multi_mode` must be TRUE, FALSE, or NULL.")
  if (!multi_mode && ragged)
    stop("multi_mode = FALSE, but data[['", dep, "']] has dimensions ",
         paste(d, collapse = " x "),
         ". One-mode data must be square (or cubic).")
  multi_mode
}


#' Normalise the permutation grouping
#'
#' @param groups Vector; an unnamed list of vectors (one per network); a
#'   named list of per-dimension groupings (multi-mode); or NULL.
#' @param y The dependent variable, used for the expected node counts.
#' @param large Logical; is y a list of networks?
#'
#' @return A list with \code{groups} and \code{grouped} (whether the user
#'   supplied any).
#' @keywords internal

normalise_groups <- function(groups, y, large = FALSE) {
  grouped <- !is.null(groups)
  nodes   <- function(a) dim(a)[1]

  if (!grouped) return(list(groups = NULL, grouped = FALSE))

  # A NAMED list is one grouping per dimension (multi-mode).  It applies to
  # every network, so it is passed through untouched; rmperm_dim_groups()
  # validates the per-dimension lengths against the actual dimensions.
  if (is.list(groups) && !is.null(names(groups)))
    return(list(groups = groups, grouped = TRUE))

  if (!large) {
    n <- nodes(y)
    if (length(groups) != n)
      stop("groups has length ", length(groups), " but the network has ",
           n, " nodes.")
    return(list(groups = as.factor(groups), grouped = TRUE))
  }

  ns <- vapply(y, nodes, integer(1))

  if (is.list(groups)) {
    if (length(groups) != length(y))
      stop("groups is a list of length ", length(groups), " but there are ",
           length(y), " networks.")
    for (i in seq_along(groups)) {
      if (length(groups[[i]]) != ns[i])
        stop("groups[[", i, "]] has length ", length(groups[[i]]),
             " but network ", i, " has ", ns[i], " nodes.")
    }
    return(list(groups = lapply(groups, as.factor), grouped = TRUE))
  }

  if (any(ns != length(groups)))
    stop("groups has length ", length(groups),
         " but the networks do not all have that many nodes; ",
         "supply one grouping vector per network instead.")
  list(groups = as.factor(groups), grouped = TRUE)
}


#' Vectorise network or CSS data into a model data frame
#'
#' Single dispatch point for the matrix and CSS builders, and for the
#' single- and multi-network cases.
#'
#' @param data Named list of matrices or arrays.
#' @param dep Character; dependent variable name.
#' @param data_vars Character; variables to vectorise.
#' @param groups Normalised grouping.
#' @param diag Logical; include the diagonal?
#' @param mode Character; "directed" or "undirected".
#' @param css Logical; CSS data?
#' @param large Logical; a list of networks?
#' @param multi_mode Logical; each dimension its own node set?
#'
#' @return A list with \code{pred}, and \code{valid}/\code{valid_list} for CSS.
#' @keywords internal

qap_build_pred <- function(data, dep, data_vars, groups, diag, mode,
                           css = FALSE, large = FALSE, multi_mode = FALSE) {
  mode_internal <- if (mode == "directed") "digraph" else "graph"

  one <- function(y, xs, g, net) {
    if (css) {
      make_css_data(y = y, x = xs, nets = net, diag = diag, mode = mode,
                    multi_mode = multi_mode)
    } else {
      list(pred = make_qap_data(y = y, x = xs, g = g, diag = diag,
                                mode = mode_internal, net = net,
                                perm = FALSE, xi = NULL,
                                multi_mode = multi_mode),
           valid = NULL)
    }
  }

  if (!large) {
    xs <- data[data_vars]
    b  <- one(data[[dep]], xs, groups, 1)
    pred <- b$pred
    names(pred)[names(pred) == "yv"] <- dep
    return(list(pred = pred, valid = b$valid, valid_list = NULL))
  }

  nets       <- seq_along(data[[dep]])
  pred_list  <- vector("list", length(nets))
  valid_list <- vector("list", length(nets))
  for (net in nets) {
    xs <- lapply(data_vars, function(v) data[[v]][[net]])
    names(xs) <- data_vars
    g  <- if (is.list(groups)) groups[[net]] else groups
    b  <- one(data[[dep]][[net]], xs, g, net)
    pred_list[[net]]  <- b$pred
    valid_list[[net]] <- b$valid
  }
  pred <- do.call(rbind, pred_list)
  names(pred)[names(pred) == "yv"] <- dep
  list(pred = pred, valid = NULL, valid_list = valid_list)
}


#' Fit the baseline model, with or without comparisons
#'
#' @inheritParams fit_qap_model
#' @param dep Character; dependent variable name.
#' @param comparison Comparison list or NULL.
#'
#' @return A fit list, or a named list of them when comparisons are used.
#' @keywords internal

qap_fit_baseline <- function(mod, pred, dep, family, estimator, use_fixest,
                             fixest_se_cluster, use_robust_errors, main_vars,
                             has_random, reference, comparison) {
  fit1 <- function(p) {
    f <- fit_qap_model(mod = mod, pred = p, family = family,
                       estimator = estimator, use_fixest = use_fixest,
                       fixest_se_cluster = fixest_se_cluster,
                       use_robust_errors = use_robust_errors,
                       main_vars = main_vars, has_random = has_random,
                       reference = reference, baseline = TRUE)
    warn_aliased_coefs(f)
    f
  }

  if (is.null(comparison)) return(fit1(pred))

  out <- vector("list", length(comparison))
  names(out) <- names(comparison)
  for (k in seq_along(comparison)) {
    pk <- pred[pred[[dep]] %in% comparison[[k]], ]
    pk[[dep]] <- ifelse(pk[[dep]] == comparison[[k]][1], 0, 1)
    out[[k]] <- fit1(pk)
  }
  out
}
