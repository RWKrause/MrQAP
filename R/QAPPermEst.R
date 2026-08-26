#' Estimate one permutation
#'
#' Called by \code{future.apply::future_lapply()} from
#' \code{qap_cpu_perms()}.  Permutes the data, rebuilds the vectorised model
#' frame, refits, and compares the result to the baseline.
#'
#' Handles both network matrices and CSS arrays: the only differences are
#' which builder vectorises the data and whether \code{RMPerm()} applies one
#' node order across all three dimensions.
#'
#' @param i Integer; iteration index (unused; required by
#'   \code{future_lapply}).
#' @param data. Named list of matrices or 3D arrays.
#' @param perm_var. NULL for qapy (permute the outcome), or the name of the
#'   already-residualised variable to permute for qapspp.
#' @param mode. Character; \code{"digraph"}/\code{"graph"} for matrices,
#'   \code{"directed"}/\code{"undirected"} for CSS.
#' @param diag. Logical; include the diagonal?
#' @param mod. Formula for the model.
#' @param groups. Permutation grouping.
#' @param fit. Baseline fit, or list of fits when comparisons are used.
#' @param family. Model family.
#' @param estimator. Estimator type.
#' @param use_fixest. Logical.
#' @param fixest_se_cluster. Cluster variable or NULL.
#' @param use_robust_errors. Logical.
#' @param has_random. Logical.
#' @param main_vars. Character vector of main predictor names.
#' @param data_vars. Character vector of all data variable names.
#' @param parsed. Output of \code{parse_qap_formula()}.
#' @param comp. Comparison list or NULL.
#' @param reference. Reference category or NULL.
#' @param css. Logical; CSS data?
#' @param multi_mode. Logical; each dimension its own node set?
#'
#' @return A list with lower/larger/abs comparison results, or NULL if this
#'   permutation could not be used.
#' @keywords internal

QAPPermEst <- function(i,
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
                       reference.,
                       css. = FALSE,
                       multi_mode. = FALSE) {

  dep   <- parsed.$dependent
  large <- is.list(data.[[dep]])

  # Categories present in the unpermuted outcome, for the sufficiency check.
  y_cat <- na.omit(unique(as.vector(unlist(data.[[dep]]))))

  build <- function(d) {
    if (!large) {
      xs <- lapply(data_vars., function(v) d[[v]])
      names(xs) <- data_vars.
      p <- if (css.) {
        make_css_data(y = d[[dep]], x = xs, nets = 1,
                      diag = diag., mode = mode.,
                      multi_mode = multi_mode.)$pred
      } else {
        make_qap_data(y = d[[dep]], x = xs, g = groups., diag = diag.,
                      mode = mode., net = 1, perm = FALSE, xi = NULL,
                      multi_mode = multi_mode.)
      }
    } else {
      parts <- vector("list", length(d[[dep]]))
      for (net in seq_along(d[[dep]])) {
        xs <- lapply(data_vars., function(v) d[[v]][[net]])
        names(xs) <- data_vars.
        g <- if (is.list(groups.)) groups.[[net]] else groups.
        parts[[net]] <- if (css.) {
          make_css_data(y = d[[dep]][[net]], x = xs, nets = net,
                        diag = diag., mode = mode.,
                        multi_mode = multi_mode.)$pred
        } else {
          make_qap_data(y = d[[dep]][[net]], x = xs, g = g, diag = diag.,
                        mode = mode., net = net, perm = FALSE, xi = NULL,
                        multi_mode = multi_mode.)
        }
      }
      p <- do.call(rbind, parts)
    }
    names(p)[names(p) == "yv"] <- dep
    p
  }

  # A permuted sample can be degenerate -- a constant outcome, a constant
  # predictor, or perfect collinearity inside a comparison.  Redraw rather
  # than feed the model something it cannot fit.  Well-conditioned data pass
  # on the first try, so this consumes no extra randomness.
  max_trials <- 10000L
  sufficient <- FALSE
  pred <- NULL

  for (trial in seq_len(max_trials)) {
    d <- data.
    target <- if (is.null(perm_var.)) dep else perm_var.
    d[[target]] <- perm_networks(d[[target]], groups., CSS = css.,
                                 multi_mode = multi_mode.)

    pred <- build(d)

    if (family. != "multinom" && is.null(comp.)) {
      y_ok <- length(na.omit(unique(pred[[dep]]))) > 1
    } else {
      y2_cat    <- na.omit(unique(pred[[dep]]))
      y_present <- all(y_cat %in% y2_cat)
      y_mult    <- all(table(pred[[dep]]) > 2)
      y_ok      <- y_present && y_mult
    }

    x_ok <- TRUE
    num_preds <- pred[, data_vars.[data_vars. %in% names(pred)], drop = FALSE]
    num_preds <- num_preds[, vapply(num_preds, is.numeric, logical(1)),
                           drop = FALSE]
    if (ncol(num_preds) > 0) {
      x_ok <- all(vapply(num_preds,
                         function(col) length(unique(col)) > 1, logical(1)))
    }

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

    sufficient <- y_ok && x_ok
    if (sufficient) break
  }

  # Every other failure mode here returns NULL and is filtered by
  # qap_aggregate(), which reports how many permutations were dropped.
  if (!sufficient) return(NULL)

  fit_one <- function(p) {
    tryCatch(
      fit_qap_model(mod = mod., pred = p, family = family.,
                    estimator = estimator., use_fixest = use_fixest.,
                    fixest_se_cluster = fixest_se_cluster.,
                    use_robust_errors = use_robust_errors.,
                    main_vars = main_vars., has_random = has_random.,
                    reference = reference.),
      error = function(e) NULL
    )
  }

  xi_arg <- perm_var.

  if (is.null(comp.)) {
    pf <- fit_one(pred)
    if (is.null(pf)) return(NULL)
    return(compare_perm_to_baseline(pf$coefficients, pf$t, fit., xi = xi_arg))
  }

  out <- vector("list", length(comp.))
  names(out) <- names(comp.)
  for (k in seq_along(comp.)) {
    pk <- pred[pred[[dep]] %in% comp.[[k]], ]
    pk[[dep]] <- ifelse(pk[[dep]] == comp.[[k]][1], 0, 1)

    pf <- fit_one(pk)
    if (is.null(pf)) return(NULL)

    out[[k]] <- compare_perm_to_baseline(pf$coefficients, pf$t,
                                         fit.[[k]], xi = xi_arg)
  }
  out
}
