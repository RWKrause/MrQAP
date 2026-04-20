#' Print QAPglm() results for generalized linear regressions
#'
#' @param x An object of class \code{QAPGLM}.
#' @param print_b Logical; also show parameter-based p-values?
#' @param print_random Logical; show random intercepts?
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export

print.QAPGLM <- function(x, ...,
                          print_b = FALSE,
                          print_random = FALSE) {
  stopifnot(inherits(x, "QAPGLM"))

  # --- header ---
  if (is.null(x$random.intercepts)) {
    cat("\nGeneralized Linear Network Model\n")
  } else {
    cat("\nGeneralized Linear Mixed Network Model fit by REML\n")
  }
  if (!is.null(x$estimator) && x$estimator == "gmm")
    cat("\nEstimator: Generalized Method-of-Moments.")
  if (!is.null(x$theta))
    cat("\nNegative binomial dispersion (theta):", format(round(x$theta, 4)))
  if (!is.null(x$zi_coefficients)) {
    cat("\nZero-inflation coefficients:")
    cat("\n  ", paste(names(x$zi_coefficients),
                      format(round(x$zi_coefficients, 4)),
                      sep = " = ", collapse = ", "))
  }
  if (!is.null(x$groups))
    cat("\nPermutations were performed within groups only.")

  if (x$nullhyp == "qapy")
    cat("\nThe outcome matrix Y was permuted", format(x$reps), "times.")
  if (x$nullhyp == "qapspp") {
    cat("\nSignificance was estimated using Dekker's")
    cat("\n  'semi-partialling plus' procedure with",
        format(x$reps), "permutations.")
  }

  if (x$diag) {
    cat("\nDiagonal values (loops) were used in the estimation.")
  } else {
    cat("\nDiagonal values (loops) were ignored.")
  }
  cat("\nThe outcome was treated as", format(paste0(x$mode, ".")))
  cat("\nModel family:", format(x$family))

  # --- comparisons ---
  if (!is.null(x$comp)) {
    for (k in seq_along(x$comp)) {
      cat("\n\n--- Comparison:", names(x$comp)[k], "---")
      cat("\n   ", x$comp[[k]][1], "vs", x$comp[[k]][2])
      .print_glm_table(x$base[[k]], x$lower[[k]], x$larger[[k]], x$abs[[k]],
                        x$nullhyp, print_b)
    }
  } else {
    cat("\n\nCoefficients:\n")
    if (print_b) {
      cmat <- matrix(NA, nrow = length(x$coefficients), ncol = 5)
      cmat[, 1] <- format(as.numeric(x$coefficients))
      cmat[, 2] <- format(exp(as.numeric(x$coefficients)))
      cmat[, 3] <- format(x$lower[1, ])
      cmat[, 4] <- format(x$larger[1, ])
      cmat[, 5] <- format(x$abs[1, ])
      if (x$nullhyp == "qapspp") cmat[1, 3:5] <- "*"
      colnames(cmat) <- c("Estimate", "Exp(b)", "Pr(<=b)", "Pr(>=b)", "Pr(>=|b|)")
      rownames(cmat) <- names(x$coefficients)
      print.table(cmat)
      cat("--------------\n")
    }

    cmat <- matrix(NA, nrow = length(x$coefficients), ncol = 5)
    cmat[, 1] <- format(as.numeric(x$coefficients))
    cmat[, 2] <- format(exp(as.numeric(x$coefficients)))
    cmat[, 3] <- format(x$lower[2, ])
    cmat[, 4] <- format(x$larger[2, ])
    cmat[, 5] <- format(x$abs[2, ])
    if (x$nullhyp == "qapspp") cmat[1, 3:5] <- "*"
    colnames(cmat) <- c("Estimate", "Exp(b)", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
    rownames(cmat) <- names(x$coefficients)
    print.table(cmat)
  }

  if (x$nullhyp == "qapspp")
    cat("\n* Significance test for the intercept is not defined with qapspp.\n")

  cat("--------------\n")

  if (!is.null(x$random.intercepts) && print_random) {
    cmat <- matrix(round(x$random.intercepts, 3),
                   nrow = length(x$random.intercepts), ncol = 1)
    colnames(cmat) <- "Random Intercepts:"
    print.table(cmat)
    cat("--------------\n")
  }

  if (!is.null(x$simple_fit) && !is.null(x$estimator) &&
      x$estimator != "gmm") {
    cat("\nAIC:", format(AIC(x$simple_fit)))
    cat("\nBIC:", format(BIC(x$simple_fit)))
  }

  # --- confusion matrix ---
  if (!is.null(x$confusion_matrix)) {
    cat("\n\n")
    print(x$confusion_matrix)
  }
  cat("\n")
  invisible(x)
}

# Helper for comparison tables
.print_glm_table <- function(base, lower, larger, abs_mat, nullhyp, print_b) {
  cat("\n\nCoefficients:\n")
  nc <- length(base$coefficients)
  cmat <- matrix(NA, nrow = nc, ncol = 4)
  cmat[, 1] <- format(round(as.numeric(base$coefficients), 4))
  cmat[, 2] <- format(lower[2, ])
  cmat[, 3] <- format(larger[2, ])
  cmat[, 4] <- format(abs_mat[2, ])
  if (nullhyp == "qapspp") cmat[1, 2:4] <- "*"
  colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
  rownames(cmat) <- names(base$coefficients)
  print.table(cmat)
}
