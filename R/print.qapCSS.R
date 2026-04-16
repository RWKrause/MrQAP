#' Print QAPcss() results
#'
#' @param x An object of class \code{QAPCSS}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export

print.QAPCSS <- function(x, ...) {

  # --- header ---
  if (x$family != "multinom") {
    if (!any(x$random)) {
      cat("\nGeneralized Linear Network Model for CSS\n\n")
    } else {
      cat("\nGeneralized Linear Mixed Network Model for CSS fit by REML\n\n")
    }
  } else {
    cat("\nMultinomial Choice Network Model for CSS\n\n")
    cat("The reference group was", format(paste0(x$reference, ".")), "\n")
  }

  if (!is.null(x$estimator) && x$estimator == "gmm")
    cat("Estimator: Generalized Method-of-Moments.\n")
  if (!is.null(x$theta))
    cat("Negative binomial dispersion (theta):", format(round(x$theta, 4)), "\n")
  if (!is.null(x$zi_coefficients)) {
    cat("Zero-inflation coefficients:\n")
    cat("  ", paste(names(x$zi_coefficients),
                    format(round(x$zi_coefficients, 4)),
                    sep = " = ", collapse = ", "), "\n")
  }

  if (!is.null(x$groups))
    cat("Permutations were performed within groups only.\n")

  if (x$nullhyp == "qapy")
    cat("The outcome array Y was permuted", format(x$reps), "times.\n")
  if (x$nullhyp == "qapspp") {
    cat("Significance was estimated using Dekker's\n")
    cat("  'semi-partialling plus' procedure with",
        format(x$reps), "permutations.\n")
  }

  if (x$robust_se)
    cat("T-values are based on robust standard errors.\n")

  if (x$diag) {
    cat("Diagonal values (loops) were used in the estimation.\n",
        "  Results may be biased because of that.\n")
  } else {
    cat("Diagonal values (loops) were ignored.\n")
  }
  cat("The outcome was treated as", format(paste0(x$mode, ".")), "\n")

  # --- results ---
  if (x$family != "multinom") {
    if (is.null(x$comp)) {
      glm_tab(x, comp = x$comp)
    } else {
      for (mod in seq_along(x$comp)) {
        glm_tab(x, comp = names(x$comp)[[mod]])
      }
    }
  } else {
    cat("\nCoefficients:\n\n")
    for (option in seq_len(nrow(x$base$coefficients))) {
      cat(format(paste0("-- ", rownames(x$base$coefficients)[option], "\n")))

      nc <- ncol(x$base$coefficients)
      cmat <- matrix(NA, nrow = nc, ncol = 4)
      row_idx <- option + nrow(x$base$coefficients)
      cmat[, 1] <- format(as.numeric(x$base$coefficients[option, ]))
      cmat[, 2] <- format(x$lower[row_idx, ])
      cmat[, 3] <- format(x$larger[row_idx, ])
      cmat[, 4] <- format(x$abs[row_idx, ])
      if (x$nullhyp == "qapspp") cmat[1, 2:4] <- "*"
      colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
      rownames(cmat) <- colnames(x$base$coefficients)
      print.table(cmat)
      cat("\n\n")
    }

    if (x$nullhyp == "qapspp")
      cat("* Significance test for the intercept is undefined with qapspp.\n")

    cat("\nAIC of base model:", format(AIC(x$base$base_model)))
    cat("\nBIC of base model:", format(BIC(x$base$base_model)))
    cat("\n")
  }

  # --- confusion matrix ---
  if (!is.null(x$confusion_matrix)) {
    cat("\n")
    print(x$confusion_matrix)
  }

  cat("\n")
  invisible(x)
}
