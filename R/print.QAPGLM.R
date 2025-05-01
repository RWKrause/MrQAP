#' Print QAPglm() results for generalized linear regressions
#'
#' This function prints results for a \code{QAPglm()} call where \code{family} is not \code{'gaussian'}, for \code{family = 'gaussian'} use \code{print.QAPRegression()} instead.
#'
#' @docType methods
#'
#' @param x list; an \code{R} object of class \code{QAPGLM} returned by \code{QAPglm(...), where is not 'gaussian'}.
#'
#' @param print_b logical; Shall p-values derived from parameter comparisons also be returned or only those derived from T-values (which is more accurate in most cases). Default FALSE
#'
#' @param print_random logical; Shall all random intercepts be returned? Default FALSE
#' @param ... Potential other parameters to be passed to lower level functions
#'
#' @returns Prints a results table for a \code{QAPGLM} model.
#' @export
#'
#'
print.QAPGLM <- function(x,
                         ...,
                         print_b = FALSE,
                         print_random = FALSE) {
  stopifnot(inherits(x, "QAPGLM"))

  if (is.null(x$random.intercepts)) {
    cat("\nGeneralized Linear Network Model\n")
  } else {
    cat("\nGeneralized Linear Mixed Network Model fit by REML\n")
  }
  if (x$estimator == 'gmm') {
    cat("\nEstimator: Generalized Method-of-Moments.")
  }
  if (!is.null(x$groups)) {
    cat("\nPermutations were performed within groups only.")
  }

  if (x$nullhyp == 'qapy') {
    cat("\nThe outcome matrix Y was permuted",format(x$reps),'times.')
  }
  if (x$nullhyp == 'qapspp') {
    cat("\nSignificance was estimated using Dekker's")
    cat("\n  'semi-partialling plus' procedure with",
        format(x$reps),'permutations.')
  }

  if (x$diag) {
    cat("\nDiagonal values (loops) were used in the estimation")
  } else {
    cat("\nDiagonal values (loops) were ignored.")
  }

  cat('\nThe outcome was treated as', format(paste0(x$mode,'.')))



  cat("\nModel family:", format(x$family))
  cat("\n\nCoefficients:\n")
  if (print_b) {
    cmat <- matrix(NA, nrow = length(x$coefficients), ncol = 5)
    cmat[,1] <- as.vector(format(as.numeric(x$coefficients)))
    cmat[,2] <- as.vector(format(exp(as.numeric(x$coefficients))))
    cmat[,3] <- as.vector(format(x$lower[1,]))
    cmat[,4] <- as.vector(format(x$larger[1,]))
    cmat[,5] <- as.vector(format(x$abs[1,]))
    if (x$nullhyp == 'qapspp') {
      cmat[1,3:5] <- '*'
    }
    colnames(cmat) <- c("Estimate", "Exp(b)", "Pr(<=b)", "Pr(>=b)", "Pr(>=|b|)")
    rownames(cmat) <- names(x$coefficients)
    print.table(cmat)

    cat('--------------\n')
  }

  cmat <- matrix(NA, nrow = length(x$coefficients), ncol = 5)
  cmat[,1] <- as.vector(format(as.numeric(x$coefficients)))
  cmat[,2] <- as.vector(format(exp(as.numeric(x$coefficients))))
  cmat[,3] <- as.vector(format(x$lower[2,]))
  cmat[,4] <- as.vector(format(x$larger[2,]))
  cmat[,5] <- as.vector(format(x$abs[2,]))
  if (x$nullhyp == 'qapspp') {
    cmat[1,3:5] <- '*'
  }
  colnames(cmat) <- c("Estimate", "Exp(b)", "Pr(<=b)", "Pr(>=b)", "Pr(>=|b|)")
  rownames(cmat) <- names(x$coefficients)
  print.table(cmat)

  if (x$nullhyp == 'qapspp') {
    cat("\n* Significance test for the intercept is not defined with qapspp.\n")
  }

  cat('--------------\n')

  if (!is.null(x$random.intercepts) && print_random) {
    cmat <- matrix(round(x$random.intercepts,3),
                   nrow = length(x$random.intercepts),
                   ncol = 1)
    colnames(cmat) <- c("Random Intercepts:")
    print.table(cmat)
    cat('--------------\n')
  }
  if (x$estimator != 'gmm') {
    cat("\n\nAIC:", format(AIC(x$simple_fit)))
    cat("\nBIC:",   format(BIC(x$simple_fit)))
    cat("\n\n")
  }
}
