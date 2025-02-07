#' Print QAPglm() results for linear regression or linear probability model
#'
#' This function prints results for a \code{QAPglm(..., family = 'gaussian')} call.
#'
#' @docType methods
#'
#' @param x an \code{R} object of class \code{QAPRegression} returned by \code{QAPglm(..., family = 'gaussian')}.
#'
#' @param print_b Shall p-values derived from parameter comparisons also be returned or only those derived from T-values (which is more accurate in most cases). Default FALSE
#'
#' @param print_random Shall all random intercepts be returned? By default FALSE
#' @param ... Potential other parameters to be passed to lower level functions
#'
#' @returns Prints a results table for a \code{QAPRegression} model.
#' @export


print.QAPRegression <- function(x,
                                ...,
                                print_b = FALSE,
                                print_random = FALSE) {
  stopifnot(inherits(x, "QAPRegression"))
  if (is.null(x$random.intercepts)) {
    cat("\nOLS Network Model\n\n")
  } else {
    cat("\nLinear Mixed Network Model fit by REML\n\n")
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

  cat("\n\nCoefficients:\n")
  if (print_b) {
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(x$lower[1,])))
   cmat <- cbind(cmat, as.vector(format(x$larger[1,])))
   cmat <- cbind(cmat, as.vector(format(x$abs[1,])))
   colnames(cmat) <- c("Estimate", "Pr(<=b)", "Pr(>=b)", "Pr(>=|b|)")
   rownames(cmat) <- names(x$coefficients)
   if (x$nullhyp == 'qapspp') {
     cmat[1,2:4] <- '*'
   }
   print.table(cmat)

   cat('\n--------------\n')
  }
  cmat <- as.vector(format(as.numeric(x$coefficients)))
  cmat <- cbind(cmat, as.vector(format(x$lower[2,])))
  cmat <- cbind(cmat, as.vector(format(x$larger[2,])))
  cmat <- cbind(cmat, as.vector(format(x$abs[2,])))
  colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
  rownames(cmat) <- names(x$coefficients)
  if (x$nullhyp == 'qapspp') {
    cmat[1,2:4] <- '*'
  }
  print.table(cmat)

  if (x$nullhyp == 'qapspp') {
    cat("\n* Significance test for the intercept is not defined with qapspp.\n")
  }
  cat('\n--------------\n')

  if (!is.null(x$random.intercepts) && print_random) {
    cmat <- matrix(round(x$random.intercepts,3),
                   nrow = length(x$random.intercepts),
                                 ncol = 1)
    colnames(cmat) <- c("Random Intercepts:")
    if (length(unique(x$groups)) != length(x$random.intercepts)) {
      cat('Some random intercepts were not estimated.\n',
          'Group specific intercepts will not be named.\n\n')
    } else {
      rownames(cmat) <- as.character(x$groups[order(x$groups)])
    }
    print.table(cmat)
    cat('\n--------------\n')
  }
}
