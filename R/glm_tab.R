#' Internal auxiliary function to help plot QAPcss() output
#'
#' @param x list; an \code{R} object of class \code{QAPCSS} returned by \code{QAPcss(...)}.
#' @param print_b logical; Shall p-values derived from parameter comparisons also be returned or only those derived from T-values (which is more accurate in most cases). Default FALSE
#' @param print_random logical; Shall all random intercepts be returned? Default FALSE
#' @param comp list; see \code{QAPcss(...)}
#' @param nullhyp character or numeric; see \code{QAPcss(...)}
#'
#' @returns Prints the results table.

glm_tab <- function(x,
                    print_b = print_b,
                    print_random = print_random,
                    comp,
                    nullhyp) {

  if (!is.null(comp)) {
    cat("\n\nComparison between ",
        comp[1], ' and ', comp[2])
  }

  cat("\n\nCoefficients:\n")
  if (print_b) {
    cmat <- matrix(NA, nrow = length(x$coefficients), ncol = 4)
    cmat[,1] <- as.vector(format(round(as.numeric(x$coefficients),3)))
    cmat[,2] <- as.vector(format(x$lower[1,]))
    cmat[,3] <- as.vector(format(x$larger[1,]))
    cmat[,4] <- as.vector(format(x$abs[1,]))
    if (x$nullhyp == 'qapspp') {
      cmat[1,2:4] <- '*'
    }
    colnames(cmat) <- c("Estimate", "Pr(<b)", "Pr(>b)", "Pr(>|b|)")
    rownames(cmat) <- names(x$coefficients)
    print.table(cmat)

    cat('--------------\n')
  }

  cmat <- matrix(NA, nrow = length(x$coefficients), ncol = 4)
  cmat[,1] <- as.vector(format(round(as.numeric(x$coefficients),3)))
  cmat[,2] <- as.vector(format(x$lower[2,]))
  cmat[,3] <- as.vector(format(x$larger[2,]))
  cmat[,4] <- as.vector(format(x$abs[2,]))
  if (nullhyp == 'qapspp') {
    cmat[1,2:4] <- '*'
  }
  colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
  rownames(cmat) <- names(x$coefficients)
  print.table(cmat)

  if (nullhyp == 'qapspp') {
    cat("\n* Significance test for the intercept is undefined with qapspp.\n")
  }

  cat('--------------\n')

  if (!is.null(x$random.intercepts) && print_random) {
    cat('\nRandom Intercepts\n')

    cmat <- round(x$random.intercepts,3)

    print.table(cmat)
    cat('--------------\n')
  }

  cat("\nAIC:", format(AIC(x$base_model)))
  cat("\nBIC:", format(BIC(x$base_model)))
  cat("\n")
}

