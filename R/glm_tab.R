#' Internal auxiliary function to help plot QAPcss() output
#'
#' @param x list; an \code{R} object of class \code{QAPCSS} returned by \code{QAPcss(...)}.
#' @param comp list; see \code{QAPcss(...)}
#'
#' @returns Prints the results table.

glm_tab <- function(x,
                    comp) {

  if (!is.null(comp)) {
    cat("\n\nComparison between ",
        x$comp[[comp]][1], ' and ', x$comp[[comp]][2])
    cat("\n\nCoefficients:\n")

    cmat <- matrix(NA, nrow = length(x$base[[comp]]$coefficients), ncol = 4)
    cmat[,1] <- as.vector(format(round(as.numeric(x$base[[comp]]$coefficients),3)))
    cmat[,2] <- as.vector(format(x$lower[[comp]][2,]))
    cmat[,3] <- as.vector(format(x$larger[[comp]][2,]))
    cmat[,4] <- as.vector(format(x$abs[[comp]][2,]))
    if (x$nullhyp == 'qapspp') {
      cmat[1,2:4] <- '*'
    }
    colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
    rownames(cmat) <- names(x$base[[comp]]$coefficients)
    print.table(cmat)

    if (x$nullhyp == 'qapspp') {
      cat("\n* Significance test for the intercept is undefined with qapspp.\n")
    }


    cat("\nAIC of base model:", format(AIC(x$base[[comp]]$base_model)))
    cat("\nBIC of base model:", format(BIC(x$base[[comp]]$base_model)))
    cat("\n")
  } else {
    cat("\n\nCoefficients:\n")

    cmat <- matrix(NA, nrow = length(x$base$coefficients), ncol = 4)
    cmat[,1] <- as.vector(format(round(as.numeric(x$base$coefficients),3)))
    cmat[,2] <- as.vector(format(x$lower[2,]))
    cmat[,3] <- as.vector(format(x$larger[2,]))
    cmat[,4] <- as.vector(format(x$abs[2,]))
    if (x$nullhyp == 'qapspp') {
      cmat[1,2:4] <- '*'
    }
    colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
    rownames(cmat) <- names(x$base$coefficients)
    print.table(cmat)

    if (x$nullhyp == 'qapspp') {
      cat("\n* Significance test for the intercept is undefined with qapspp.\n")
    }

    cat("\nAIC of base model:", format(AIC(x$base$base_model)))
    cat("\nBIC of base model:", format(BIC(x$base$base_model)))
    cat("\n")
  }
}

