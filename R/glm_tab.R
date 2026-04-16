#' Internal auxiliary function to print QAPcss() coefficient tables
#'
#' @param x An object of class \code{QAPCSS}.
#' @param comp Character or NULL; name of comparison to print.
#'
#' @return Prints the results table.
#' @keywords internal

glm_tab <- function(x, comp) {
  if (!is.null(comp)) {
    cat("\n\nComparison between",
        x$comp[[comp]][1], "and", x$comp[[comp]][2])
    cat("\n\nCoefficients:\n")

    nc <- length(x$base[[comp]]$coefficients)
    cmat <- matrix(NA, nrow = nc, ncol = 4)
    cmat[, 1] <- format(round(as.numeric(x$base[[comp]]$coefficients), 3))
    cmat[, 2] <- format(x$lower[[comp]][2, ])
    cmat[, 3] <- format(x$larger[[comp]][2, ])
    cmat[, 4] <- format(x$abs[[comp]][2, ])
    if (x$nullhyp == "qapspp") cmat[1, 2:4] <- "*"
    colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
    rownames(cmat) <- names(x$base[[comp]]$coefficients)
    print.table(cmat)

    if (x$nullhyp == "qapspp")
      cat("\n* Significance test for the intercept is undefined with qapspp.\n")

    if (!is.null(x$base[[comp]]$base_model)) {
      cat("\nAIC of base model:", format(AIC(x$base[[comp]]$base_model)))
      cat("\nBIC of base model:", format(BIC(x$base[[comp]]$base_model)))
    }
    cat("\n")
  } else {
    cat("\n\nCoefficients:\n")

    nc <- length(x$base$coefficients)
    cmat <- matrix(NA, nrow = nc, ncol = 4)
    cmat[, 1] <- format(round(as.numeric(x$base$coefficients), 3))
    cmat[, 2] <- format(x$lower[2, ])
    cmat[, 3] <- format(x$larger[2, ])
    cmat[, 4] <- format(x$abs[2, ])
    if (x$nullhyp == "qapspp") cmat[1, 2:4] <- "*"
    colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
    rownames(cmat) <- names(x$base$coefficients)
    print.table(cmat)

    if (x$nullhyp == "qapspp")
      cat("\n* Significance test for the intercept is undefined with qapspp.\n")

    if (!is.null(x$base$base_model)) {
      cat("\nAIC of base model:", format(AIC(x$base$base_model)))
      cat("\nBIC of base model:", format(BIC(x$base$base_model)))
    }
    cat("\n")
  }
}
