#' Print QAPglm() results for linear regression
#'
#' @param x An object of class \code{QAPRegression}.
#' @param print_b Logical; also show parameter-based p-values?
#' @param print_random Logical; show random intercepts?
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export

print.QAPRegression <- function(x, ...,
                                print_b = FALSE,
                                print_random = FALSE) {
  stopifnot(inherits(x, "QAPRegression"))

  if (is.null(x$random.intercepts)) {
    cat("\nOLS Network Model\n\n")
  } else {
    cat("\nLinear Mixed Network Model fit by REML\n\n")
  }

  if (!is.null(x$groups))
    cat("Permutations were performed within groups only.\n")

  if (x$nullhyp == "qapy")
    cat("The outcome matrix Y was permuted", format(x$reps), "times.\n")
  if (x$nullhyp == "qapspp") {
    cat("Significance was estimated using Dekker's\n")
    cat("  'semi-partialling plus' procedure with",
        format(x$reps), "permutations.\n")
  }

  if (x$diag) {
    cat("Diagonal values (loops) were used in the estimation.\n")
  } else {
    cat("Diagonal values (loops) were ignored.\n")
  }
  cat("The outcome was treated as", format(paste0(x$mode, ".")), "\n")

  if (!is.null(x$r.squared)) {
    cat("\nR-squared:    ", format(round(x$r.squared, 4)))
    cat("\nAdj R-squared:", format(round(x$adj.r.squared, 4)), "\n")
  }

  cat("\nCoefficients:\n")
  if (print_b) {
    cmat <- cbind(format(as.numeric(x$coefficients)),
                  format(x$lower[1, ]),
                  format(x$larger[1, ]),
                  format(x$abs[1, ]))
    colnames(cmat) <- c("Estimate", "Pr(<=b)", "Pr(>=b)", "Pr(>=|b|)")
    rownames(cmat) <- names(x$coefficients)
    if (x$nullhyp == "qapspp") cmat[1, 2:4] <- "*"
    print.table(cmat)
    cat("\n--------------\n")
  }

  cmat <- cbind(format(as.numeric(x$coefficients)),
                format(x$lower[2, ]),
                format(x$larger[2, ]),
                format(x$abs[2, ]))
  colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
  rownames(cmat) <- names(x$coefficients)
  if (x$nullhyp == "qapspp") cmat[1, 2:4] <- "*"
  print.table(cmat)

  if (x$nullhyp == "qapspp")
    cat("\n* Significance test for the intercept is not defined with qapspp.\n")

  cat("\n--------------\n")

  if (!is.null(x$random.intercepts) && print_random) {
    cmat <- matrix(round(x$random.intercepts, 3),
                   nrow = length(x$random.intercepts), ncol = 1)
    colnames(cmat) <- "Random Intercepts:"
    print.table(cmat)
    cat("\n--------------\n")
  }
  cat("\n")
  invisible(x)
}
