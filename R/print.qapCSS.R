#' Print QAPcss() results
#'
#' This function prints results for \code{QAPcsss()}.
#'
#' @docType methods
#'
#' @param x list; an \code{R} object of class \code{QAPCSS} returned by \code{QAPcss()}.
#'
#' @param print_b logical; Shall p-values derived from parameter comparisons also be returned or only those derived from T-values (which is more accurate in most cases). Default FALSE
#'
#' @param print_random logical; Shall all random intercepts be returned? Default FALSE
#' @param ... Potential other parameters to be passed to lower level functions
#'
#' @returns Prints a results table for a \code{QAPcss()} model.
#' @export
#'

print.QAPCSS <- function(x,
                          print_b = FALSE,
                          print_random = FALSE,
                          ...) {

  source('glm_tab.R')
  if (x$family != 'multinom') {
    if (!any(x$random)) {
      cat("\nGeneralized Linear Network Model for CSS\n\n")
    } else {
      cat("\nGeneralized Linear Mixed Network Model for CSS fit by REML\n\n")
    }

    if (!is.null(x$groups)) {
      cat("\nPermutations were performed within groups only.")
    }

    if (x$nullhyp == 'qapy') {
      cat("\nThe outcome array Y was permuted",format(x$reps),'times.')
    }
    if (x$nullhyp == 'qapspp') {
      cat("\nSignificance was estimated using Dekker's")
      cat("\n  'semi-partialling plus' procedure with",
          format(x$reps),'permutations.')
    }

    if (x$robust_se) {
      cat("\nT-values are based on robust standard errors.")
    }

    if (x$diag) {
      cat("\nDiagonal values (loops) were used in the estimation.",
          "\n  Results may be biased because of that.")
    } else {
      cat("\nDiagonal values (loops) were ignored.")
    }

    cat('\nThe outcome was treated as', format(paste0(x$mode,'.')))

    if (is.null(x$comp)) {
      glm_tab(x,
              print_b = print_b,
              print_random = print_random,
              nullhyp = x$nullhyp,
              groups = x$groups,
              comp = x$comp)
    } else {
      nn <- names(x)[!(names(x) %in% c("nullhyp",
                                       "family",
                                       "groups",
                                       "diag",
                                       "mode",
                                       "reps",
                                       "comp",
                                       "random",
                                       "robust_se"))]
      for (mod in 1:length(x$comp)) {
        glm_tab(x[[nn[mod]]],
                print_b = print_b,
                print_random = print_random,
                nullhyp = x$nullhyp,
                groups = x$groups,
                comp = x$comp[[mod]])
      }
    }

  } else {
    cat("\nMultinominal Choice Network Model for CSS\n\n")


    if (!is.null(x$groups)) {
      cat("\nPermutations were performed within groups only.")
    }

    if (x$nullhyp == 'qapy') {
      cat("\nThe outcome array Y was permuted",format(x$reps),'times.')
    }
    if (x$nullhyp == 'qapspp') {
      cat("\nSignificance was estimated using Dekker's")
      cat("\n  'semi-partialling plus' procedure with",
          format(x$reps),'permutations.')
    }

    if (x$diag) {
      cat("\nDiagonal values (loops) were used in the estimation.",
          "\n  Results may be biased because of that.")
    } else {
      cat("\nDiagonal values (loops) were ignored.")
    }

    cat('\nThe outcome was treated as', format(paste0(x$mode,'.')))

    cat('\nThe reference group was', format(paste0(x$reference,'.')))


    cat("\n\nCoefficients:\n\n")

    for (option in 1:nrow(x$coefficients)) {
      cat(format(paste0('-- ',rownames(x$coefficients)[option],'\n')))

      cmat <- matrix(NA, nrow = ncol(x$coefficients), ncol = 4)
      cmat[,1] <- as.vector(format(as.numeric(x$coefficients[option,])))
      cmat[,2] <- as.vector(format(x$lower[option + nrow(x$coefficients),]))
      cmat[,3] <- as.vector(format(x$larger[option + nrow(x$coefficients),]))
      cmat[,4] <- as.vector(format(x$abs[option + nrow(x$coefficients),]))
      if (x$nullhyp == 'qapspp') {
        cmat[1,2:4] <- '*'
      }
      colnames(cmat) <- c("Estimate", "Pr(<=t)", "Pr(>=t)", "Pr(>=|t|)")
      rownames(cmat) <- colnames(x$coefficients)
      print.table(cmat)
      cat("\n\n")

    }


    if (x$nullhyp == 'qapspp') {
      cat("\n* Significance test for the intercept is undefined with qapspp.\n")
    }

    cat("\nAIC:", format(AIC(x$simple_fit)))
    cat("\nBIC:", format(BIC(x$simple_fit)))
    cat("\n")
  }

}
