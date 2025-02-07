#' Combine multiple estimates of the same QAP model
#'
#' The function combines multiple outputs of \code{QAPglm()} or \code{QAPcss()} estimated on the same data with the same predictors. Models are weighted by \code{reps} when combining, thus, \code{reps} does not need to be the same for all models in \code{res} and \code{res2}.
#'
#' @param res list; either a \code{list} of \code{class} \code{QAPcss}, \code{QAPregression}, or \code{QAPglm}, or a \code{list} of multiple model outputs that should all be combined.
#' @param res2 list; if \code{res} is a single model output, \code{res2} is the other model output that should be combined with the first.
#'
#' @returns Returns a \code{list} of the \code{class} of the combined model output.
#' @export

combine_qap_estimates <- function(res, res2 = NULL) {
  if (class(res) %in% c("QAPcss","QAPregression","QAPglm") &&
      !is.null(res2)) {
    res <- list(res, res2)
  }
  n_res <- length(res)

  return_res <- res[[1]]
  if (is.null(res$comp)) {
    for (i in 1:(n_res - 1)) {
      return_res$lower <- (return_res$lower  *  return_res$reps +
                             res[[i + 1]]$lower * res[[i + 1]]$reps) /
        (return_res$reps + res[[i + 1]]$reps)

      return_res$larger <- (return_res$larger  *  return_res$reps +
                              res[[i + 1]]$larger * res[[i + 1]]$reps) /
        (return_res$reps + res[[i + 1]]$reps)

      return_res$abs <- (return_res$abs  *  return_res$reps +
                           res[[i + 1]]$abs * res[[i + 1]]$reps) /
        (return_res$reps + res[[i + 1]]$reps)


      return_res$reps <- return_res$reps + res[[i + 1]]$reps

    }
  } else {
    for (i in 1:(n_res - 1)) {
      for (com in names(res$comp)) {
        return_res[[com]]$lower <- (return_res[[com]]$lower *
                                      return_res$reps +
                               res[[i + 1]][[com]]$lower *
                                 res[[i + 1]]$reps) /
          (return_res$reps + res[[i + 1]]$reps)

        return_res$larger <- (return_res[[com]]$larger *
                                return_res$reps +
                                res[[i + 1]][[com]]$larger *
                                res[[i + 1]]$reps) /
          (return_res$reps + res[[i + 1]]$reps)

        return_res$abs <- (return_res[[com]]$abs *
                             return_res$reps +
                             res[[i + 1]][[com]]$abs *
                             res[[i + 1]]$reps) /
          (return_res$reps + res[[i + 1]]$reps)


        return_res$reps <- return_res$reps + res[[i + 1]]$reps
      }
    }
  }

  return(return_res)
}
