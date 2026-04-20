#' Probabilistic Confusion Matrix for Binary Classification
#'
#' Creates a confusion matrix from predicted probabilities using Monte Carlo
#' sampling.  Instead of rounding predicted probabilities at 0.5, this function
#' draws from the Bernoulli distribution implied by each predicted probability,
#' producing a confusion matrix that properly accounts for prediction
#' uncertainty.
#'
#' @param actual Integer/numeric vector of observed binary outcomes (0/1).
#' @param predicted_prob Numeric vector of predicted probabilities (same length
#'   as \code{actual}).
#' @param n_draws Integer; number of Monte Carlo draws (default 1000).
#' @param seed Optional random seed for reproducibility.
#' @param threshold_comparison Logical; also compute the traditional threshold
#'   (0.5) confusion matrix for comparison (default TRUE).
#'
#' @return An object of class \code{QAPConfusionMatrix} containing:
#'   \describe{
#'     \item{confusion_matrix}{Mean confusion matrix across draws.}
#'     \item{confusion_matrix_sd}{Standard deviation of each cell.}
#'     \item{confusion_matrix_ci_lower}{2.5th percentile of each cell.}
#'     \item{confusion_matrix_ci_upper}{97.5th percentile of each cell.}
#'     \item{accuracy, sensitivity, specificity, precision, f1}{Mean metrics.}
#'     \item{accuracy_sd}{Standard deviation of accuracy.}
#'     \item{threshold_confusion_matrix}{Traditional threshold-based matrix.}
#'   }
#'
#' @export

probabilistic_confusion_matrix <- function(actual,
                                           predicted_prob,
                                           n_draws = 1000,
                                           seed = NULL,
                                           threshold_comparison = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  n <- length(actual)
  stopifnot(length(predicted_prob) == n)
  stopifnot(all(actual %in% c(0, 1)))
  stopifnot(all(predicted_prob >= 0 & predicted_prob <= 1, na.rm = TRUE))

  # Monte Carlo draws: n_obs x n_draws matrix of 0/1

  draws <- matrix(rbinom(n * n_draws, size = 1, prob = predicted_prob),
                  nrow = n, ncol = n_draws)

  actual_1 <- actual == 1
  actual_0 <- actual == 0
  n_actual_1 <- sum(actual_1)
  n_actual_0 <- sum(actual_0)

  tp <- colSums(draws[actual_1, , drop = FALSE])
  fn <- n_actual_1 - tp
  fp <- colSums(draws[actual_0, , drop = FALSE])
  tn <- n_actual_0 - fp

  make_mat <- function(v_tn, v_fp, v_fn, v_tp, FUN) {
    matrix(c(FUN(v_tn), FUN(v_fp), FUN(v_fn), FUN(v_tp)),
           nrow = 2, ncol = 2,
           dimnames = list(Actual = c("0", "1"),
                           Predicted = c("0", "1")))
  }

  # Handle division-by-zero for metrics
  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    mean(x)
  }

  result <- list(
    confusion_matrix       = make_mat(tn, fp, fn, tp, mean),
    confusion_matrix_sd    = make_mat(tn, fp, fn, tp, sd),
    confusion_matrix_ci_lower = make_mat(tn, fp, fn, tp,
                                          function(x) quantile(x, 0.025)),
    confusion_matrix_ci_upper = make_mat(tn, fp, fn, tp,
                                          function(x) quantile(x, 0.975)),
    accuracy       = mean((tp + tn) / n),
    accuracy_sd    = sd((tp + tn) / n),
    sensitivity    = safe_mean(tp / (tp + fn)),
    specificity    = safe_mean(tn / (tn + fp)),
    precision      = safe_mean(tp / (tp + fp)),
    f1             = safe_mean(2 * tp / (2 * tp + fp + fn)),
    n_draws        = n_draws,
    n              = n
  )

  if (threshold_comparison) {
    pred_class <- as.integer(predicted_prob >= 0.5)
    result$threshold_confusion_matrix <- table(
      Actual    = factor(actual, levels = c(0, 1)),
      Predicted = factor(pred_class, levels = c(0, 1))
    )
    result$threshold_accuracy <- mean(pred_class == actual)
  }

  class(result) <- "QAPConfusionMatrix"
  return(result)
}


#' @export
print.QAPConfusionMatrix <- function(x, ...) {
  cat("\nProbabilistic Confusion Matrix\n")
  cat("(based on", x$n_draws,
      "Monte Carlo draws from predicted probabilities)\n\n")

  cat("Mean Confusion Matrix:\n")
  print(round(x$confusion_matrix, 1))

  cat("\n95% Credible Intervals:\n")
  ci <- paste0("[", round(x$confusion_matrix_ci_lower, 1), ", ",
               round(x$confusion_matrix_ci_upper, 1), "]")
  ci_mat <- matrix(ci, nrow = 2, ncol = 2,
                   dimnames = dimnames(x$confusion_matrix))
  print(noquote(ci_mat))

  cat("\nDerived Metrics (mean across draws):\n")
  cat("  Accuracy:    ", round(x$accuracy, 4),
      " (SD:", round(x$accuracy_sd, 4), ")\n")
  cat("  Sensitivity: ", round(x$sensitivity, 4), "\n")
  cat("  Specificity: ", round(x$specificity, 4), "\n")
  cat("  Precision:   ", round(x$precision, 4), "\n")
  cat("  F1 Score:    ", round(x$f1, 4), "\n")

  if (!is.null(x$threshold_confusion_matrix)) {
    cat("\nFor comparison - Threshold (0.5) Confusion Matrix:\n")
    print(x$threshold_confusion_matrix)
    cat("  Accuracy:", round(x$threshold_accuracy, 4), "\n")
  }
  cat("\n")
}
