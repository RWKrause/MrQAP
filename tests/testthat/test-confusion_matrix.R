# ============================================================
# Tests for the probabilistic confusion matrix
# ============================================================

test_that("probabilistic_confusion_matrix produces correct structure", {
  set.seed(42)
  actual <- c(rep(0, 50), rep(1, 50))
  probs  <- c(runif(50, 0, 0.4), runif(50, 0.6, 1))

  cm <- probabilistic_confusion_matrix(actual, probs, n_draws = 500, seed = 1)

  expect_s3_class(cm, "QAPConfusionMatrix")
  expect_true(is.matrix(cm$confusion_matrix))
  expect_equal(dim(cm$confusion_matrix), c(2, 2))
  expect_true(cm$accuracy > 0.5)
  expect_true(cm$sensitivity > 0.5)
  expect_true(cm$specificity > 0.5)
  expect_equal(cm$n, 100)
  expect_equal(cm$n_draws, 500)
})

test_that("confusion matrix sums equal N", {
  set.seed(42)
  actual <- c(0, 0, 1, 1, 1)
  probs  <- c(0.1, 0.2, 0.8, 0.9, 0.7)
  cm <- probabilistic_confusion_matrix(actual, probs, n_draws = 1000, seed = 1)

  total <- sum(cm$confusion_matrix)
  expect_equal(total, 5, tolerance = 0.01)
})

test_that("threshold comparison works", {
  actual <- c(0, 0, 1, 1)
  probs  <- c(0.1, 0.3, 0.7, 0.9)
  cm <- probabilistic_confusion_matrix(actual, probs,
                                        threshold_comparison = TRUE)
  expect_true(!is.null(cm$threshold_confusion_matrix))
  expect_equal(cm$threshold_accuracy, 1.0)
})

test_that("perfect predictions give high accuracy", {
  actual <- c(rep(0, 100), rep(1, 100))
  probs  <- c(rep(0.01, 100), rep(0.99, 100))
  cm <- probabilistic_confusion_matrix(actual, probs, n_draws = 500, seed = 1)
  expect_true(cm$accuracy > 0.95)
})

test_that("print method works without error", {
  actual <- c(0, 0, 1, 1)
  probs  <- c(0.2, 0.3, 0.7, 0.8)
  cm <- probabilistic_confusion_matrix(actual, probs, n_draws = 100)
  expect_output(print(cm), "Probabilistic Confusion Matrix")
})
