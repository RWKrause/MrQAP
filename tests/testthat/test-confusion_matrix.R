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

test_that("threshold metrics are computed at the 0.5 cut-off", {
  actual <- c(0, 0, 0, 1, 1, 1)
  probs  <- c(0.1, 0.6, 0.2, 0.7, 0.4, 0.9)
  cm <- probabilistic_confusion_matrix(actual, probs, n_draws = 10, seed = 1)
  # tp = 2, fp = 1, fn = 1, tn = 2
  expect_equal(cm$threshold_accuracy, 4 / 6)
  expect_equal(cm$threshold_sensitivity, 2 / 3)
  expect_equal(cm$threshold_specificity, 2 / 3)
  expect_equal(cm$threshold_precision, 2 / 3)

  none <- probabilistic_confusion_matrix(actual, rep(0.2, 6), n_draws = 10,
                                         seed = 1)
  expect_true(is.na(none$threshold_precision))
  expect_equal(none$threshold_sensitivity, 0)
})

# --- QAP() integration ------------------------------------------------------

cm_data <- function() {
  set.seed(3)
  n  <- 15
  x1 <- matrix(rnorm(n^2), n)
  x2 <- matrix(rnorm(n^2), n)
  y  <- matrix(rbinom(n^2, 1, plogis(-0.3 + 1.2 * x1 + 0.2 * x2)), n)
  list(y = y, x1 = x1, x2 = x2)
}

test_that("binomial fits carry a metric_change table, full minus reduced", {
  d   <- cm_data()
  fit <- QAP(y ~ x1 + x2, data = d, family = "binomial", reps = 5, seed = 7)
  mc  <- fit$metric_change
  expect_s3_class(mc, "data.frame")
  expect_equal(ncol(mc), 10)
  expect_equal(mc$term, c("x1", "x2"))

  # The x2 row is the full model's metrics minus those of y ~ x1.
  red <- glm(y ~ x1, family = binomial,
             data = data.frame(y = d$y[row(d$y) != col(d$y)],
                               x1 = d$x1[row(d$x1) != col(d$x1)]))
  red_cm <- probabilistic_confusion_matrix(red$y, fitted(red), seed = 7)
  expect_equal(unlist(mc[2, 2:9]),
               qap_class_metrics(fit$confusion_matrix) -
                 qap_class_metrics(red_cm))

  out <- capture.output(summary(fit))
  expect_true(any(grepl("Change in classification", out)))
  expect_true(any(grepl("Sens(draw)", out, fixed = TRUE)))
})

test_that("draws use 350415 when no seed is given", {
  d <- cm_data()
  a <- QAP(y ~ x1 + x2, data = d, family = "binomial", reps = 5)
  b <- QAP(y ~ x1 + x2, data = d, family = "binomial", reps = 5)
  expect_identical(a$confusion_matrix$accuracy, b$confusion_matrix$accuracy)
  expect_identical(a$metric_change, b$metric_change)
  o <- qap_binary_outcome(list(base_model = a$simple_fit))
  ref <- probabilistic_confusion_matrix(o$actual, o$prob, seed = 350415)
  expect_identical(a$confusion_matrix$accuracy, ref$accuracy)
})

test_that("less_mem keeps the classification metrics", {
  d <- cm_data()
  a <- QAP(y ~ x1 + x2, data = d, family = "binomial", reps = 5, seed = 7)
  b <- QAP(y ~ x1 + x2, data = d, family = "binomial", reps = 5, seed = 7,
           less_mem = TRUE)
  expect_equal(b$metric_change, a$metric_change)
  expect_s3_class(b$confusion_matrix, "QAPConfusionMatrix")
})

test_that("comparison models get one matrix and one table per comparison", {
  d  <- cm_data()
  set.seed(4)
  yc <- matrix(sample(c("a", "b", "c"), length(d$y), TRUE), nrow(d$y))
  yc[d$x1 > 0.5] <- "c"
  fit <- QAP(y ~ x1 + x2, data = list(y = yc, x1 = d$x1, x2 = d$x2),
             family = "binomial", reps = 5, seed = 7,
             comparison = list(ab = c("a", "b"), ac = c("a", "c")))
  expect_named(fit$confusion_matrix, c("ab", "ac"))
  expect_s3_class(fit$confusion_matrix$ac, "QAPConfusionMatrix")
  expect_named(fit$metric_change, c("ab", "ac"))
  expect_equal(fit$metric_change$ac$term, c("x1", "x2"))
  # n of each matrix is the size of that comparison's subset.
  expect_equal(fit$confusion_matrix$ab$n,
               sum(yc[row(yc) != col(yc)] %in% c("a", "b")))
  expect_output(print(fit), "Comparison: ac")
  expect_output(print(summary(fit)), "Change in classification")
})

test_that("terms are dropped from fixest and mixed formulas", {
  expect_equal(qap_fixed_terms(y ~ x1 + x2 | fe, TRUE), c("x1", "x2"))
  expect_equal(deparse(qap_drop_term(y ~ x1 + x2 | fe, "x1", TRUE)),
               "y ~ x2 | fe")
  expect_equal(deparse(qap_drop_term(y ~ x1 | fe, "x1", TRUE)), "y ~ 1 | fe")
  expect_equal(qap_fixed_terms(y ~ x1 * x2 + (1 | sv), FALSE),
               c("x1", "x2", "x1:x2"))
  expect_equal(deparse(qap_drop_term(y ~ x1 + x2 + (1 | sv), "x1", FALSE)),
               "y ~ x2 + (1 | sv)")
})

test_that("qap_is_binary counts distinct non-NA values of any type", {
  expect_true(qap_is_binary(c(0, 1, 1, NA)))
  expect_true(qap_is_binary(c("a", "b", "a")))
  expect_true(qap_is_binary(c(1, 2)))
  expect_false(qap_is_binary(c(0, 1, 2)))
  expect_false(qap_is_binary(c(1, 1, NA)))
  expect_false(qap_is_binary(numeric(0)))
})

# --- linear probability models ------------------------------------------------

lpm_data <- function() {
  set.seed(3)
  n  <- 15
  x1 <- matrix(rnorm(n^2), n)
  x2 <- matrix(rnorm(n^2), n)
  y  <- matrix(rbinom(n^2, 1, plogis(-0.5 + 2 * x1 + 0.2 * x2)), n)
  list(y = y, x1 = x1, x2 = x2)
}
off <- function(m) m[row(m) != col(m)]

test_that("a gaussian fit on a 0/1 outcome is treated as an LPM", {
  .qap_state$lpm_warned <- FALSE
  d <- lpm_data()
  expect_warning(
    fit <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 7),
    "Horrace & Oaxaca")
  cm <- fit$confusion_matrix
  expect_s3_class(cm, "QAPConfusionMatrix")
  expect_equal(fit$metric_change$term, c("x1", "x2"))

  fv <- fitted(lm(off(d$y) ~ off(d$x1) + off(d$x2)))
  expect_equal(cm$n_out, sum(fv < 0 | fv > 1))
  expect_gt(cm$n_out, 0)
  # Clipping cannot move a value across 0.5.
  expect_equal(cm$threshold_accuracy, mean((fv >= 0.5) == off(d$y)))

  # Once per session: a second fit with values out of range is quiet.
  expect_no_warning(QAP(y ~ x1 + x2, data = d, reps = 5, seed = 7))

  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("lay outside [0, 1]", out, fixed = TRUE)))
  expect_true(any(grepl("Horrace & Oaxaca", out)))
  expect_output(print(fit), "Probabilistic Confusion Matrix")
})

test_that("a two-valued LPM outcome coded 1/2 gives the 0/1 metrics", {
  d <- lpm_data()
  d2 <- d
  d2$y <- d$y + 1
  a <- suppressWarnings(QAP(y ~ x1 + x2, data = d,  reps = 5, seed = 7))
  b <- suppressWarnings(QAP(y ~ x1 + x2, data = d2, reps = 5, seed = 7))
  expect_equal(b$metric_change, a$metric_change)
  expect_equal(b$confusion_matrix$n_out, a$confusion_matrix$n_out)
})

test_that("mixed LPMs and gaussian comparisons get confusion matrices", {
  d <- lpm_data()
  r <- suppressWarnings(QAP(y ~ x1 + x2, data = d, reps = 5, seed = 7,
                            random_intercept_sender = TRUE))
  expect_s3_class(r$confusion_matrix, "QAPConfusionMatrix")
  expect_false(is.null(r$confusion_matrix$n_out))

  set.seed(4)
  yc <- matrix(sample(c("a", "b", "c"), length(d$y), TRUE), nrow(d$y))
  yc[d$x1 > 0.5] <- "c"
  fc <- suppressWarnings(QAP(y ~ x1 + x2,
                             data = list(y = yc, x1 = d$x1, x2 = d$x2),
                             reps = 5, seed = 7,
                             comparison = list(ab = c("a", "b"),
                                               ac = c("a", "c"))))
  expect_named(fc$confusion_matrix, c("ab", "ac"))
  expect_false(is.null(fc$confusion_matrix$ac$n_out))
  expect_named(fc$metric_change, c("ab", "ac"))
})

test_that("non-binary outcomes get no confusion matrix", {
  d <- lpm_data()
  g <- QAP(y ~ x1, data = list(y = d$x2, x1 = d$x1), reps = 5)
  expect_null(g$confusion_matrix)
  expect_null(g$metric_change)
  set.seed(5)
  yp <- matrix(rpois(length(d$y), 2), nrow(d$y))
  p <- QAP(y ~ x1, data = list(y = yp, x1 = d$x1), family = "poisson",
           reps = 5)
  expect_null(p$confusion_matrix)
  expect_false(any(grepl("Change in classification",
                         capture.output(print(summary(p))))))
})

# --- GMM and test_parameters ---------------------------------------------------

test_that("binomial GMM fits get probabilities rebuilt from the coefficients", {
  d <- lpm_data()
  fit <- QAP(y ~ x1 + x2, data = d, family = "binomial", estimator = "gmm",
             reps = 5, seed = 7)
  b <- fit$coefficients
  p <- plogis(b[1] + b[2] * off(d$x1) + b[3] * off(d$x2))
  ref <- probabilistic_confusion_matrix(off(d$y), p, seed = 7)
  expect_equal(fit$confusion_matrix$accuracy, ref$accuracy)
  expect_equal(fit$confusion_matrix$threshold_precision,
               ref$threshold_precision)
  expect_equal(fit$metric_change$term, c("x1", "x2"))
  expect_false(anyNA(fit$metric_change$accuracy_draw))
})

test_that("metric_change covers only terms involving a tested predictor", {
  d  <- lpm_data()
  d$x3 <- matrix(rnorm(length(d$y)), nrow(d$y))
  fit <- QAP(y ~ x1 * x2 + x3, data = d, family = "binomial", reps = 5,
             seed = 7, test_parameters = "x2")
  expect_equal(fit$metric_change$term, c("x2", "x1:x2"))
  all_t <- QAP(y ~ x1 * x2 + x3, data = d, family = "binomial", reps = 5,
               seed = 7)
  expect_equal(all_t$metric_change$term, c("x1", "x2", "x3", "x1:x2"))
  expect_equal(fit$metric_change,
               all_t$metric_change[all_t$metric_change$term %in%
                                     c("x2", "x1:x2"), ],
               ignore_attr = TRUE)
})

test_that("one-class reduced models are starred, with a footnote", {
  d   <- cm_data()
  # Without x1 only the intercept is left: one probability, one class.
  fit <- QAP(y ~ x1, data = d, family = "binomial", reps = 5, seed = 7)
  expect_true(fit$metric_change$one_class)
  tab <- capture.output(qap_print_metric_change(fit$metric_change))
  row <- grep("^x1 ", tab, value = TRUE)
  expect_length(gregexpr("*", row, fixed = TRUE)[[1]], 4)
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("same side of 0.5", out, fixed = TRUE)))

  # A reduced model predicting both classes gets neither.
  two <- QAP(y ~ x1 + x2, data = d, family = "binomial", reps = 5, seed = 7)
  expect_false(two$metric_change$one_class[two$metric_change$term == "x2"])
  mc <- two$metric_change[two$metric_change$term == "x2", ]
  out2 <- capture.output(qap_print_metric_change(mc))
  expect_false(any(grepl("*", out2, fixed = TRUE)))
})
