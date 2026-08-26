# ============================================================
# Tests for QAPcss (formula-based interface)
# ============================================================

make_css_test_data <- function(n = 5, seed = 42) {
  set.seed(seed)
  y  <- array(rnorm(n^3), dim = c(n, n, n))
  x1 <- array(rnorm(n^3), dim = c(n, n, n))
  x2 <- array(rnorm(n^3), dim = c(n, n, n))
  for (k in 1:n) {
    diag(y[, , k])  <- NA
    diag(x1[, , k]) <- NA
    diag(x2[, , k]) <- NA
  }
  list(y = y, x1 = x1, x2 = x2)
}

make_css_test_binary <- function(n = 5, seed = 42) {
  set.seed(seed)
  y  <- array(rbinom(n^3, 1, 0.3), dim = c(n, n, n))
  x1 <- array(rnorm(n^3), dim = c(n, n, n))
  for (k in 1:n) {
    diag(y[, , k])  <- NA
    diag(x1[, , k]) <- NA
  }
  list(y = y, x1 = x1)
}

test_that("QAPcss runs with gaussian family and qapy", {
  d <- make_css_test_data()
  fit <- QAPcss(y ~ x1 + x2, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1, ncores = 1)
  expect_s3_class(fit, "QAPCSS")
  expect_true(!is.null(fit$base$coefficients))
})

test_that("QAPcss runs with gaussian and qapspp", {
  d <- make_css_test_data()
  fit <- QAPcss(y ~ x1 + x2, data = d,
                family = "gaussian", nullhyp = "qapspp",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPCSS")
})

test_that("QAPcss runs with binomial family", {
  d <- make_css_test_binary()
  fit <- QAPcss(y ~ x1, data = d,
                family = "binomial", nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPCSS")
  expect_true(!is.null(fit$confusion_matrix))
})

test_that("QAPcss works with random intercept for perceiver", {
  d <- make_css_test_data()
  fit <- QAPcss(y ~ x1, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1,
                random_intercept_perceiver = TRUE)
  expect_s3_class(fit, "QAPCSS")
})

test_that("QAPcss print method works", {
  d <- make_css_test_data()
  fit <- QAPcss(y ~ x1, data = d, reps = 20, seed = 1)
  expect_output(print(fit), "CSS")
})

test_that("QAPcss works with comparison argument", {
  n <- 5
  set.seed(42)
  y  <- array(sample(c("TP", "FP", "TN", "FN"), n^3, replace = TRUE),
              dim = c(n, n, n))
  x1 <- array(rnorm(n^3), dim = c(n, n, n))
  for (k in 1:n) {
    diag(y[, , k])  <- NA
    diag(x1[, , k]) <- NA
  }
  d <- list(y = y, x1 = x1)

  comp <- list(commission = c("FP", "TN"))
  fit <- QAPcss(y ~ x1, data = d, family = "gaussian",
                comparison = comp, nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPCSS")
  expect_true(!is.null(fit$comp))
})
