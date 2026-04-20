# ============================================================
# Tests for QAPglm (formula-based interface)
# ============================================================

# --- Helper: generate test data ---
make_test_data <- function(n = 10, seed = 42) {
  set.seed(seed)
  y  <- matrix(rnorm(n^2), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_test_data_binary <- function(n = 10, seed = 42) {
  set.seed(seed)
  y  <- matrix(rbinom(n^2, 1, 0.3), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

# ---- basic gaussian ----
test_that("QAPglm runs with gaussian family and qapy", {
  d <- make_test_data()
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1, ncores = 1)
  expect_s3_class(fit, "QAPRegression")
  expect_true(!is.null(fit$coefficients))
  expect_true(!is.null(fit$lower))
  expect_length(fit$coefficients, 3)  # intercept + 2 predictors
})

test_that("QAPglm runs with gaussian family and qapspp", {
  d <- make_test_data()
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "gaussian", nullhyp = "qapspp",
                reps = 20, seed = 1, ncores = 1)
  expect_s3_class(fit, "QAPRegression")
  expect_true(!is.null(fit$lower))
  expect_equal(ncol(fit$lower), 3)
})

test_that("QAPglm runs with single predictor (qapspp -> qapy fallback)", {
  d <- make_test_data()
  fit <- QAPglm(y ~ x1, data = d,
                family = "gaussian", nullhyp = "qapspp",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPRegression")
  expect_equal(fit$nullhyp, "qapy")
})

# ---- binomial ----
test_that("QAPglm runs with binomial family", {
  d <- make_test_data_binary()
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "binomial", nullhyp = "qapy",
                reps = 20, seed = 1, ncores = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$confusion_matrix))
  expect_s3_class(fit$confusion_matrix, "QAPConfusionMatrix")
})

# ---- robust errors ----
test_that("QAPglm works with robust errors", {
  d <- make_test_data()
  fit <- QAPglm(y ~ x1, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1, use_robust_errors = TRUE)
  expect_s3_class(fit, "QAPRegression")
  expect_true(fit$robust_se)
})

# ---- random intercepts ----
test_that("QAPglm works with random intercept for sender", {
  d <- make_test_data()
  fit <- QAPglm(y ~ x1, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1,
                random_intercept_sender = TRUE)
  expect_s3_class(fit, "QAPRegression")
})

# ---- multiple networks ----
test_that("QAPglm works with list of matrices", {
  n <- 8
  set.seed(42)
  y1 <- matrix(rnorm(n^2), n, n);  diag(y1) <- NA
  y2 <- matrix(rnorm(n^2), n, n);  diag(y2) <- NA
  x1_1 <- matrix(rnorm(n^2), n, n); diag(x1_1) <- NA
  x1_2 <- matrix(rnorm(n^2), n, n); diag(x1_2) <- NA
  d <- list(y = list(y1, y2), x1 = list(x1_1, x1_2))

  fit <- QAPglm(y ~ x1, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPRegression")
})

# ---- print methods ----
test_that("print.QAPRegression works", {
  d <- make_test_data()
  fit <- QAPglm(y ~ x1, data = d, reps = 20, seed = 1)
  expect_output(print(fit), "OLS Network Model")
})

test_that("print.QAPGLM works", {
  d <- make_test_data_binary()
  fit <- QAPglm(y ~ x1, data = d, family = "binomial",
                reps = 20, seed = 1)
  expect_output(print(fit), "Generalized Linear Network Model")
})

# ---- comparison argument ----
test_that("QAPglm works with comparisons", {
  n <- 10
  set.seed(42)
  # Create a categorical outcome
  y <- matrix(sample(c("TP", "FP", "TN", "FN"), n^2, replace = TRUE), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  diag(y) <- NA; diag(x1) <- NA
  d <- list(y = y, x1 = x1)

  comp <- list(commission = c("FP", "TN"))
  fit <- QAPglm(y ~ x1, data = d, family = "gaussian",
                comparison = comp, nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$comp))
})

# ---- helper: count data ----
make_test_data_count <- function(n = 10, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  lambda <- exp(0.5 * x1 + 0.3 * x2)
  y  <- matrix(rpois(n^2, as.vector(lambda)), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_test_data_negbin <- function(n = 10, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  mu <- exp(0.5 * x1 + 0.3 * x2)
  theta <- 2
  y  <- matrix(rnbinom(n^2, mu = as.vector(mu), size = theta), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_test_data_zip <- function(n = 10, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  lambda <- exp(0.5 * x1 + 0.3 * x2)
  pi_z <- 0.3  # zero-inflation probability
  y <- matrix(NA, n, n)
  for (i in seq_len(n^2)) {
    if (runif(1) < pi_z) {
      y[i] <- 0
    } else {
      y[i] <- rpois(1, lambda[i])
    }
  }
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

# ---- negative binomial ----
test_that("QAPglm runs with negbin family (standard)", {
  skip_if_not_installed("MASS")
  d <- make_test_data_negbin()
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "negbin", nullhyp = "qapy",
                reps = 20, seed = 1, ncores = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$theta))
  expect_true(!is.null(fit$coefficients))
  expect_length(fit$coefficients, 3)
})

test_that("QAPglm runs with negbin family (GMM)", {
  skip_if_not_installed("gmm")
  d <- make_test_data_negbin(n = 15, seed = 99)
  # Some permutations may fail numerically; expect warnings
  fit <- suppressWarnings(
    QAPglm(y ~ x1 + x2, data = d,
           family = "negbin", estimator = "gmm",
           nullhyp = "qapy",
           reps = 20, seed = 1, ncores = 1)
  )
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$coefficients))
})

test_that("print.QAPGLM shows dispersion for negbin", {
  skip_if_not_installed("MASS")
  d <- make_test_data_negbin()
  fit <- QAPglm(y ~ x1, data = d, family = "negbin",
                reps = 20, seed = 1)
  expect_output(print(fit), "theta")
})

# ---- zero-inflated Poisson ----
test_that("QAPglm runs with zip family (standard)", {
  skip_if_not_installed("pscl")
  d <- make_test_data_zip()
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "zip", nullhyp = "qapy",
                reps = 20, seed = 1, ncores = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$zi_coefficients))
  expect_true(!is.null(fit$coefficients))
})

test_that("QAPglm runs with zip family (GMM)", {
  skip_if_not_installed("gmm")
  d <- make_test_data_zip(n = 15, seed = 99)
  fit <- suppressWarnings(
    QAPglm(y ~ x1 + x2, data = d,
           family = "zip", estimator = "gmm",
           nullhyp = "qapy",
           reps = 20, seed = 1, ncores = 1)
  )
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$coefficients))
})

test_that("print.QAPGLM shows zero-inflation for zip", {
  skip_if_not_installed("pscl")
  d <- make_test_data_zip()
  fit <- QAPglm(y ~ x1, data = d, family = "zip",
                reps = 20, seed = 1)
  expect_output(print(fit), "Zero-inflation")
})

# ---- negbin with random intercepts (requires glmmTMB) ----
test_that("QAPglm runs with negbin + random intercepts", {
  skip_if_not_installed("glmmTMB")
  d <- make_test_data_negbin()
  fit <- QAPglm(y ~ x1, data = d,
                family = "negbin", nullhyp = "qapy",
                reps = 20, seed = 1,
                random_intercept_sender = TRUE)
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$theta))
})

# ---- zip with random intercepts (requires glmmTMB) ----
test_that("QAPglm runs with zip + random intercepts", {
  skip_if_not_installed("glmmTMB")
  d <- make_test_data_zip()
  fit <- QAPglm(y ~ x1, data = d,
                family = "zip", nullhyp = "qapy",
                reps = 20, seed = 1,
                random_intercept_sender = TRUE)
  expect_s3_class(fit, "QAPGLM")
})
