# ============================================================
# Tests for GMM estimation (moment functions, residuals, QAPglm GMM)
# ============================================================

# --- Helper: generate test data ---
make_gmm_binary <- function(n = 15, seed = 42) {
  set.seed(seed)
  y  <- matrix(rbinom(n^2, 1, 0.3), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_gmm_count <- function(n = 15, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  lambda <- exp(0.3 * x1 + 0.2 * x2)
  y  <- matrix(rpois(n^2, as.vector(lambda)), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_gmm_negbin <- function(n = 15, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  mu <- exp(0.3 * x1 + 0.2 * x2)
  y  <- matrix(rnbinom(n^2, mu = as.vector(mu), size = 2), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_gmm_zip <- function(n = 15, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  lambda <- exp(0.3 * x1 + 0.2 * x2)
  pi_z <- 0.3
  y <- matrix(NA, n, n)
  for (i in seq_len(n^2)) {
    y[i] <- if (runif(1) < pi_z) 0 else rpois(1, lambda[i])
  }
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}


# ============================================================
# Unit tests for moment functions
# ============================================================

test_that("logit_moments returns n x p matrix", {
  n <- 50; p <- 3
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Y <- rbinom(n, 1, 0.4)
  theta <- rnorm(p)
  g <- logit_moments(theta, list(y = Y, x = X))
  expect_equal(nrow(g), n)
  expect_equal(ncol(g), p)
})

test_that("logit_moments returns zero-mean at true parameters", {
  set.seed(1)
  n <- 5000; p <- 2
  X <- cbind(1, rnorm(n))
  beta_true <- c(-0.5, 1.0)
  prob <- 1 / (1 + exp(-X %*% beta_true))
  Y <- rbinom(n, 1, prob)
  g <- logit_moments(beta_true, list(y = Y, x = X))
  # At the true parameters, mean of moment conditions should be ~0
  expect_true(all(abs(colMeans(g)) < 0.1))
})

test_that("poisson_moments returns n x p matrix", {
  n <- 50; p <- 3
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Y <- rpois(n, 2)
  theta <- rnorm(p)
  g <- poisson_moments(theta, list(y = Y, x = X))
  expect_equal(nrow(g), n)
  expect_equal(ncol(g), p)
})

test_that("negbin_moments returns n x (p+1) matrix", {
  n <- 50; p <- 3
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Y <- rnbinom(n, mu = 3, size = 2)
  theta <- rnorm(p + 1)  # p regression + log(alpha)
  g <- negbin_moments(theta, list(y = Y, x = X))
  expect_equal(nrow(g), n)
  expect_equal(ncol(g), p + 1)
})

test_that("zip_moments returns n x (p+1) matrix", {
  n <- 50; p <- 3
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Y <- rpois(n, 2); Y[sample(n, 15)] <- 0
  theta <- rnorm(p + 1)
  g <- zip_moments(theta, list(y = Y, x = X))
  expect_equal(nrow(g), n)
  expect_equal(ncol(g), p + 1)
})


# ============================================================
# Unit tests for residual functions
# ============================================================

test_that("logit_resid returns correct length", {
  n <- 30; p <- 3
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Y <- rbinom(n, 1, 0.4)
  gmmo <- list(dat = list(y = Y, x = X),
               coefficients = rnorm(p))
  r <- logit_resid(gmmo)
  expect_length(r, n)
  expect_true(all(abs(r) <= 1))  # bounded for logistic
})

test_that("poisson_resid returns correct length", {
  n <- 30; p <- 3
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Y <- rpois(n, 2)
  gmmo <- list(dat = list(y = Y, x = X),
               coefficients = rep(0, p))
  r <- poisson_resid(gmmo)
  expect_length(r, n)
})

test_that("negbin_resid returns correct length", {
  n <- 30; p <- 3
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Y <- rnbinom(n, mu = 3, size = 2)
  gmmo <- list(dat = list(y = Y, x = X),
               coefficients = c(rep(0, p), 0))  # p + log(alpha)
  r <- negbin_resid(gmmo)
  expect_length(r, n)
})

test_that("zip_resid returns correct length", {
  n <- 30; p <- 3
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Y <- rpois(n, 2); Y[sample(n, 10)] <- 0
  gmmo <- list(dat = list(y = Y, x = X),
               coefficients = c(rep(0, p), 0))  # p + logit(pi)
  r <- zip_resid(gmmo)
  expect_length(r, n)
})


# ============================================================
# Integration tests: QAPglm with GMM estimator
# ============================================================

test_that("QAPglm GMM works for binomial", {
  skip_if_not_installed("gmm")
  d <- make_gmm_binary()
  fit <- suppressWarnings(
    QAPglm(y ~ x1 + x2, data = d,
           family = "binomial", estimator = "gmm",
           nullhyp = "qapy", reps = 20, seed = 1, ncores = 1)
  )
  expect_s3_class(fit, "QAPGLM")
  expect_equal(fit$estimator, "gmm")
  expect_length(fit$coefficients, 3)
  expect_true(all(!is.na(fit$coefficients)))
})

test_that("QAPglm GMM works for poisson", {
  skip_if_not_installed("gmm")
  d <- make_gmm_count()
  fit <- suppressWarnings(
    QAPglm(y ~ x1 + x2, data = d,
           family = "poisson", estimator = "gmm",
           nullhyp = "qapy", reps = 20, seed = 1, ncores = 1)
  )
  expect_s3_class(fit, "QAPGLM")
  expect_equal(fit$estimator, "gmm")
  expect_length(fit$coefficients, 3)
})

test_that("QAPglm GMM works for negbin", {
  skip_if_not_installed("gmm")
  d <- make_gmm_negbin(n = 15, seed = 99)
  fit <- suppressWarnings(
    QAPglm(y ~ x1 + x2, data = d,
           family = "negbin", estimator = "gmm",
           nullhyp = "qapy", reps = 20, seed = 1, ncores = 1)
  )
  expect_s3_class(fit, "QAPGLM")
  expect_equal(fit$estimator, "gmm")
  expect_length(fit$coefficients, 3)
})

test_that("QAPglm GMM works for zip", {
  skip_if_not_installed("gmm")
  d <- make_gmm_zip(n = 15, seed = 99)
  fit <- suppressWarnings(
    QAPglm(y ~ x1 + x2, data = d,
           family = "zip", estimator = "gmm",
           nullhyp = "qapy", reps = 20, seed = 1, ncores = 1)
  )
  expect_s3_class(fit, "QAPGLM")
  expect_equal(fit$estimator, "gmm")
  expect_length(fit$coefficients, 3)
})

test_that("QAPglm GMM with robust errors", {
  skip_if_not_installed("gmm")
  d <- make_gmm_binary()
  fit <- suppressWarnings(
    QAPglm(y ~ x1 + x2, data = d,
           family = "binomial", estimator = "gmm",
           use_robust_errors = TRUE,
           nullhyp = "qapy", reps = 20, seed = 1, ncores = 1)
  )
  expect_s3_class(fit, "QAPGLM")
  expect_true(fit$robust_se)
  expect_true(all(!is.na(fit$t)))
})

test_that("QAPglm GMM rejects unsupported family", {
  skip_if_not_installed("gmm")
  d <- list(y = matrix(rnorm(100), 10, 10),
            x1 = matrix(rnorm(100), 10, 10))
  diag(d$y) <- diag(d$x1) <- NA
  expect_error(
    QAPglm(y ~ x1, data = d,
           family = "gaussian", estimator = "gmm",
           reps = 5, seed = 1),
    "GMM estimator"
  )
})

test_that("QAPglm GMM print shows estimator", {
  skip_if_not_installed("gmm")
  d <- make_gmm_binary()
  fit <- suppressWarnings(
    QAPglm(y ~ x1, data = d,
           family = "binomial", estimator = "gmm",
           reps = 20, seed = 1)
  )
  expect_output(print(fit), "Method-of-Moments")
})
