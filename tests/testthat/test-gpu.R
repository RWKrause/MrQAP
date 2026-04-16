# ============================================================
# Tests for GPU batch OLS (QAPglm and QAPcss)
# ============================================================

# All tests skip if torch is not installed

# --- Helper: generate test data ---
make_gpu_data <- function(n = 10, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- 0.5 * x1 + 0.3 * x2 + matrix(rnorm(n^2), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_gpu_css_data <- function(n = 6, seed = 42) {
  set.seed(seed)
  x1 <- array(rnorm(n^3), dim = c(n, n, n))
  x2 <- array(rnorm(n^3), dim = c(n, n, n))
  y  <- 0.5 * x1 + 0.3 * x2 + array(rnorm(n^3), dim = c(n, n, n))
  for (i in 1:n) {
    diag(y[,,i]) <- diag(x1[,,i]) <- diag(x2[,,i]) <- NA
  }
  list(y = y, x1 = x1, x2 = x2)
}

# ---- QAPglm GPU qapy ----
test_that("QAPglm GPU works for gaussian qapy", {
  skip_if_not_installed("torch")
  d <- make_gpu_data()
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1, use_gpu = TRUE)
  expect_s3_class(fit, "QAPRegression")
  expect_length(fit$coefficients, 3)
  expect_true(!is.null(fit$lower))
  expect_true(!is.null(fit$larger))
  expect_true(!is.null(fit$abs))
})

# ---- QAPglm GPU qapspp ----
test_that("QAPglm GPU works for gaussian qapspp", {
  skip_if_not_installed("torch")
  d <- make_gpu_data()
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "gaussian", nullhyp = "qapspp",
                reps = 20, seed = 1, use_gpu = TRUE)
  expect_s3_class(fit, "QAPRegression")
  expect_length(fit$coefficients, 3)
  # qapspp: lower/larger/abs are matrices with named columns
  expect_true(all(c("x1", "x2") %in% colnames(fit$lower)))
})

# ---- QAPcss GPU qapy ----
test_that("QAPcss GPU works for gaussian qapy", {
  skip_if_not_installed("torch")
  d <- make_gpu_css_data()
  fit <- QAPcss(y ~ x1 + x2, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1, use_gpu = TRUE)
  expect_s3_class(fit, "QAPCSS")
  expect_true(!is.null(fit$lower))
})

# ---- QAPcss GPU qapspp ----
test_that("QAPcss GPU works for gaussian qapspp", {
  skip_if_not_installed("torch")
  d <- make_gpu_css_data()
  fit <- QAPcss(y ~ x1 + x2, data = d,
                family = "gaussian", nullhyp = "qapspp",
                reps = 20, seed = 1, use_gpu = TRUE)
  expect_s3_class(fit, "QAPCSS")
  expect_true(all(c("x1", "x2") %in% colnames(fit$lower)))
})

# ---- GPU single predictor falls back to qapy ----
test_that("QAPglm GPU with single predictor and qapspp falls back to qapy", {
  skip_if_not_installed("torch")
  d <- make_gpu_data()
  fit <- QAPglm(y ~ x1, data = d,
                family = "gaussian", nullhyp = "qapspp",
                reps = 20, seed = 1, use_gpu = TRUE)
  expect_s3_class(fit, "QAPRegression")
  expect_equal(fit$nullhyp, "qapy")
})
