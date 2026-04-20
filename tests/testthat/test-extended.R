# ============================================================
# Extended tests: coverage for fixest, multinomial, comparisons,
# combine with comparisons, robust errors
# ============================================================

# --- Helpers ---
make_ext_data <- function(n = 10, seed = 42) {
  set.seed(seed)
  y  <- matrix(rnorm(n^2), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_ext_binary <- function(n = 10, seed = 42) {
  set.seed(seed)
  y  <- matrix(rbinom(n^2, 1, 0.3), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_ext_count <- function(n = 10, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  lambda <- exp(0.3 * x1 + 0.2 * x2)
  y  <- matrix(rpois(n^2, as.vector(lambda)), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}

make_ext_multinom <- function(n = 10, seed = 42) {
  set.seed(seed)
  y  <- matrix(sample(c("A", "B", "C"), n^2, replace = TRUE), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  diag(y) <- diag(x1) <- NA
  list(y = y, x1 = x1)
}

make_ext_comparison <- function(n = 10, seed = 42) {
  set.seed(seed)
  y <- matrix(sample(c("TP", "FP", "TN", "FN"), n^2, replace = TRUE), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  diag(y) <- diag(x1) <- diag(x2) <- NA
  list(y = y, x1 = x1, x2 = x2)
}


# ============================================================
# Fixest tests
# ============================================================

test_that("QAPglm works with fixest FE (gaussian)", {
  skip_if_not_installed("fixest")
  d <- make_ext_data()
  fit <- QAPglm(y ~ x1 | sv, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPRegression")
  expect_true(!is.null(fit$coefficients))
})

test_that("QAPglm works with fixest FE (binomial)", {
  skip_if_not_installed("fixest")
  d <- make_ext_binary()
  fit <- QAPglm(y ~ x1 | sv, data = d,
                family = "binomial", nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
})


# ============================================================
# Multinomial tests
# ============================================================

test_that("QAPglm works with multinomial family", {
  skip_if_not_installed("nnet")
  d <- make_ext_multinom()
  fit <- QAPglm(y ~ x1, data = d,
                family = "multinom", nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$coefficients))
})

test_that("QAPglm multinomial with reference category", {
  skip_if_not_installed("nnet")
  d <- make_ext_multinom()
  fit <- QAPglm(y ~ x1, data = d,
                family = "multinom", nullhyp = "qapy",
                reference = "A",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
})


# ============================================================
# Comparison + qapspp tests
# ============================================================

test_that("QAPglm comparison works with qapspp", {
  d <- make_ext_comparison()
  comp <- list(commission = c("FP", "TN"))
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "gaussian", comparison = comp,
                nullhyp = "qapspp",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$comp))
  # p-values should be list with named elements
  expect_true("commission" %in% names(fit$lower))
})

test_that("QAPglm comparison works with qapy", {
  d <- make_ext_comparison()
  comp <- list(commission = c("FP", "TN"))
  fit <- QAPglm(y ~ x1 + x2, data = d,
                family = "gaussian", comparison = comp,
                nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true("commission" %in% names(fit$lower))
})


# ============================================================
# Robust errors with non-gaussian
# ============================================================

test_that("QAPglm works with robust errors for binomial", {
  d <- make_ext_binary()
  fit <- QAPglm(y ~ x1, data = d,
                family = "binomial", nullhyp = "qapy",
                use_robust_errors = TRUE,
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(fit$robust_se)
})

test_that("QAPglm works with robust errors for poisson", {
  d <- make_ext_count()
  fit <- QAPglm(y ~ x1, data = d,
                family = "poisson", nullhyp = "qapy",
                use_robust_errors = TRUE,
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(fit$robust_se)
})


# ============================================================
# combine_qap_estimates with comparisons
# ============================================================

test_that("combine_qap_estimates works with comparisons", {
  comp <- list(commission = c("FP", "TN"))
  r1 <- list(
    commission = list(
      lower  = matrix(0.1, 2, 3),
      larger = matrix(0.9, 2, 3),
      abs    = matrix(0.2, 2, 3)
    ),
    reps = 100,
    comp = comp
  )
  r2 <- list(
    commission = list(
      lower  = matrix(0.3, 2, 3),
      larger = matrix(0.7, 2, 3),
      abs    = matrix(0.4, 2, 3)
    ),
    reps = 100,
    comp = comp
  )
  class(r1) <- class(r2) <- "QAPGLM"
  combined <- combine_qap_estimates(list(r1, r2))
  expect_equal(combined$reps, 200)
  expect_equal(combined$commission$lower[1, 1], 0.2, tolerance = 0.001)
  expect_equal(combined$commission$larger[1, 1], 0.8, tolerance = 0.001)
  expect_equal(combined$commission$abs[1, 1], 0.3, tolerance = 0.001)
})

test_that("combine_qap_estimates works with res + res2 syntax", {
  r1 <- list(lower = matrix(0.1, 2, 2), larger = matrix(0.9, 2, 2),
             abs = matrix(0.2, 2, 2), reps = 50, comp = NULL)
  r2 <- list(lower = matrix(0.3, 2, 2), larger = matrix(0.7, 2, 2),
             abs = matrix(0.4, 2, 2), reps = 50, comp = NULL)
  class(r1) <- class(r2) <- "QAPRegression"
  combined <- combine_qap_estimates(r1, r2)
  expect_equal(combined$reps, 100)
})
