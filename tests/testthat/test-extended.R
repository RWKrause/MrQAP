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
  fit <- QAP(y ~ x1 | sv, data = d,
                family = "gaussian", nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPRegression")
  expect_true(!is.null(fit$coefficients))
})

test_that("QAPglm works with fixest FE (binomial)", {
  skip_if_not_installed("fixest")
  d <- make_ext_binary()
  fit <- QAP(y ~ x1 | sv, data = d,
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
  fit <- QAP(y ~ x1, data = d,
                family = "multinom", nullhyp = "qapy",
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(!is.null(fit$coefficients))
})

test_that("QAPglm multinomial with reference category", {
  skip_if_not_installed("nnet")
  d <- make_ext_multinom()
  fit <- QAP(y ~ x1, data = d,
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
  fit <- QAP(y ~ x1 + x2, data = d,
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
  fit <- QAP(y ~ x1 + x2, data = d,
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
  fit <- QAP(y ~ x1, data = d,
                family = "binomial", nullhyp = "qapy",
                use_robust_errors = TRUE,
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(fit$robust_se)
})

test_that("QAPglm works with robust errors for poisson", {
  d <- make_ext_count()
  fit <- QAP(y ~ x1, data = d,
                family = "poisson", nullhyp = "qapy",
                use_robust_errors = TRUE,
                reps = 20, seed = 1)
  expect_s3_class(fit, "QAPGLM")
  expect_true(fit$robust_se)
})


# ============================================================
# combine_qap_estimates with comparisons
# ============================================================

# QAP()/QAP() store comparison results as fit$lower[[comparison]] --
# a list of K matrices, one per comparison -- NOT as fit[[comparison]]$lower.
# The fixtures below mirror the real layout.
# Stub fits for the pooling arithmetic. These must carry the same fields a
# real fit does -- including the shared "QAP" parent class, the settings
# combine_qap_estimates() checks, and n_valid -- or they stop testing the
# code path QAP() actually produces. See the drift guard below.
make_comp_fit <- function(lo, la, ab, reps, comps = c("commission", "omission")) {
  mk <- function(v) setNames(lapply(comps, function(i) matrix(v, 2, 3)), comps)
  structure(list(lower = mk(lo), larger = mk(la), abs = mk(ab),
                 reps = reps, n_valid = reps, nullhyp = "qapy",
                 family = "binomial", estimator = "standard",
                 mode = "directed", diag = FALSE, css = FALSE,
                 multi_mode = FALSE, robust_se = FALSE, reference = NULL,
                 coefficients = c("(Intercept)" = 0.1, x1 = 0.2, x2 = 0.3),
                 comp = setNames(lapply(comps, function(i) c("FP", "TN")),
                                 comps)),
            class = c("QAPGLM", "QAP"))
}

test_that("combine_qap_estimates pools reps within each comparison", {
  r1 <- make_comp_fit(0.1, 0.9, 0.2, 100)
  r2 <- make_comp_fit(0.3, 0.7, 0.4, 100)

  combined <- combine_qap_estimates(list(r1, r2))

  expect_equal(combined$reps, 200)
  # Each of the K comparisons is pooled separately and all K survive.
  expect_named(combined$lower, c("commission", "omission"))
  for (cm in c("commission", "omission")) {
    expect_equal(combined$lower[[cm]][1, 1],  0.2, tolerance = 0.001)
    expect_equal(combined$larger[[cm]][1, 1], 0.8, tolerance = 0.001)
    expect_equal(combined$abs[[cm]][1, 1],    0.3, tolerance = 0.001)
  }
  # No spurious top-level element named after a comparison.
  expect_null(combined$commission)
})

test_that("combine_qap_estimates weights by reps", {
  r1 <- make_comp_fit(0.0, 1.0, 0.0, 300)
  r2 <- make_comp_fit(1.0, 0.0, 1.0, 100)

  combined <- combine_qap_estimates(list(r1, r2))

  expect_equal(combined$reps, 400)
  expect_equal(combined$lower[["commission"]][1, 1], 0.25, tolerance = 1e-8)
})

test_that("combine_qap_estimates rejects incompatible models", {
  r1 <- make_comp_fit(0.1, 0.9, 0.2, 100)
  r2 <- make_comp_fit(0.3, 0.7, 0.4, 100, comps = c("commission", "other"))
  expect_error(combine_qap_estimates(list(r1, r2)), "same comparisons")

  r3 <- make_comp_fit(0.3, 0.7, 0.4, 100)
  r3$nullhyp <- "qapspp"
  expect_error(combine_qap_estimates(list(r1, r3)), "same nullhyp")

  expect_error(combine_qap_estimates(list(r1)), "at least two")
})

test_that("combine_qap_estimates matches the layout QAPglm produces", {
  # Guards against the fixtures above drifting from the real object.
  set.seed(11)
  n <- 8
  y <- matrix(sample(c("FP", "TN", "TP"), n^2, replace = TRUE), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
                family = "binomial", reps = 5, seed = 1,
                comparison = list(commission = c("FP", "TN")))

  expect_true(is.list(fit$lower))
  expect_named(fit$lower, "commission")
  expect_true(is.matrix(fit$lower[["commission"]]))
  expect_null(fit$commission)
  # The stubs above must carry the same class and the same fields that
  # combine_qap_estimates() reads, or they test nothing real.
  expect_true(inherits(fit, "QAP"))
  for (f in c("reps", "n_valid", "nullhyp", "family", "estimator", "mode",
              "diag", "css", "multi_mode", "robust_se"))
    expect_false(is.null(fit[[f]]), info = f)
})

test_that("combine_qap_estimates works with res + res2 syntax", {
  mk <- function(v1, v2, v3) structure(
    list(lower = matrix(v1, 2, 2), larger = matrix(v2, 2, 2),
         abs = matrix(v3, 2, 2), reps = 50, n_valid = 50,
         nullhyp = "qapy", family = "gaussian", estimator = "standard",
         mode = "directed", diag = FALSE, css = FALSE, multi_mode = FALSE,
         robust_se = FALSE, reference = NULL,
         coefficients = c("(Intercept)" = 1, x1 = 2), comp = NULL),
    class = c("QAPRegression", "QAP"))
  combined <- combine_qap_estimates(mk(0.1, 0.9, 0.2), mk(0.3, 0.7, 0.4))
  expect_equal(combined$reps, 100)
})
