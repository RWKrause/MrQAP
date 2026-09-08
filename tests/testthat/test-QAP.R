# ============================================================
# The unified QAP() entry point: detection, assertions, and the
# guarantee that the returned object has the same shape whatever
# the input shape was.
# ============================================================

mat_data <- function(n = 10, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(y = x1 - 0.5 * x2 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
}

css_data <- function(n = 6, seed = 2) {
  set.seed(seed)
  mk <- function() array(rnorm(n^3), dim = c(n, n, n))
  x1 <- mk(); x2 <- mk()
  list(y = x1 - 0.5 * x2 + array(rnorm(n^3), dim = c(n, n, n)),
       x1 = x1, x2 = x2)
}

# ---- CSS detection -----------------------------------------------------

test_that("QAP detects matrix data and CSS data from the dimensions", {
  m <- QAP(y ~ x1 + x2, data = mat_data(), reps = 5, seed = 1)
  expect_false(m$css)
  expect_s3_class(m, "QAPRegression")

  s <- QAP(y ~ x1 + x2, data = css_data(), reps = 5, seed = 1)
  expect_true(s$css)
  expect_s3_class(s, "QAPCSS")
})

test_that("detection looks inside a list of networks", {
  d <- mat_data(n = 6)
  dl <- list(y = list(d$y, d$y), x1 = list(d$x1, d$x1),
             x2 = list(d$x2, d$x2))
  expect_false(QAP(y ~ x1 + x2, data = dl, reps = 5, seed = 1)$css)

  c1 <- css_data(n = 5)
  cl <- list(y = list(c1$y, c1$y), x1 = list(c1$x1, c1$x1),
             x2 = list(c1$x2, c1$x2))
  expect_true(QAP(y ~ x1 + x2, data = cl, reps = 5, seed = 1)$css)
})

test_that("css = TRUE/FALSE asserts rather than guesses", {
  expect_error(QAP(y ~ x1 + x2, data = mat_data(), css = TRUE, reps = 5),
               "css = TRUE")
  expect_error(QAP(y ~ x1 + x2, data = css_data(), css = FALSE, reps = 5),
               "css = FALSE")
  # and agrees with detection when it is right
  expect_s3_class(QAP(y ~ x1 + x2, data = mat_data(), css = FALSE,
                      reps = 5, seed = 1), "QAPRegression")
  expect_s3_class(QAP(y ~ x1 + x2, data = css_data(), css = TRUE,
                      reps = 5, seed = 1), "QAPCSS")
})

test_that("detect_css rejects data that are neither", {
  expect_error(detect_css(1:10, dep = "y"), "dimension")
  expect_error(detect_css(matrix(1, 2, 2), css = NA), "must be TRUE")
})

# ---- one object shape --------------------------------------------------

test_that("matrix and CSS fits carry the same top-level fields", {
  m <- QAP(y ~ x1 + x2, data = mat_data(), reps = 5, seed = 1)
  s <- QAP(y ~ x1 + x2, data = css_data(), reps = 5, seed = 1)

  # CSS used to expose results only under $base.
  for (f in c("coefficients", "t", "lower", "larger", "abs",
              "nullhyp", "family", "mode", "reps", "css", "random",
              "robust_se", "estimator", "base")) {
    expect_false(is.null(m[[f]]), info = paste("matrix fit lacks", f))
    expect_false(is.null(s[[f]]), info = paste("CSS fit lacks", f))
  }
  expect_equal(s$coefficients, s$base$coefficients)
  expect_equal(s$t, s$base$t)
})

test_that("random-intercept flags are reported for both shapes", {
  m <- QAP(y ~ x1 + x2, data = mat_data(), reps = 5, seed = 1)
  expect_named(m$random, c("sender", "receiver", "perceiver", "nets"))
  expect_false(any(m$random))
})

# ---- argument handling -------------------------------------------------

test_that("nullhyp defaults to qapspp for both shapes", {
  expect_equal(QAP(y ~ x1 + x2, data = mat_data(), reps = 5, seed = 1)$nullhyp,
               "qapspp")
  expect_equal(QAP(y ~ x1 + x2, data = css_data(), reps = 5, seed = 1)$nullhyp,
               "qapspp")
})

test_that("a single predictor falls back to qapy", {
  expect_equal(QAP(y ~ x1, data = mat_data(), reps = 5, seed = 1)$nullhyp,
               "qapy")
})

test_that("bad argument values are rejected by name", {
  d <- mat_data()
  expect_error(QAP(y ~ x1 + x2, data = d, family = "gausian", reps = 5))
  expect_error(QAP(y ~ x1 + x2, data = d, mode = "diriected", reps = 5))
  expect_error(QAP(y ~ x1 + x2, data = d, nullhyp = "qap", reps = 5))
  expect_error(QAP(y ~ x1 + x2, data = d, estimator = "ols", reps = 5))
})

test_that("random_intercept_perceiver is ignored for matrix data", {
  expect_warning(
    QAP(y ~ x1 + x2, data = mat_data(), random_intercept_perceiver = TRUE,
        reps = 5, seed = 1),
    "CSS data only"
  )
})

test_that("a missing dependent variable is reported clearly", {
  expect_error(QAP(nope ~ x1, data = mat_data(), reps = 5), "not found")
})

# ---- groups ------------------------------------------------------------

test_that("group vectors are length-checked for both shapes", {
  expect_error(QAP(y ~ x1 + x2, data = mat_data(n = 10),
                   groups = c(1, 1, 2), reps = 5), "10 nodes")
  expect_error(QAP(y ~ x1 + x2, data = css_data(n = 6),
                   groups = c(1, 1, 2), reps = 5), "6 nodes")
})

test_that("per-network group lists are length-checked", {
  d <- mat_data(n = 6)
  dl <- list(y = list(d$y, d$y), x1 = list(d$x1, d$x1),
             x2 = list(d$x2, d$x2))
  expect_error(QAP(y ~ x1 + x2, data = dl,
                   groups = list(rep(1:2, each = 3)), reps = 5),
               "list of length 1")
  expect_error(QAP(y ~ x1 + x2, data = dl,
                   groups = list(rep(1:2, each = 3), c(1, 2)), reps = 5),
               "network 2")
})

test_that("fit$groups is NULL unless the user supplied groups", {
  expect_null(QAP(y ~ x1 + x2, data = css_data(), reps = 5, seed = 1)$groups)
  expect_null(QAP(y ~ x1 + x2, data = mat_data(), reps = 5, seed = 1)$groups)
  g <- QAP(y ~ x1 + x2, data = mat_data(n = 10),
           groups = rep(1:2, each = 5), reps = 5, seed = 1)
  expect_false(is.null(g$groups))
})

# ---- the RNG state is left alone --------------------------------------

test_that("QAP restores the caller's RNG state", {
  set.seed(4242)
  before <- .Random.seed
  invisible(QAP(y ~ x1 + x2, data = mat_data(), reps = 5, seed = 99))
  expect_identical(.Random.seed, before)
})

test_that("seeded runs are reproducible", {
  a <- QAP(y ~ x1 + x2, data = mat_data(), reps = 10, seed = 7)
  b <- QAP(y ~ x1 + x2, data = mat_data(), reps = 10, seed = 7)
  expect_equal(a$abs, b$abs)
  expect_equal(a$coefficients, b$coefficients)
})

# ---- the future plan is left alone ------------------------------------

test_that("QAP does not change the active future plan when ncores is NULL", {
  future::plan(future::sequential)
  before <- class(future::plan())
  invisible(QAP(y ~ x1 + x2, data = mat_data(), reps = 5, seed = 1))
  expect_identical(class(future::plan()), before)
})
