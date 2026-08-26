# ============================================================
# The fast OLS path must be an optimisation, never a different model.
#
# Every case below is run twice with the same seed -- once through
# qap_ols_perms(), once through the general formula-based engine -- and
# the two must agree exactly. option(MrQAP.disable_fast_ols) forces the
# general path.
# ============================================================

both_paths <- function(expr) {
  e <- substitute(expr)
  pe <- parent.frame()
  run <- function(disable) {
    old <- options(MrQAP.disable_fast_ols = disable)
    on.exit(options(old))
    suppressWarnings(suppressMessages(eval(e, pe)))
  }
  list(fast = run(FALSE), slow = run(TRUE))
}

strip <- function(fit) {
  fit[c("coefficients", "t", "lower", "larger", "abs",
        "r.squared", "adj.r.squared")]
}

wk <- function(n = 12, seed = 31, na = FALSE) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- 0.12 * x1 + matrix(rnorm(n^2), n, n)   # weak: p-values informative
  if (na) y[3, 4] <- NA
  list(y = y, x1 = x1, x2 = x2)
}

test_that("fast and general paths agree for qapy", {
  r <- both_paths(QAP(y ~ x1 + x2, data = wk(), nullhyp = "qapy",
                      reps = 30, seed = 1))
  expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12)
  # and the case is informative, not saturated at 0
  expect_gt(max(r$fast$abs[2, ]), 0)
  expect_lt(min(r$fast$abs[2, ]), 1)
})

test_that("fast and general paths agree for qapspp", {
  r <- both_paths(QAP(y ~ x1 + x2, data = wk(), nullhyp = "qapspp",
                      reps = 30, seed = 2))
  expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12)
})

test_that("fast and general paths agree with robust standard errors", {
  # Regression test: the fast path once ignored use_robust_errors and
  # compared ordinary permutation t-values against an HC3 baseline.
  for (nh in c("qapy", "qapspp")) {
    r <- both_paths(QAP(y ~ x1 + x2, data = wk(), nullhyp = nh,
                        use_robust_errors = TRUE, reps = 30, seed = 3))
    expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12,
                 info = nh)
  }
})

test_that("fast and general paths agree for undirected data", {
  set.seed(4)
  n <- 12
  sym <- function() {
    m <- matrix(rnorm(n^2), n, n)
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  x1 <- sym(); x2 <- sym()
  d <- list(y = 0.12 * x1 + sym(), x1 = x1, x2 = x2)

  r <- both_paths(QAP(y ~ x1 + x2, data = d, mode = "undirected",
                      nullhyp = "qapspp", reps = 30, seed = 5))
  expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12)
})

test_that("fast and general paths agree with the diagonal included", {
  r <- both_paths(QAP(y ~ x1 + x2, data = wk(), diag = TRUE,
                      nullhyp = "qapy", reps = 30, seed = 6))
  expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12)
})

test_that("fast and general paths agree with grouped permutation", {
  r <- both_paths(QAP(y ~ x1 + x2, data = wk(), groups = rep(1:3, each = 4),
                      nullhyp = "qapspp", reps = 30, seed = 7))
  expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12)
})

test_that("fast and general paths agree for multiple networks", {
  set.seed(8)
  n <- 8
  mk <- function() matrix(rnorm(n^2), n, n)
  x1 <- list(mk(), mk()); x2 <- list(mk(), mk())
  y  <- Map(function(a) 0.12 * a + mk(), x1)
  d  <- list(y = y, x1 = x1, x2 = x2)

  r <- both_paths(QAP(y ~ x1 + x2, data = d, nullhyp = "qapspp",
                      reps = 30, seed = 9))
  expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12)
})

test_that("fast and general paths agree for CSS", {
  set.seed(10)
  n <- 7
  mk <- function() array(rnorm(n^3), dim = c(n, n, n))
  x1 <- mk(); x2 <- mk()
  d <- list(y = 0.12 * x1 + mk(), x1 = x1, x2 = x2)

  for (nh in c("qapy", "qapspp")) {
    r <- both_paths(QAP(y ~ x1 + x2, data = d, nullhyp = nh,
                        reps = 20, seed = 11))
    expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12, info = nh)
  }
})

test_that("fast and general paths agree for two-mode data", {
  set.seed(12)
  nr <- 9; nc <- 6
  x1 <- matrix(rnorm(nr * nc), nr, nc)
  x2 <- matrix(rnorm(nr * nc), nr, nc)
  d <- list(y = 0.15 * x1 + matrix(rnorm(nr * nc), nr, nc), x1 = x1, x2 = x2)

  r <- both_paths(QAP(y ~ x1 + x2, data = d, nullhyp = "qapspp",
                      reps = 30, seed = 13))
  expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12)
})

# ---- eligibility -------------------------------------------------------

test_that("the fast path declines cases it cannot handle", {
  d <- wk()
  ok <- function(...) qap_ols_eligible(..., data = d, vars = c("y", "x1", "x2"))

  expect_true(ok("gaussian", "standard", FALSE, FALSE, NULL))
  expect_false(ok("binomial", "standard", FALSE, FALSE, NULL))
  expect_false(ok("gaussian", "gmm",      FALSE, FALSE, NULL))
  expect_false(ok("gaussian", "standard", TRUE,  FALSE, NULL))   # random
  expect_false(ok("gaussian", "standard", FALSE, TRUE,  NULL))   # fixest
  expect_false(ok("gaussian", "standard", FALSE, FALSE,
                  list(a = c("x", "y"))))                        # comparison
})

test_that("the fast path declines data containing NAs", {
  # With NAs the set of valid cells shifts under permutation, so the
  # precomputed index would be wrong.
  dn <- wk(na = TRUE)
  expect_false(qap_ols_eligible("gaussian", "standard", FALSE, FALSE, NULL,
                                dn, c("y", "x1", "x2")))
  dl <- list(y = list(wk()$y, wk(na = TRUE)$y),
             x1 = list(wk()$x1, wk()$x1))
  expect_false(qap_ols_eligible("gaussian", "standard", FALSE, FALSE, NULL,
                                dl, c("y", "x1")))
})

test_that("NA data still give identical results on the general path", {
  r <- both_paths(QAP(y ~ x1 + x2, data = wk(na = TRUE), nullhyp = "qapspp",
                      reps = 30, seed = 14))
  expect_equal(strip(r$fast), strip(r$slow), tolerance = 1e-12)
})

# ---- the solver itself -------------------------------------------------

test_that("qap_ols_solve matches lm() and its HC3 standard errors", {
  set.seed(15)
  N <- 300
  X <- cbind("(Intercept)" = 1, a = rnorm(N), b = rnorm(N))
  y <- 1 + 0.5 * X[, "a"] - 0.3 * X[, "b"] + rnorm(N, sd = abs(X[, "a"]))

  ref <- lm(y ~ X[, "a"] + X[, "b"])
  f   <- qap_ols_solve(X, y)
  expect_equal(unname(f$coefficients), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(unname(f$t), unname(summary(ref)$coefficients[, 3]),
               tolerance = 1e-10)

  fr <- qap_ols_solve(X, y, robust = TRUE)
  expect_equal(unname(fr$t),
               unname(coef(ref) / HC3(X[, -1, drop = FALSE], residuals(ref))),
               tolerance = 1e-10)
})
