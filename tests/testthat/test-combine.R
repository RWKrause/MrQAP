# ============================================================
# Pooling repeated runs, and the checks that stop unrelated runs being
# pooled in the first place.
# ============================================================

cmb_data <- function(n = 10, seed = 1, signal = 0.5) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(y = signal * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
}


# ---- what pooling produces ---------------------------------------------

test_that("pooling sums reps and averages the proportions", {
  d <- cmb_data()
  a <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1, nullhyp = "qapy")
  b <- QAP(y ~ x1 + x2, data = d, reps = 60, seed = 2, nullhyp = "qapy")
  ab <- combine_qap_estimates(a, b)

  expect_equal(ab$reps, 80)
  expect_equal(ab$n_valid, 80)
  # Weighted by the permutations each run contributed, not by a flat mean.
  expect_equal(ab$abs, (a$abs * 20 + b$abs * 60) / 80)
  # The baseline fit is untouched by pooling.
  expect_identical(ab$coefficients, a$coefficients)
})

test_that("pooling stacks the permutation draws for confint()", {
  d <- cmb_data()
  a <- QAP(y ~ x1 + x2, data = d, reps = 40, seed = 1, nullhyp = "qapy")
  b <- QAP(y ~ x1 + x2, data = d, reps = 40, seed = 2, nullhyp = "qapy")
  ab <- combine_qap_estimates(a, b)

  # Previously only run 1's draws survived, while reps claimed the total.
  expect_equal(nrow(ab$null_dist$b), 80L)
  expect_equal(nrow(ab$null_dist$t), 80L)
  expect_equal(ab$null_dist$b, rbind(a$null_dist$b, b$null_dist$b))
  expect_true(all(is.finite(confint(ab)[c("x1", "x2"), ])))
})

test_that("pooling weights by permutations that contributed, not requested", {
  # A run whose permutations partly failed has a smaller denominator; the
  # pooled value must use it. Build that case directly.
  d <- cmb_data()
  a <- QAP(y ~ x1 + x2, data = d, reps = 40, seed = 1, nullhyp = "qapy")
  b <- QAP(y ~ x1 + x2, data = d, reps = 40, seed = 2, nullhyp = "qapy")
  b$n_valid <- 10                      # as if 30 of 40 had failed

  ab <- combine_qap_estimates(a, b)
  expect_equal(ab$abs, (a$abs * 40 + b$abs * 10) / 50)
  expect_equal(ab$n_valid, 50)
})

test_that("qapspp pools each predictor with its own weight", {
  d <- cmb_data()
  a <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1, nullhyp = "qapspp")
  b <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 2, nullhyp = "qapspp")

  expect_named(a$n_valid, c("x1", "x2"))
  ab <- combine_qap_estimates(a, b)

  for (v in c("x1", "x2"))
    expect_equal(unname(ab$abs["t", v]),
                 unname((a$abs["t", v] * a$n_valid[v] +
                         b$abs["t", v] * b$n_valid[v]) /
                        (a$n_valid[v] + b$n_valid[v])))
  # The intercept is never tested under qapspp and stays NA.
  expect_true(is.na(ab$abs["t", "(Intercept)"]))
})

test_that("pooling many runs approximates one long run", {
  d <- cmb_data(signal = 0.35)
  runs <- lapply(1:5, function(s)
    QAP(y ~ x1 + x2, data = d, reps = 100, seed = s, nullhyp = "qapy"))
  pooled <- combine_qap_estimates(runs)
  long   <- QAP(y ~ x1 + x2, data = d, reps = 500, seed = 99,
                nullhyp = "qapy")

  expect_equal(pooled$reps, 500)
  # Monte Carlo slack: both estimate the same p-value from 500 draws.
  expect_equal(unname(pooled$abs["t", "x1"]), unname(long$abs["t", "x1"]),
               tolerance = 0.1)
})


# ---- what pooling refuses ----------------------------------------------

test_that("pooling refuses runs fitted to different data", {
  a <- QAP(y ~ x1 + x2, data = cmb_data(seed = 1), reps = 20, seed = 1)
  b <- QAP(y ~ x1 + x2, data = cmb_data(seed = 2), reps = 20, seed = 1)
  # Same shape, same names, same settings -- different networks. Only the
  # coefficients reveal it.
  expect_error(combine_qap_estimates(a, b), "not fitted to the same data")
})

test_that("pooling refuses mismatched settings", {
  d <- cmb_data()
  base <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1, nullhyp = "qapy")

  other_null <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1,
                    nullhyp = "qapspp")
  expect_error(combine_qap_estimates(base, other_null), "same nullhyp")

  robust <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1,
                nullhyp = "qapy", use_robust_errors = TRUE)
  expect_error(combine_qap_estimates(base, robust), "same robust_se")

  undirected <- suppressWarnings(
    QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1, nullhyp = "qapy",
        mode = "undirected"))
  expect_error(combine_qap_estimates(base, undirected), "same mode")
})

test_that("pooling refuses different families and classes", {
  set.seed(3); n <- 10
  x1 <- matrix(rnorm(n^2), n, n)
  db <- list(y = matrix(rbinom(n^2, 1, 0.4), n, n), x1 = x1)
  dg <- list(y = 0.5 * x1 + matrix(rnorm(n^2), n, n), x1 = x1)

  g <- QAP(y ~ x1, data = dg, reps = 15, seed = 1)
  b <- QAP(y ~ x1, data = db, family = "binomial", reps = 15, seed = 1)
  expect_error(combine_qap_estimates(g, b), "same class")
})

test_that("pooling still requires at least two models", {
  a <- QAP(y ~ x1, data = cmb_data(), reps = 10, seed = 1)
  expect_error(combine_qap_estimates(a), "at least two")
})

test_that("a model without n_valid pools with a warning", {
  d <- cmb_data()
  a <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1, nullhyp = "qapy")
  b <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 2, nullhyp = "qapy")
  b$n_valid <- NULL                    # as if fitted by an older version
  expect_warning(combine_qap_estimates(a, b), "before n_valid was recorded")
})


# ---- dimnames validation -----------------------------------------------

test_that("mismatched node orderings are rejected", {
  set.seed(4); n <- 8
  nm <- paste0("a", seq_len(n))
  y  <- matrix(rnorm(n^2), n, n, dimnames = list(nm, nm))
  x1 <- matrix(rnorm(n^2), n, n, dimnames = list(rev(nm), nm))

  # Same actors, different order: the commonest way to get a silently
  # meaningless MRQAP result.
  expect_error(QAP(y ~ x1, data = list(y = y, x1 = x1), reps = 5, seed = 1),
               "different order")
})

test_that("different node sets are rejected", {
  set.seed(4); n <- 8
  y  <- matrix(rnorm(n^2), n, n,
               dimnames = list(paste0("a", 1:n), paste0("a", 1:n)))
  x1 <- matrix(rnorm(n^2), n, n,
               dimnames = list(paste0("z", 1:n), paste0("a", 1:n)))
  expect_error(QAP(y ~ x1, data = list(y = y, x1 = x1), reps = 5, seed = 1),
               "same set of nodes")
})

test_that("matching dimnames pass, and absent ones are not required", {
  set.seed(4); n <- 8
  nm <- paste0("a", seq_len(n))
  y  <- matrix(rnorm(n^2), n, n, dimnames = list(nm, nm))
  x1 <- matrix(rnorm(n^2), n, n, dimnames = list(nm, nm))
  x2 <- matrix(rnorm(n^2), n, n)          # no dimnames at all

  expect_no_error(QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
                      reps = 5, seed = 1))
})

test_that("dimnames are checked per network and per dimension", {
  set.seed(5); n <- 6
  nm <- paste0("n", seq_len(n))
  mk <- function(dn) matrix(rnorm(n^2), n, n, dimnames = dn)

  y  <- list(mk(list(nm, nm)), mk(list(nm, nm)))
  x1 <- list(mk(list(nm, nm)), mk(list(nm, rev(nm))))  # 2nd net, columns

  expect_error(QAP(y ~ x1, data = list(y = y, x1 = x1), reps = 5, seed = 1),
               "column names")
})

test_that("CSS arrays are checked on all three dimensions", {
  set.seed(6); n <- 5
  nm <- paste0("p", seq_len(n))
  y  <- array(rnorm(n^3), c(n, n, n), dimnames = list(nm, nm, nm))
  x1 <- array(rnorm(n^3), c(n, n, n), dimnames = list(nm, nm, rev(nm)))
  expect_error(QAP(y ~ x1, data = list(y = y, x1 = x1), reps = 3, seed = 1),
               "perceiver names")
})
