# ============================================================
# The Frisch-Waugh-Lovell path used by qapspp on the fast OLS route.
#
# It computes one coefficient and its t-value from a decomposition of the
# untested columns, instead of re-solving the whole design each time. That
# is exact algebra, so these tests hold it to lm() and to the full solver
# at machine precision -- not to a tolerance that would hide a real error.
# ============================================================

fwl_case <- function(n = 120, p = 4, seed = 1) {
  set.seed(seed)
  X <- cbind("(Intercept)" = 1,
             matrix(rnorm(n * p), n, p,
                    dimnames = list(NULL, paste0("x", seq_len(p)))))
  y <- as.vector(X %*% c(0.5, 0.8, -0.4, 0.2, 0.1)[seq_len(p + 1)] +
                   rnorm(n, sd = 0.7))
  list(X = X, y = y)
}


test_that("FWL matches lm() for every column", {
  d <- fwl_case()
  ref <- summary(lm(d$y ~ d$X[, -1]))$coefficients

  for (col in 2:ncol(d$X)) {
    fw  <- qap_fwl_template(d$X, d$y, col, robust = FALSE)
    got <- qap_fwl_solve(fw, d$X[, col])
    expect_equal(got$b, unname(ref[col, 1]), tolerance = 1e-10,
                 info = colnames(d$X)[col])
    expect_equal(got$t, unname(ref[col, 3]), tolerance = 1e-10,
                 info = colnames(d$X)[col])
  }
})

test_that("FWL matches the full solver on a permuted column", {
  d <- fwl_case()
  set.seed(7)

  for (robust in c(FALSE, TRUE)) {
    for (col in 2:ncol(d$X)) {
      fw <- qap_fwl_template(d$X, d$y, col, robust = robust)
      for (i in 1:5) {
        xp <- d$X[sample.int(nrow(d$X)), col]

        Xp <- d$X; Xp[, col] <- xp
        ref <- qap_ols_solve(Xp, d$y, robust = robust)
        got <- qap_fwl_solve(fw, xp)

        expect_equal(got$b, unname(ref$coefficients[col]), tolerance = 1e-10)
        expect_equal(got$t, unname(ref$t[col]), tolerance = 1e-10,
                     info = paste("robust =", robust))
      }
    }
  }
})

test_that("the FWL robust t-value is the HC3 sandwich, not an approximation", {
  d <- fwl_case()
  col <- 3L
  fw  <- qap_fwl_template(d$X, d$y, col, robust = TRUE)
  got <- qap_fwl_solve(fw, d$X[, col])

  # HC3() builds the whole sandwich; FWL specialises it to one coefficient.
  b  <- qr.coef(qr(d$X), d$y)
  se <- HC3(d$X[, -1, drop = FALSE], d$y - as.vector(d$X %*% b))
  expect_equal(got$t, unname(b[col] / se[col]), tolerance = 1e-10)
})

test_that("the leverage decomposition is exact", {
  # h_full = h_Z + x_tilde^2 / sxx is what lets the robust path avoid
  # rebuilding the hat matrix each replication.
  d  <- fwl_case()
  col <- 2L
  fw <- qap_fwl_template(d$X, d$y, col, robust = TRUE)

  x_t <- qr.resid(fw$qrZ, d$X[, col])
  h   <- fw$h_Z + x_t^2 / sum(x_t^2)
  expect_equal(h, rowSums(qr.Q(qr(d$X))^2), tolerance = 1e-10)
})

test_that("qap_compare_one matches compare_perm_to_baseline", {
  # The FWL path never forms the other coefficients, so it uses a
  # single-column comparison. It must produce exactly what the general
  # comparison produces for that column.
  base <- list(coefficients = c("(Intercept)" = 1, x1 = 0.8, x2 = -0.3),
               t            = c("(Intercept)" = 4, x1 = 2.5, x2 = -1.1))
  perm_b <- c("(Intercept)" = 0.9, x1 = 0.4, x2 = -0.9)
  perm_t <- c("(Intercept)" = 3.1, x1 = 1.2, x2 = -3.0)

  for (xi in c("x1", "x2")) {
    ref <- compare_perm_to_baseline(perm_b, perm_t, base, xi = xi)
    got <- qap_compare_one(perm_b[[xi]], perm_t[[xi]],
                           base$coefficients[[xi]], base$t[[xi]])
    expect_equal(got$lower,  ref$lower,  info = xi)
    expect_equal(got$larger, ref$larger, info = xi)
    expect_equal(got$abs,    ref$abs,    info = xi)
    expect_equal(got$draw,   ref$draw,   info = xi)
  }
})

test_that("the fast qapspp path agrees with the general path end to end", {
  set.seed(3); n <- 14
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  d  <- list(y = 0.3 * x1 - 0.2 * x2 + matrix(rnorm(n^2), n, n),
             x1 = x1, x2 = x2)

  run <- function(robust, disable) {
    old <- options(MrQAP.disable_fast_ols = disable)
    on.exit(options(old), add = TRUE)
    QAP(y ~ x1 + x2, data = d, nullhyp = "qapspp", reps = 40,
        seed = 5, use_robust_errors = robust)
  }

  for (robust in c(FALSE, TRUE)) {
    fast <- run(robust, FALSE)
    slow <- run(robust, TRUE)

    # Same permutations, same arithmetic: the p-values must be identical,
    # not merely close.
    expect_equal(fast$abs,    slow$abs,    tolerance = 1e-12)
    expect_equal(fast$lower,  slow$lower,  tolerance = 1e-12)
    expect_equal(fast$larger, slow$larger, tolerance = 1e-12)
    expect_equal(fast$null_dist$b, slow$null_dist$b, tolerance = 1e-10)
    expect_equal(fast$null_dist$t, slow$null_dist$t, tolerance = 1e-10)
  }
})
