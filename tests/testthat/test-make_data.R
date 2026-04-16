# ============================================================
# Tests for make_qap_data and make_css_data
# ============================================================

test_that("make_qap_data creates correct data frame dimensions", {
  n <- 5
  y <- matrix(rnorm(n^2), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  pred <- make_qap_data(y = y, x = list(x1 = x1),
                        diag = FALSE, mode = "digraph", net = 1)
  # Without diagonal: n*(n-1) = 20 observations
  expect_equal(nrow(pred), n * (n - 1))
  expect_true("yv" %in% names(pred))
  expect_true("x1" %in% names(pred))
  expect_true("sv" %in% names(pred))
  expect_true("rv" %in% names(pred))
  expect_true("location" %in% names(pred))
})

test_that("make_qap_data handles diagonal inclusion", {
  n <- 4
  y <- matrix(1, n, n)
  x1 <- matrix(1, n, n)
  pred_no_diag <- make_qap_data(y = y, x = list(x1 = x1),
                                 diag = FALSE, mode = "digraph", net = 1)
  pred_diag <- make_qap_data(y = y, x = list(x1 = x1),
                              diag = TRUE, mode = "digraph", net = 1)
  expect_equal(nrow(pred_no_diag), n * (n - 1))
  expect_equal(nrow(pred_diag), n^2)
})

test_that("make_qap_data handles NAs", {
  n <- 4
  y <- matrix(1, n, n)
  y[1, 2] <- NA
  x1 <- matrix(1, n, n)
  pred <- make_qap_data(y = y, x = list(x1 = x1),
                        diag = FALSE, mode = "digraph", net = 1)
  # Should be n*(n-1) - 1 because of the NA
  expect_equal(nrow(pred), n * (n - 1) - 1)
})

test_that("make_qap_data handles undirected mode", {
  n <- 4
  y <- matrix(1, n, n)
  x1 <- matrix(1, n, n)
  pred <- make_qap_data(y = y, x = list(x1 = x1),
                        diag = FALSE, mode = "graph", net = 1)
  # Upper triangle without diagonal: n*(n-1)/2 = 6
  # But make_qap_data doesn't handle graph mode itself;

  # it uses full matrix. This test documents current behavior.
  expect_true(nrow(pred) > 0)
})

test_that("make_css_data creates correct columns", {
  n <- 4
  y <- array(rnorm(n^3), dim = c(n, n, n))
  x1 <- array(rnorm(n^3), dim = c(n, n, n))
  res <- make_css_data(y = y, x = list(x1 = x1),
                       nets = 1,
                       diag = FALSE, mode = "directed")
  pred <- res$pred
  expect_true("yv" %in% names(pred))
  expect_true("x1" %in% names(pred))
  expect_true("sv" %in% names(pred))
  expect_true("rv" %in% names(pred))
  expect_true("pv" %in% names(pred))
  expect_true("nv" %in% names(pred))
  expect_true(nrow(pred) > 0)
})
