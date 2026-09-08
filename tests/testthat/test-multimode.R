# ============================================================
# Multi-mode data: each array dimension is its own node set.
#
# The permutation must then draw an independent node order per
# dimension, and everything that assumes a square/cubic shape --
# self-ties, reciprocity, the index vectors -- must adapt or be
# rejected.
# ============================================================

# A matrix whose value encodes its own position, so a permutation can be
# read back off the result: m[i, j] == i * 100 + j.
tag_matrix <- function(nr, nc) outer(seq_len(nr), seq_len(nc),
                                     function(i, j) i * 100 + j)
row_order <- function(p) p[, 1] %/% 100
col_order <- function(p) p[1, ] %% 100

# Same idea in 3D: a[i, j, k] == i * 10000 + j * 100 + k.
tag_array <- function(d) {
  array(slice.index(array(0, d), 1) * 10000 +
        slice.index(array(0, d), 2) * 100 +
        slice.index(array(0, d), 3), dim = d)
}

# ---- RMPerm ------------------------------------------------------------

test_that("one-mode permutation applies the same order to rows and columns", {
  set.seed(1)
  m <- tag_matrix(6, 6)
  for (i in 1:30) {
    p <- RMPerm(m)
    expect_identical(row_order(p), col_order(p))
  }
})

test_that("multi-mode permutation draws rows and columns independently", {
  set.seed(2)
  m <- tag_matrix(6, 6)   # square, but declared multi-mode
  same <- logical(50)
  for (i in seq_along(same)) {
    p <- RMPerm(m, multi_mode = TRUE)
    same[i] <- identical(row_order(p), col_order(p))
  }
  # Independent orders coincide only by chance (1/720 per draw).
  expect_false(all(same))
})

test_that("multi-mode permutation handles rectangular matrices", {
  set.seed(3)
  m <- tag_matrix(3, 5)
  for (i in 1:30) {
    p <- RMPerm(m, multi_mode = TRUE)
    expect_equal(dim(p), c(3L, 5L))
    expect_setequal(as.vector(p), as.vector(m))
    expect_setequal(row_order(p), 1:3)
    expect_setequal(col_order(p), 1:5)
  }
})

test_that("one-mode CSS shares one order across all three dimensions", {
  set.seed(4)
  a <- tag_array(c(5, 5, 5))
  for (i in 1:20) {
    p  <- RMPerm(a, CSS = TRUE)
    o1 <- p[, 1, 1] %/% 10000
    o2 <- (p[1, , 1] %/% 100) %% 100
    o3 <- p[1, 1, ] %% 100
    expect_identical(o1, o2)
    expect_identical(o2, o3)
  }
})

test_that("multi-mode CSS uses three independent orders", {
  set.seed(5)
  a <- tag_array(c(4, 5, 6))
  for (i in 1:20) {
    p <- RMPerm(a, CSS = TRUE, multi_mode = TRUE)
    expect_equal(dim(p), c(4L, 5L, 6L))
    expect_setequal(as.vector(p), as.vector(a))
    expect_setequal(p[, 1, 1] %/% 10000, 1:4)
    expect_setequal((p[1, , 1] %/% 100) %% 100, 1:5)
    expect_setequal(p[1, 1, ] %% 100, 1:6)
  }
})

test_that("one-mode still requires a square or cubic shape", {
  expect_error(RMPerm(matrix(1, 3, 4)), "square")
  expect_error(RMPerm(array(1, c(2, 3, 4)), CSS = TRUE), "square")
  expect_silent(RMPerm(matrix(1, 3, 4), multi_mode = TRUE))
})

test_that("multi-mode groups are given per dimension", {
  set.seed(6)
  m  <- tag_matrix(4, 6)
  gr <- c(1, 1, 2, 2)
  gc <- c(1, 1, 1, 2, 2, 2)

  for (i in 1:30) {
    p <- RMPerm(m, groups = list(sender = gr, receiver = gc),
                multi_mode = TRUE)
    expect_equal(gr[row_order(p)], gr)
    expect_equal(gc[col_order(p)], gc)
  }
})

test_that("group specifications are validated", {
  m <- tag_matrix(4, 6)
  expect_error(RMPerm(m, groups = list(c(1, 1, 2, 2)), multi_mode = TRUE),
               "must be named")
  expect_error(RMPerm(m, groups = list(sender = c(1, 1, 2, 2),
                                       nonsense = 1:6),
                      multi_mode = TRUE),
               "Unknown grouping name")
  expect_error(RMPerm(m, groups = list(sender = c(1, 2)), multi_mode = TRUE),
               "dimension 1 has 4 nodes")
  expect_error(RMPerm(m, groups = c(1, 1, 2, 2), multi_mode = TRUE),
               "one grouping per dimension")
})

# ---- array_to_vector ---------------------------------------------------

test_that("array_to_vector keeps every cell of a non-cubic array", {
  a <- array(1:24, dim = c(2, 3, 4))
  v <- array_to_vector(a, mode. = "directed", diag. = FALSE)
  expect_length(v, 24)
  expect_equal(v, as.vector(a))
})

test_that("array_to_vector undirected takes the upper triangle per slice", {
  n <- 4
  a <- tag_array(c(n, n, n))
  v <- array_to_vector(a, mode. = "undirected", diag. = FALSE)
  expect_length(v, n * (n - 1) / 2 * n)
  # every retained cell has sender < receiver
  sen <- v %/% 10000
  rec <- (v %/% 100) %% 100
  expect_true(all(sen < rec))
})

# ---- detection ---------------------------------------------------------

test_that("multi_mode is inferred only when the dimensions differ", {
  expect_false(detect_multi_mode(matrix(1, 5, 5)))
  expect_true(detect_multi_mode(matrix(1, 3, 5)))
  expect_false(detect_multi_mode(array(1, c(4, 4, 4))))
  expect_true(detect_multi_mode(array(1, c(4, 5, 4))))
  # a list of networks is probed by its first element
  expect_true(detect_multi_mode(list(matrix(1, 3, 5), matrix(1, 3, 5))))
})

test_that("multi_mode = FALSE on ragged data is an error, not a guess", {
  expect_error(detect_multi_mode(matrix(1, 3, 5), multi_mode = FALSE),
               "must be square")
  expect_true(detect_multi_mode(matrix(1, 5, 5), multi_mode = TRUE))
  expect_error(detect_multi_mode(matrix(1, 5, 5), multi_mode = NA),
               "must be TRUE")
})

# ---- QAP() end to end --------------------------------------------------

two_mode_data <- function(nr = 8, nc = 5, seed = 11) {
  set.seed(seed)
  x1 <- matrix(rnorm(nr * nc), nr, nc)
  x2 <- matrix(rnorm(nr * nc), nr, nc)
  list(y = 1 + 1.5 * x1 - 0.7 * x2 + matrix(rnorm(nr * nc, sd = .4), nr, nc),
       x1 = x1, x2 = x2)
}

test_that("a two-mode fit uses every cell: no diagonal is dropped", {
  d <- two_mode_data(nr = 8, nc = 5)
  pred <- make_qap_data(y = d$y, x = list(x1 = d$x1), diag = FALSE,
                        mode = "digraph", net = 1, multi_mode = TRUE)
  expect_equal(nrow(pred), 8 * 5)
  expect_setequal(pred$sv, 1:8)
  expect_setequal(pred$rv, 1:5)
})

test_that("one-mode data still drops the diagonal", {
  n <- 6
  m <- matrix(rnorm(n^2), n, n)
  pred <- make_qap_data(y = m, x = list(x1 = m), diag = FALSE,
                        mode = "digraph", net = 1, multi_mode = FALSE)
  expect_equal(nrow(pred), n * (n - 1))
})

test_that("QAP fits two-mode data and matches lm()", {
  d <- two_mode_data()
  fit <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1)

  expect_true(fit$multi_mode)
  ref <- lm(as.vector(d$y) ~ as.vector(d$x1) + as.vector(d$x2))
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
  expect_false(anyNA(fit$abs[2, c("x1", "x2")]))
})

test_that("QAP rejects what multi-mode makes meaningless", {
  d <- two_mode_data()
  expect_error(QAP(y ~ x1 + x2, data = d, diag = TRUE, reps = 5),
               "no self-ties")
  expect_error(QAP(y ~ x1 + x2, data = d, mode = "undirected", reps = 5),
               "no reciprocal dyad")
})

test_that("QAP accepts per-dimension groups for two-mode data", {
  d <- two_mode_data(nr = 8, nc = 6)
  fit <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1,
             groups = list(sender = rep(1:2, each = 4),
                           receiver = rep(1:3, each = 2)))
  expect_true(fit$multi_mode)
  expect_false(is.null(fit$groups))
  expect_false(anyNA(fit$abs[2, c("x1", "x2")]))
})

test_that("QAP handles a three-mode CSS array", {
  set.seed(12)
  d <- c(4, 5, 6)
  mk <- function() array(rnorm(prod(d)), dim = d)
  x1 <- mk(); x2 <- mk()
  dat <- list(y = x1 - 0.5 * x2 + array(rnorm(prod(d), sd = .4), dim = d),
              x1 = x1, x2 = x2)

  fit <- QAP(y ~ x1 + x2, data = dat, reps = 10, seed = 1)
  expect_true(fit$css)
  expect_true(fit$multi_mode)

  ref <- lm(as.vector(dat$y) ~ as.vector(x1) + as.vector(x2))
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
})

test_that("square two-mode data must be declared", {
  d <- two_mode_data(nr = 6, nc = 6)
  # Square: detection cannot tell, so it assumes one-mode and drops the
  # diagonal.
  auto <- QAP(y ~ x1 + x2, data = d, reps = 10, seed = 1)
  expect_false(auto$multi_mode)

  told <- QAP(y ~ x1 + x2, data = d, multi_mode = TRUE, reps = 10, seed = 1)
  expect_true(told$multi_mode)
  # Declaring it keeps the 6 diagonal cells that one-mode mode discards.
  expect_false(isTRUE(all.equal(auto$coefficients, told$coefficients)))
})
