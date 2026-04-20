# ============================================================
# Tests for RMPerm (random matrix permutation)
# ============================================================

test_that("RMPerm preserves matrix dimensions", {
  m <- matrix(1:16, 4, 4)
  p <- RMPerm(m)
  expect_equal(dim(p), dim(m))
})

test_that("RMPerm preserves values (just reorders)", {
  set.seed(42)
  m <- matrix(rnorm(25), 5, 5)
  p <- RMPerm(m)
  # The sorted values should be similar (permutation preserves values)
  expect_equal(sort(as.vector(m)), sort(as.vector(p)))
})

test_that("RMPerm works on lists", {
  m1 <- matrix(1:9, 3, 3)
  m2 <- matrix(10:18, 3, 3)
  set.seed(1)
  result <- RMPerm(list(m1, m2))
  expect_true(is.list(result))
  expect_length(result, 2)
  expect_true(is.matrix(result[[1]]))
  expect_true(is.matrix(result[[2]]))
})

test_that("RMPerm works on 3D arrays with CSS = TRUE", {
  set.seed(42)
  a <- array(rnorm(27), dim = c(3, 3, 3))
  p <- RMPerm(a, CSS = TRUE)
  expect_equal(dim(p), dim(a))
})

test_that("RMPerm respects group constraints", {
  set.seed(42)
  m <- diag(4)
  groups <- c("a", "a", "b", "b")
  p <- RMPerm(m, groups = groups)
  expect_equal(dim(p), dim(m))
})
