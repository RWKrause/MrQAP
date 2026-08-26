# ============================================================
# Tests for helper functions: HC3, df_to_mat, combine_qap_estimates,
# array_to_vector, gpu_available
# ============================================================

# ---- HC3 ----
test_that("HC3 returns correct length", {
  X <- matrix(rnorm(30), 10, 3)
  e <- rnorm(10)
  result <- HC3(X, e)
  expect_length(result, 4)  # intercept + 3 columns
})

test_that("HC3 produces positive SEs", {
  X <- matrix(rnorm(50), 10, 5)
  e <- rnorm(10)
  result <- HC3(X, e)
  expect_true(all(result > 0))
})

# ---- array_to_vector ----
test_that("array_to_vector works for directed", {
  a <- array(1:27, dim = c(3, 3, 3))
  v <- array_to_vector(a, mode. = "directed", diag. = TRUE)
  expect_length(v, 27)
})

test_that("array_to_vector works for undirected", {
  a <- array(1:27, dim = c(3, 3, 3))
  v <- array_to_vector(a, mode. = "undirected", diag. = FALSE)
  # 3 perceivers, each upper tri without diag = 3 elements
  expect_length(v, 9)
})

# ---- df_to_mat ----
test_that("df_to_mat creates correct matrix from dataframe", {
  df <- data.frame(
    sender   = c("A", "A", "B"),
    receiver = c("B", "C", "C"),
    weight   = c(1, 2, 3)
  )
  result <- df_to_mat(df, sender = "sender", receiver = "receiver")
  expect_true(is.list(result))
  expect_true("weight" %in% names(result))
  expect_true(is.matrix(result$weight))
  expect_equal(result$weight["A", "B"], 1)
  expect_equal(result$weight["A", "C"], 2)
})

test_that("df_to_mat handles undirected mode", {
  df <- data.frame(
    sender   = c("A", "A"),
    receiver = c("B", "C"),
    weight   = c(1, 2)
  )
  result <- df_to_mat(df, sender = "sender", receiver = "receiver",
                      mode = "undirected")
  expect_equal(result$weight["A", "B"], result$weight["B", "A"])
})

# ---- combine_qap_estimates ----
test_that("combine_qap_estimates works with two runs", {
  # The class must match what QAP() sets -- "QAPRegression" plus the shared
  # "QAP" parent -- not the lowercase spelling no real fit ever carries.
  mk <- function(v1, v2, v3) structure(
    list(lower = matrix(v1, 2, 3), larger = matrix(v2, 2, 3),
         abs = matrix(v3, 2, 3), reps = 100, n_valid = 100,
         nullhyp = "qapy", family = "gaussian", estimator = "standard",
         mode = "directed", diag = FALSE, css = FALSE, multi_mode = FALSE,
         robust_se = FALSE, reference = NULL,
         coefficients = c("(Intercept)" = 1, x1 = 2, x2 = 3), comp = NULL),
    class = c("QAPRegression", "QAP"))
  r1 <- mk(0.1, 0.9, 0.2); r2 <- mk(0.3, 0.7, 0.4)
  combined <- combine_qap_estimates(list(r1, r2))
  expect_equal(combined$reps, 200)
  expect_equal(combined$lower[1, 1], 0.2, tolerance = 0.001)
})

# ---- gpu_available ----
test_that("gpu_available returns logical", {
  result <- gpu_available()
  expect_true(is.logical(result))
})

# ---- GMM moment functions ----
test_that("negbin_moments returns correct dimensions", {
  n <- 50; p <- 3
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Y <- rpois(n, 2)
  theta <- rnorm(p + 1)  # p regression + 1 log(alpha)
  g <- negbin_moments(theta, list(y = Y, x = X))
  expect_equal(nrow(g), n)
  expect_equal(ncol(g), p + 1)  # p score columns + 1 variance moment
})

test_that("zip_moments returns correct dimensions", {
  n <- 50; p <- 3
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Y <- rpois(n, 2)
  Y[sample(n, 15)] <- 0  # add zeros
  theta <- rnorm(p + 1)  # p regression + 1 logit(pi)
  g <- zip_moments(theta, list(y = Y, x = X))
  expect_equal(nrow(g), n)
  expect_equal(ncol(g), p + 1)  # p score columns + 1 zero-indicator moment
})

test_that("negbin_resid returns correct length", {
  n <- 50; p <- 3
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Y <- rpois(n, 2)
  # Fake gmm output
  gmmo <- list(dat = list(y = Y, x = X),
               coefficients = rnorm(p + 1))
  r <- negbin_resid(gmmo)
  expect_length(r, n)
})

test_that("zip_resid returns correct length", {
  n <- 50; p <- 3
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  Y <- rpois(n, 2)
  gmmo <- list(dat = list(y = Y, x = X),
               coefficients = rnorm(p + 1))
  r <- zip_resid(gmmo)
  expect_length(r, n)
})
