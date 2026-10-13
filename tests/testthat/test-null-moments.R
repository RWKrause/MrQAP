# ============================================================
# Mean and SD of the permutation distribution (fit$null_moments).
#
# The standard: on every path the moments are exactly colMeans() / sd()
# of the draws in fit$null_dist, so they can never describe a different
# distribution from the one the p-values come from.
# ============================================================

nm_data <- function(n = 10, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(y = 0.4 * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
}

nm_binary <- function(n = 12, seed = 2) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  eta <- 0.8 * x1
  list(y = matrix(rbinom(n^2, 1, 1 / (1 + exp(-eta))), n, n),
       x1 = x1, x2 = x2)
}

# Recompute the moments from the draws and compare.
expect_moments_match_draws <- function(mom, nd) {
  expect_false(is.null(mom))
  expect_equal(rownames(mom), c("b_mean", "b_sd", "t_mean", "t_sd"))
  expect_equal(colnames(mom), colnames(nd$b))
  for (st in c("b", "t")) {
    m  <- nd[[st]]
    ok <- colSums(!is.na(m)) > 0
    expect_equal(unname(mom[paste0(st, "_mean"), ok]),
                 unname(colMeans(m[, ok, drop = FALSE], na.rm = TRUE)))
    expect_equal(unname(mom[paste0(st, "_sd"), ok]),
                 unname(apply(m[, ok, drop = FALSE], 2, sd, na.rm = TRUE)))
    expect_true(all(is.na(mom[paste0(st, "_mean"), !ok])))
  }
}


# ---- every permutation path --------------------------------------------

test_that("moments match the draws on the fast OLS path", {
  d <- nm_data()
  for (nh in c("qapy", "qapspp")) {
    f <- QAP(y ~ x1 + x2, data = d, reps = 30, seed = 1, nullhyp = nh)
    expect_moments_match_draws(f$null_moments, f$null_dist)
  }
})

test_that("moments match the draws on the general CPU path", {
  d <- nm_binary()
  for (nh in c("qapy", "qapspp")) {
    f <- QAP(y ~ x1 + x2, data = d, family = "binomial", reps = 20,
             seed = 1, nullhyp = nh)
    expect_moments_match_draws(f$null_moments, f$null_dist)
  }
})

test_that("moments match the draws on the GPU path", {
  skip_if_not(torch_ready(), "torch has no usable LibTorch back end")
  d <- nm_data()
  f <- suppressMessages(
    QAP(y ~ x1 + x2, data = d, reps = 30, seed = 1, use_gpu = TRUE))
  expect_moments_match_draws(f$null_moments, f$null_dist)
})

test_that("moments match the draws for CSS data", {
  set.seed(3)
  a1 <- array(rnorm(6^3), dim = c(6, 6, 6))
  a2 <- array(rnorm(6^3), dim = c(6, 6, 6))
  ya <- 0.5 * a1 + array(rnorm(6^3), dim = c(6, 6, 6))
  f <- QAP(y ~ x1 + x2, data = list(y = ya, x1 = a1, x2 = a2),
           reps = 15, seed = 1)
  expect_moments_match_draws(f$null_moments, f$null_dist)
})

test_that("a comparison model keeps one set of moments per comparison", {
  set.seed(11)
  n <- 10
  y <- matrix(sample(c("FP", "TN", "TP"), n^2, replace = TRUE), n, n)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  comp <- list(commission = c("FP", "TN"), omission = c("TP", "TN"))
  f <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
           family = "binomial", reps = 10, seed = 1, comparison = comp)
  expect_named(f$null_moments, names(comp))
  for (cn in names(comp))
    expect_moments_match_draws(f$null_moments[[cn]], f$null_dist[[cn]])

  tab <- summary(f)$coefficients
  expect_true(all(c("perm_t_mean", "perm_t_sd") %in% names(tab$commission)))
})

test_that("multinomial fits keep their draws and moments", {
  skip_if_not_installed("nnet")
  set.seed(4)
  n  <- 12
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- matrix(sample(c("a", "b", "c"), n^2, replace = TRUE), n, n)
  d  <- list(y = y, x1 = x1, x2 = x2)
  for (nh in c("qapy", "qapspp")) {
    f <- QAP(y ~ x1 + x2, data = d, family = "multinom", reps = 10,
             seed = 1, nullhyp = nh)
    tab <- summary(f)$coefficients
    # Flattened exactly as the coefficient table is.
    expect_equal(colnames(f$null_dist$t), tab$term)
    expect_moments_match_draws(f$null_moments, f$null_dist)
    slopes <- !grepl("(Intercept)", tab$term, fixed = TRUE)
    expect_true(all(is.finite(tab$perm_t_mean[slopes])))
    expect_true(all(is.finite(tab$perm_t_sd[slopes])))
  }
})

test_that("multinomial qapspp draws land in the tested predictor's columns", {
  skip_if_not_installed("nnet")
  set.seed(4)
  n  <- 12
  d  <- list(y = matrix(sample(c("a", "b", "c"), n^2, TRUE), n, n),
             x1 = matrix(rnorm(n^2), n, n), x2 = matrix(rnorm(n^2), n, n))
  f <- QAP(y ~ x1 + x2, data = d, family = "multinom", reps = 10, seed = 1)
  # The observed t of b:x1 under its own qapspp draws must reproduce the
  # stored two-sided p-value; a column mix-up would not.
  tab <- summary(f)$coefficients
  for (nm in c("b:x1", "c:x1", "b:x2", "c:x2")) {
    obs <- tab$statistic[tab$term == nm]
    p   <- mean(abs(f$null_dist$t[, nm]) >= abs(obs))
    expect_equal(p, tab$p_value[tab$term == nm], info = nm)
  }
})


# ---- reporting -----------------------------------------------------------

test_that("summary() reports the moments, NA for the qapspp intercept", {
  f <- QAP(y ~ x1 + x2, data = nm_data(), reps = 30, seed = 1)
  tab <- summary(f)$coefficients
  cols <- c("perm_b_mean", "perm_b_sd", "perm_t_mean", "perm_t_sd")
  expect_true(all(cols %in% names(tab)))
  expect_true(all(is.na(tab[tab$term == "(Intercept)", cols])))
  expect_equal(tab$perm_t_sd[tab$term == "x1"],
               unname(f$null_moments["t_sd", "x1"]))
  expect_equal(as.data.frame(f), tab)

  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("Mean(t*)", out, fixed = TRUE)))
  expect_true(any(grepl("SD(t*)", out, fixed = TRUE)))
})

test_that("non-gaussian tables keep the moment columns", {
  f <- QAP(y ~ x1 + x2, data = nm_binary(), family = "binomial",
           reps = 10, seed = 1)
  tab <- summary(f)$coefficients
  expect_equal(tail(names(tab), 4),
               c("perm_b_mean", "perm_b_sd", "perm_t_mean", "perm_t_sd"))
})

test_that("less_mem drops the draws but keeps the moments", {
  d <- nm_data()
  full <- QAP(y ~ x1 + x2, data = d, reps = 30, seed = 1)
  lean <- QAP(y ~ x1 + x2, data = d, reps = 30, seed = 1, less_mem = TRUE)
  expect_null(lean$null_dist)
  expect_equal(lean$null_moments, full$null_moments)
  expect_equal(summary(lean)$coefficients$perm_t_mean,
               summary(full)$coefficients$perm_t_mean)
})

test_that("a fit without recorded moments still summarises", {
  f <- QAP(y ~ x1 + x2, data = nm_data(), reps = 10, seed = 1)
  f$null_moments <- NULL
  tab <- summary(f)$coefficients
  expect_true(all(is.na(tab$perm_t_mean)))
})


# ---- pooling ---------------------------------------------------------------

test_that("pooled moments equal the moments of the stacked draws", {
  d <- nm_data()
  for (nh in c("qapy", "qapspp")) {
    a <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1, nullhyp = nh)
    b <- QAP(y ~ x1 + x2, data = d, reps = 35, seed = 2, nullhyp = nh)
    ab <- combine_qap_estimates(a, b)
    expect_moments_match_draws(ab$null_moments, ab$null_dist)
  }
})

test_that("moments pool without the draws (less_mem)", {
  d <- nm_data()
  a <- QAP(y ~ x1 + x2, data = d, reps = 20, seed = 1)
  b <- QAP(y ~ x1 + x2, data = d, reps = 35, seed = 2)
  full <- combine_qap_estimates(a, b)

  a$null_dist <- b$null_dist <- NULL
  lean <- combine_qap_estimates(a, b)
  expect_equal(lean$null_moments, full$null_moments)
})
