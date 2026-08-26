# ============================================================
# The GPU path, exercised for real.
#
# test-gpu.R covers the pure-R scaffolding -- eligibility, batch sizing,
# permuted-column construction -- which runs without torch. This file runs
# the torch algebra itself, on the CPU device and, where the machine has
# one, on CUDA. Every test skips cleanly when torch is absent.
#
# The standard the GPU path is held to is the CPU solver: same design, same
# response, same numbers. That is a stronger claim than "it runs", and it
# is the one that matters, because a silent disagreement here would move
# p-values rather than crash.
# ============================================================

skip_no_torch <- function() skip_if_not_installed("torch")

# Devices worth testing on this machine: always the CPU torch backend,
# plus CUDA when it is actually there.
torch_devices <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) return(character(0))
  c("cpu", if (isTRUE(try(torch::cuda_is_available(), silent = TRUE))) "cuda")
}

gt_data <- function(n = 12, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  x3 <- matrix(rnorm(n^2), n, n)
  list(y = 0.5 * x1 - 0.3 * x2 + matrix(rnorm(n^2), n, n),
       x1 = x1, x2 = x2, x3 = x3)
}

gt_tmpl <- function(d, main = c("x1", "x2", "x3")) {
  qap_ols_template(d, "y", main, FALSE, "directed", FALSE, FALSE, FALSE)
}


# ---- the torch algebra against the CPU solver ---------------------------

test_that("gpu_solve_fixed_x reproduces the CPU solver", {
  skip_no_torch()
  d <- gt_data(); tmpl <- gt_tmpl(d)

  set.seed(4)
  B <- 7
  Y <- vapply(seq_len(B),
              function(j) tmpl$y[sample.int(tmpl$n_obs)], numeric(tmpl$n_obs))

  for (dev in torch_devices()) {
    sol <- gpu_solve_fixed_x(tmpl$X, Y, device = dev)
    expect_equal(dim(sol$coefficients), c(tmpl$p, B), info = dev)

    for (j in seq_len(B)) {
      ref <- qap_ols_solve(tmpl$X, Y[, j])
      expect_equal(sol$coefficients[, j], unname(ref$coefficients),
                   tolerance = 1e-9, info = paste(dev, j))
      expect_equal(sol$t[, j], unname(ref$t),
                   tolerance = 1e-9, info = paste(dev, j))
    }
  }
})

test_that("gpu_solve_varying_x reproduces the CPU solver", {
  skip_no_torch()
  d <- gt_data(); tmpl <- gt_tmpl(d)
  col <- match("x1", colnames(tmpl$X))

  set.seed(5)
  B  <- 7
  Pc <- vapply(seq_len(B),
               function(j) tmpl$X[sample.int(tmpl$n_obs), col],
               numeric(tmpl$n_obs))

  for (dev in torch_devices()) {
    sol <- gpu_solve_varying_x(tmpl$X, col, Pc, tmpl$y, device = dev)

    for (j in seq_len(B)) {
      Xp <- tmpl$X; Xp[, col] <- Pc[, j]
      ref <- qap_ols_solve(Xp, tmpl$y)
      expect_equal(sol$coefficients[, j], unname(ref$coefficients),
                   tolerance = 1e-9, info = paste(dev, j))
      expect_equal(sol$t[, j], unname(ref$t),
                   tolerance = 1e-9, info = paste(dev, j))
    }
  }
})

test_that("CUDA and CPU torch agree with each other", {
  skip_no_torch()
  skip_if_not(gpu_available(), "no CUDA device")
  d <- gt_data(n = 14); tmpl <- gt_tmpl(d)

  set.seed(6)
  Y <- vapply(1:5, function(j) tmpl$y[sample.int(tmpl$n_obs)],
              numeric(tmpl$n_obs))
  a <- gpu_solve_fixed_x(tmpl$X, Y, device = "cpu")
  b <- gpu_solve_fixed_x(tmpl$X, Y, device = "cuda")

  expect_equal(a$coefficients, b$coefficients, tolerance = 1e-9)
  expect_equal(a$t,            b$t,            tolerance = 1e-9)
})


# ---- gpu_batch_ols end to end -------------------------------------------

test_that("gpu_batch_ols returns well-formed results under qapy", {
  skip_no_torch()
  d <- gt_data(); tmpl <- gt_tmpl(d, c("x1", "x2"))
  base <- qap_ols_solve(tmpl$X, tmpl$y)

  set.seed(7)
  g <- gpu_batch_ols(tmpl, d, "y", c("x1", "x2"), groups = NULL, reps = 40,
                     baseline_fit = base, perm_var = NULL, device = "cpu")

  expect_equal(dim(g$abs), c(2L, 3L))
  expect_equal(rownames(g$abs), qap_stat_rows())
  expect_equal(colnames(g$abs), colnames(tmpl$X))
  expect_true(all(g$abs >= 0 & g$abs <= 1))
  expect_equal(g$n_valid, 40L)
  expect_equal(dim(g$draws$b), c(40L, 3L))

  # The reported proportions must be exactly what the retained draws say.
  for (v in colnames(tmpl$X))
    expect_equal(unname(g$abs["t", v]),
                 mean(abs(g$draws$t[, v]) >= abs(base$t[[v]])),
                 tolerance = 1e-12, info = v)
})

test_that("gpu_batch_ols returns well-formed results under qapspp", {
  skip_no_torch()
  d <- gt_data(); tmpl <- gt_tmpl(d, c("x1", "x2"))
  base <- qap_ols_solve(tmpl$X, tmpl$y)

  set.seed(8)
  g <- gpu_batch_ols(tmpl, d, "y", c("x1", "x2"), groups = NULL, reps = 30,
                     baseline_fit = base, perm_var = "x1", device = "cpu")
  expect_true(all(is.finite(g$draws$b[, "x1"])))
  expect_equal(g$n_valid, 30L)
})

test_that("batches that do not divide the rep count are handled", {
  skip_no_torch()
  d <- gt_data(); tmpl <- gt_tmpl(d, c("x1", "x2"))
  base <- qap_ols_solve(tmpl$X, tmpl$y)

  # 17 reps in batches of 5: three full batches and a partial one.
  set.seed(9)
  g <- gpu_batch_ols(tmpl, d, "y", c("x1", "x2"), groups = NULL, reps = 17,
                     baseline_fit = base, perm_var = NULL, batch_size = 5,
                     device = "cpu")
  expect_equal(g$n_valid, 17L)
  expect_equal(nrow(g$draws$b), 17L)
  expect_false(anyNA(g$draws$b))
})

test_that("the runtime self-check passes on a correct solver", {
  skip_no_torch()
  # gpu_batch_ols() re-solves its first batch on the CPU and stops if the
  # two disagree by more than 1e-6. Reaching a result at all means that
  # check ran and passed.
  d <- gt_data(); tmpl <- gt_tmpl(d, c("x1", "x2"))
  base <- qap_ols_solve(tmpl$X, tmpl$y)
  set.seed(10)
  expect_no_error(
    gpu_batch_ols(tmpl, d, "y", c("x1", "x2"), groups = NULL, reps = 6,
                  baseline_fit = base, perm_var = NULL, device = "cpu"))
})


# ---- through QAP() ------------------------------------------------------

test_that("use_gpu produces the same inference as the CPU, up to Monte Carlo", {
  skip_no_torch()
  # The GPU draws its permutations in batches and the CPU one at a time, so
  # the two see different permutations and cannot agree exactly. What must
  # agree is the baseline fit, exactly, and the p-values, to within the
  # sampling error of 400 draws.
  d <- gt_data(n = 16)

  for (nh in c("qapy", "qapspp")) {
    cpu <- QAP(y ~ x1 + x2, data = d, nullhyp = nh, reps = 400, seed = 1)
    gpu <- QAP(y ~ x1 + x2, data = d, nullhyp = nh, reps = 400, seed = 1,
               use_gpu = TRUE)

    expect_equal(gpu$coefficients, cpu$coefficients, tolerance = 1e-10,
                 info = nh)
    expect_equal(gpu$t, cpu$t, tolerance = 1e-10, info = nh)
    for (v in c("x1", "x2"))
      expect_equal(unname(gpu$abs["t", v]), unname(cpu$abs["t", v]),
                   tolerance = 0.08, info = paste(nh, v))
  }
})

test_that("a GPU fit supports the same extractors as any other", {
  skip_no_torch()
  d <- gt_data(n = 14)
  f <- QAP(y ~ x1 + x2, data = d, reps = 100, seed = 2, use_gpu = TRUE)

  expect_s3_class(f, "QAPRegression")
  expect_equal(nrow(summary(f)$coefficients), 3L)
  expect_true(all(is.finite(confint(f)[c("x1", "x2"), ])))
  expect_equal(nobs(f), 14 * 13)
  expect_length(predict(f), 14 * 13)
})

test_that("multi-network and grouped data work on the GPU", {
  skip_no_torch()
  set.seed(11); n <- 9
  mk <- function() matrix(rnorm(n^2), n, n)
  x1 <- list(mk(), mk())
  y  <- lapply(x1, function(a) 0.5 * a + matrix(rnorm(n^2), n, n))
  d  <- list(y = y, x1 = x1)

  f <- QAP(y ~ x1, data = d, reps = 50, seed = 3, use_gpu = TRUE)
  expect_equal(nobs(f), 2 * n * (n - 1))
  expect_true(all(is.finite(f$abs["t", ])))

  g <- QAP(y ~ x1, data = gt_data(n = 10), reps = 50, seed = 3,
           groups = rep(1:2, each = 5), use_gpu = TRUE)
  expect_true(all(is.finite(g$abs["t", ])))
})

test_that("CSS arrays work on the GPU", {
  skip_no_torch()
  set.seed(12); n <- 6
  x1 <- array(rnorm(n^3), c(n, n, n))
  d  <- list(y = 0.5 * x1 + array(rnorm(n^3), c(n, n, n)), x1 = x1)

  f <- QAP(y ~ x1, data = d, reps = 40, seed = 4, use_gpu = TRUE)
  expect_s3_class(f, "QAPCSS")
  expect_true(is.finite(f$abs["t", "x1"]))
})


# ---- what the GPU path refuses ------------------------------------------

test_that("gpu_batch_ols refuses robust standard errors", {
  skip_no_torch()
  d <- gt_data(); tmpl <- gt_tmpl(d, c("x1", "x2"))
  base <- qap_ols_solve(tmpl$X, tmpl$y)
  expect_error(
    gpu_batch_ols(tmpl, d, "y", c("x1", "x2"), groups = NULL, reps = 5,
                  baseline_fit = base, use_robust_errors = TRUE,
                  device = "cpu"),
    "does not support use_robust_errors")
})

test_that("ineligible models warn and fall back to the CPU", {
  skip_no_torch()
  d <- gt_data(n = 10)

  # Robust errors, a non-gaussian family, and NAs each rule the GPU out.
  expect_warning(QAP(y ~ x1 + x2, data = d, use_gpu = TRUE,
                     use_robust_errors = TRUE, reps = 5, seed = 1),
                 "cannot use the GPU")

  dn <- d; dn$y[2, 3] <- NA
  expect_warning(QAP(y ~ x1 + x2, data = dn, use_gpu = TRUE, reps = 5,
                     seed = 1),
                 "cannot use the GPU")
})

test_that("gpu_available reports the real state of this machine", {
  expect_type(gpu_available(), "logical")
  if (!requireNamespace("torch", quietly = TRUE))
    expect_false(gpu_available())
})
