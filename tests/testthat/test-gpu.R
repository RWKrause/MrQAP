# ============================================================
# GPU batch OLS: the pure-R scaffolding.
#
# torch is optional, so this file covers only what runs without it --
# eligibility, batch sizing, and the construction of a permuted design
# batch. The torch algebra itself is exercised in test-gpu-torch.R, which
# skips when torch cannot run and uses CUDA when the machine has it.
#
# gpu_batch_ols() additionally checks its own first batch against
# qap_ols_solve() at runtime and stops if they disagree, so a broken
# torch translation fails loudly rather than moving p-values.
# ============================================================

gd <- function(n = 10, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(y = 0.2 * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
}

test_that("gpu_available does not error when torch is absent", {
  expect_type(gpu_available(), "logical")
  expect_length(gpu_available(), 1)
})

test_that("batch size shrinks as the problem grows", {
  expect_equal(gpu_batch_size(100, 3), 1000L)                 # capped
  expect_lt(gpu_batch_size(39800, 4, budget_mb = 64), 1000L)
  expect_gte(gpu_batch_size(39800, 4, budget_mb = 64), 1L)

  # Stays within budget whenever a single rep fits...
  b <- gpu_batch_size(39800, 4, budget_mb = 64)
  expect_lte(b * 8 * 39800 * 4, 64e6)
  # ...and floors at 1 when even one rep does not.
  expect_equal(gpu_batch_size(1e6, 10, budget_mb = 64), 1L)
})

test_that("the batch budget is tunable via an option", {
  old <- options(MrQAP.gpu_batch_mb = 64)
  on.exit(options(old))
  small <- gpu_batch_size(39800, 4)
  options(MrQAP.gpu_batch_mb = 1024)
  big <- gpu_batch_size(39800, 4)
  expect_gt(big, small)
})

test_that("a permuted column batch is a rearrangement of the original", {
  d <- gd()
  tmpl <- qap_ols_template(d, "y", c("x1", "x2"), FALSE, "directed",
                           FALSE, FALSE, FALSE)
  set.seed(2)
  nb <- 5
  Pc <- gpu_perm_columns(nb, obj = d$x1, extract = tmpl$extract,
                         groups = NULL, n_obs = tmpl$n_obs)

  expect_equal(dim(Pc), c(tmpl$n_obs, nb))
  for (j in seq_len(nb))
    expect_setequal(round(Pc[, j], 10), round(tmpl$X[, "x1"], 10))
  # and successive draws differ
  expect_false(isTRUE(all.equal(Pc[, 1], Pc[, 2])))
})

test_that("use_gpu without a usable torch fails with a torch message", {
  skip_if(torch_ready(), "torch works here; this checks the no-torch path")
  d <- gd()
  # gaussian + complete data is GPU-eligible, so this reaches gpu_batch_ols
  expect_error(QAP(y ~ x1 + x2, data = d, use_gpu = TRUE, reps = 5, seed = 1),
               "torch")
})

test_that("use_gpu warns when the model cannot use the GPU at all", {
  d <- gd()
  set.seed(3)
  db <- d
  db$y <- matrix(rbinom(100, 1, 0.5), 10, 10)
  expect_warning(
    QAP(y ~ x1 + x2, data = db, family = "binomial", use_gpu = TRUE,
        reps = 5, seed = 1),
    "cannot use the GPU"
  )
  expect_warning(
    QAP(y ~ x1 + x2, data = d, use_robust_errors = TRUE, use_gpu = TRUE,
        reps = 5, seed = 1),
    "cannot use the GPU"
  )
})

test_that("GPU runs are attempted only for eligible models", {
  skip_if_not(gpu_available(), "no CUDA-capable torch")
  d <- gd()
  fit <- QAP(y ~ x1 + x2, data = d, use_gpu = TRUE, nullhyp = "qapy",
             reps = 20, seed = 1)
  ref <- QAP(y ~ x1 + x2, data = d, use_gpu = FALSE, nullhyp = "qapy",
             reps = 20, seed = 1)
  expect_equal(fit$coefficients, ref$coefficients, tolerance = 1e-8)
})


# ============================================================
# Torch algebra.
#
# These are the ONLY tests that execute the torch calls. They skip
# when torch cannot run -- the package missing, or its LibTorch back end
# never downloaded. Install both and re-run to verify.
#
# Each compares the batched tensor solve against qap_ols_solve(),
# which is itself checked against lm() in test-fastpath.R.
# ============================================================

rand_design <- function(n = 200, p = 3, seed = 1) {
  set.seed(seed)
  X <- cbind("(Intercept)" = 1,
             matrix(rnorm(n * (p - 1)), n, p - 1,
                    dimnames = list(NULL, paste0("x", seq_len(p - 1)))))
  X
}

test_that("gpu_solve_fixed_x matches the CPU solver for every column", {
  skip_if_not(torch_ready(), "torch has no usable LibTorch back end")
  X <- rand_design(n = 250, p = 4)
  set.seed(2)
  Y <- matrix(rnorm(nrow(X) * 6), nrow(X), 6)

  got <- gpu_solve_fixed_x(X, Y, device = "cpu")
  expect_equal(dim(got$coefficients), c(ncol(X), ncol(Y)))
  expect_equal(dim(got$t), c(ncol(X), ncol(Y)))

  for (j in seq_len(ncol(Y))) {
    ref <- qap_ols_solve(X, Y[, j])
    expect_equal(got$coefficients[, j], unname(ref$coefficients),
                 tolerance = 1e-10, info = paste("column", j))
    expect_equal(got$t[, j], unname(ref$t),
                 tolerance = 1e-10, info = paste("column", j))
  }
})

test_that("gpu_solve_varying_x matches the CPU solver for every design", {
  skip_if_not(torch_ready(), "torch has no usable LibTorch back end")
  X <- rand_design(n = 250, p = 4)
  set.seed(3)
  y  <- rnorm(nrow(X))
  nb <- 6
  Pc <- vapply(seq_len(nb), function(j) sample(X[, 2]), numeric(nrow(X)))

  got <- gpu_solve_varying_x(X, 2L, Pc, y, device = "cpu")
  expect_equal(dim(got$coefficients), c(ncol(X), nb))
  expect_equal(dim(got$t), c(ncol(X), nb))

  for (j in seq_len(nb)) {
    Xj <- X; Xj[, 2] <- Pc[, j]
    ref <- qap_ols_solve(Xj, y)
    expect_equal(got$coefficients[, j], unname(ref$coefficients),
                 tolerance = 1e-10, info = paste("design", j))
    expect_equal(got$t[, j], unname(ref$t),
                 tolerance = 1e-10, info = paste("design", j))
  }
})

test_that("both solvers handle a batch of one", {
  skip_if_not(torch_ready(), "torch has no usable LibTorch back end")
  X <- rand_design(n = 120, p = 3)
  set.seed(4)
  y <- rnorm(nrow(X))

  a <- gpu_solve_fixed_x(X, matrix(y, ncol = 1), device = "cpu")
  expect_equal(a$coefficients[, 1], unname(qap_ols_solve(X, y)$coefficients),
               tolerance = 1e-10)

  b <- gpu_solve_varying_x(X, 2L, X[, 2, drop = FALSE], y, device = "cpu")
  expect_equal(b$coefficients[, 1], unname(qap_ols_solve(X, y)$coefficients),
               tolerance = 1e-10)
})

test_that("the CUDA device gives the same answer as the CPU device", {
  skip_if_not(gpu_available(), "no CUDA-capable torch")
  X <- rand_design(n = 400, p = 4)
  set.seed(5)
  Y <- matrix(rnorm(nrow(X) * 8), nrow(X), 8)

  cpu <- gpu_solve_fixed_x(X, Y, device = "cpu")
  gpu <- gpu_solve_fixed_x(X, Y, device = "cuda")
  expect_equal(gpu$coefficients, cpu$coefficients, tolerance = 1e-8)
  expect_equal(gpu$t, cpu$t, tolerance = 1e-8)

  nb <- 8
  Pc <- vapply(seq_len(nb), function(j) sample(X[, 2]), numeric(nrow(X)))
  expect_equal(gpu_solve_varying_x(X, 2L, Pc, Y[, 1], device = "cuda"),
               gpu_solve_varying_x(X, 2L, Pc, Y[, 1], device = "cpu"),
               tolerance = 1e-8)
})

test_that("a full GPU run agrees statistically with the CPU run", {
  skip_if_not(torch_ready(), "torch has no usable LibTorch back end")
  # The GPU path draws permutations in the main session while the CPU path
  # uses future's L'Ecuyer streams, so the two see DIFFERENT permutations
  # and p-values agree only up to Monte Carlo error. Coefficients come from
  # the same baseline fit and must match exactly.
  set.seed(6)
  n  <- 14
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  d  <- list(y = 0.5 * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)

  for (nh in c("qapy", "qapspp")) {
    g <- QAP(y ~ x1 + x2, data = d, use_gpu = TRUE,  nullhyp = nh,
             reps = 400, seed = 1)
    c_ <- QAP(y ~ x1 + x2, data = d, use_gpu = FALSE, nullhyp = nh,
              reps = 400, seed = 1)

    expect_equal(g$coefficients, c_$coefficients, tolerance = 1e-10, info = nh)
    expect_equal(g$t, c_$t, tolerance = 1e-10, info = nh)
    expect_equal(unname(g$abs["t", "x1"]), unname(c_$abs["t", "x1"]),
                 tolerance = 0.08, info = nh)
    expect_equal(unname(g$abs["t", "x2"]), unname(c_$abs["t", "x2"]),
                 tolerance = 0.08, info = nh)
    # both paths must label the rows the same way
    expect_identical(dimnames(g$abs), dimnames(c_$abs), info = nh)
  }
})

test_that("the GPU path handles CSS and two-mode data", {
  skip_if_not(torch_ready(), "torch has no usable LibTorch back end")
  set.seed(7)
  m  <- 7
  a1 <- array(rnorm(m^3), dim = c(m, m, m))
  dc <- list(y = 0.5 * a1 + array(rnorm(m^3), dim = c(m, m, m)), x1 = a1)
  expect_s3_class(QAP(y ~ x1, data = dc, use_gpu = TRUE, reps = 50, seed = 1),
                  "QAPCSS")

  b1 <- matrix(rnorm(9 * 5), 9, 5)
  dt <- list(y = 0.5 * b1 + matrix(rnorm(9 * 5), 9, 5), x1 = b1)
  f <- QAP(y ~ x1, data = dt, use_gpu = TRUE, reps = 50, seed = 1)
  expect_true(f$multi_mode)
})
