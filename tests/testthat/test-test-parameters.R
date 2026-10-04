# ============================================================
# test_parameters: permute only the predictors a hypothesis concerns.
#
# Predictors are tested in formula order, each in its own stream of draws,
# so testing a prefix of them reproduces the full run's results for that
# prefix exactly -- which is what most of these tests lean on.
# ============================================================

tp_data <- function(n = 10, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  x3 <- matrix(rnorm(n^2), n, n)
  list(y = 0.4 * x1 - 0.2 * x2 + matrix(rnorm(n^2), n, n),
       x1 = x1, x2 = x2, x3 = x3)
}

tp_binary <- function(n = 12, seed = 2) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  x3 <- matrix(rnorm(n^2), n, n)
  eta <- 0.8 * x1
  list(y = matrix(rbinom(n^2, 1, 1 / (1 + exp(-eta))), n, n),
       x1 = x1, x2 = x2, x3 = x3)
}

# The fitted model objects capture the calling frame, test_parameters
# included, so compare everything else.
tp_strip <- function(f) {
  f$simple_fit <- NULL
  f$base$base_model <- NULL
  f
}

# The subset run must agree with the full run on x1 and leave the rest NA.
expect_prefix_matches <- function(sub, full) {
  for (m in c("lower", "larger", "abs"))
    expect_equal(sub[[m]][, "x1"], full[[m]][, "x1"], info = m)
  expect_equal(sub$null_dist$t[, "x1"], full$null_dist$t[, "x1"])
  expect_equal(sub$null_moments[, "x1"], full$null_moments[, "x1"])
  for (v in c("x2", "x3")) {
    expect_true(all(is.na(sub$abs[, v])), info = v)
    expect_true(all(is.na(sub$null_moments[, v])), info = v)
  }
  expect_identical(sub$coefficients, full$coefficients)
  expect_identical(sub$test_parameters, "x1")
}


# ---- the default is unchanged ---------------------------------------------

test_that("NULL and naming every predictor give identical fits", {
  d <- tp_data()
  a <- QAP(y ~ x1 + x2 + x3, data = d, reps = 20, seed = 1)
  b <- QAP(y ~ x1 + x2 + x3, data = d, reps = 20, seed = 1,
           test_parameters = c("x1", "x2", "x3"))
  expect_identical(tp_strip(a), tp_strip(b))
  expect_null(a$test_parameters)
})


# ---- a subset, on each path ---------------------------------------------

test_that("a subset reproduces the full run on the fast OLS path", {
  d <- tp_data()
  full <- QAP(y ~ x1 + x2 + x3, data = d, reps = 25, seed = 1)
  sub  <- QAP(y ~ x1 + x2 + x3, data = d, reps = 25, seed = 1,
              test_parameters = "x1")
  expect_prefix_matches(sub, full)
  expect_equal(names(sub$n_valid), "x1")
})

test_that("a subset reproduces the full run on the general CPU path", {
  d <- tp_binary()
  full <- QAP(y ~ x1 + x2 + x3, data = d, family = "binomial",
              reps = 15, seed = 1)
  sub  <- QAP(y ~ x1 + x2 + x3, data = d, family = "binomial",
              reps = 15, seed = 1, test_parameters = "x1")
  expect_prefix_matches(sub, full)
})

test_that("a subset reproduces the full run on the GPU path", {
  skip_if_not(torch_ready(), "torch has no usable LibTorch back end")
  d <- tp_data()
  full <- suppressMessages(
    QAP(y ~ x1 + x2 + x3, data = d, reps = 25, seed = 1, use_gpu = TRUE))
  sub  <- suppressMessages(
    QAP(y ~ x1 + x2 + x3, data = d, reps = 25, seed = 1, use_gpu = TRUE,
        test_parameters = "x1"))
  expect_prefix_matches(sub, full)
})

test_that("the order the names are given in does not matter", {
  d <- tp_data()
  a <- QAP(y ~ x1 + x2 + x3, data = d, reps = 20, seed = 1,
           test_parameters = c("x3", "x1"))
  b <- QAP(y ~ x1 + x2 + x3, data = d, reps = 20, seed = 1,
           test_parameters = c("x1", "x3"))
  expect_identical(tp_strip(a), tp_strip(b))
  expect_identical(a$test_parameters, c("x1", "x3"))
  expect_true(all(is.na(a$abs[, "x2"])))
  expect_false(anyNA(a$abs[, c("x1", "x3")]))
})

test_that("a later predictor alone is still tested correctly", {
  d <- tp_data()
  f <- QAP(y ~ x1 + x2 + x3, data = d, reps = 25, seed = 1,
           test_parameters = "x3")
  obs <- f$t[["x3"]]
  expect_equal(f$abs["t", "x3"], mean(abs(f$null_dist$t[, "x3"]) >= abs(obs)))
  expect_true(all(is.na(f$abs[, c("x1", "x2")])))
})


# ---- validation ------------------------------------------------------------

test_that("unknown names, the intercept and bad types are rejected", {
  d <- tp_data()
  fit <- function(tp) QAP(y ~ x1 + x2 + x3, data = d, reps = 5, seed = 1,
                          test_parameters = tp)
  expect_error(fit("x9"), "'x9'.*Choose from: 'x1', 'x2', 'x3'")
  expect_error(fit("(Intercept)"), "never permuted on their own")
  expect_error(fit(1), "character vector")
  expect_error(fit(character(0)), "character vector")
  expect_error(fit(NA_character_), "character vector")
})

test_that("interaction terms cannot be named", {
  d <- tp_data()
  expect_error(
    QAP(y ~ x1 * x2, data = d, reps = 5, seed = 1,
        test_parameters = "x1:x2"),
    "interaction terms")
})

test_that("qapy with a strict subset warns and tests everything", {
  d <- tp_data()
  expect_warning(
    f <- QAP(y ~ x1 + x2 + x3, data = d, reps = 10, seed = 1,
             nullhyp = "qapy", test_parameters = "x1"),
    "only reduces the work under nullhyp")
  ref <- QAP(y ~ x1 + x2 + x3, data = d, reps = 10, seed = 1,
             nullhyp = "qapy")
  expect_identical(tp_strip(f), tp_strip(ref))
})

test_that("a single predictor (forced qapy) does not warn", {
  d <- tp_data()
  expect_no_warning(
    QAP(y ~ x1, data = d, reps = 5, seed = 1, test_parameters = "x1"))
})


# ---- reporting and downstream ----------------------------------------------

test_that("print and summary say which predictors were tested", {
  d <- tp_data()
  f <- QAP(y ~ x1 + x2 + x3, data = d, reps = 10, seed = 1,
           test_parameters = "x1")
  expect_true(any(grepl("Only x1 was tested", capture.output(print(f)))))
  expect_true(any(grepl("Only x1 was tested",
                        capture.output(print(summary(f))))))
  tab <- summary(f)$coefficients
  expect_true(all(is.na(tab$p_value[tab$term %in% c("x2", "x3")])))

  g <- QAP(y ~ x1 + x2 + x3, data = tp_binary(), family = "binomial",
           reps = 5, seed = 1, test_parameters = c("x1", "x2"))
  expect_true(any(grepl("Only x1, x2 were tested", capture.output(print(g)))))
})

test_that("plot() names test_parameters for an untested predictor", {
  f <- QAP(y ~ x1 + x2 + x3, data = tp_data(), reps = 10, seed = 1,
           test_parameters = "x1")
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plot(f, which = "x2"), "test_parameters")
  expect_silent(plot(f))
})

test_that("runs with different test_parameters cannot be pooled", {
  d <- tp_data()
  a <- QAP(y ~ x1 + x2 + x3, data = d, reps = 10, seed = 1,
           test_parameters = "x1")
  b <- QAP(y ~ x1 + x2 + x3, data = d, reps = 10, seed = 2)
  expect_error(combine_qap_estimates(a, b), "same test_parameters")

  c2 <- QAP(y ~ x1 + x2 + x3, data = d, reps = 10, seed = 2,
            test_parameters = "x1")
  ac <- combine_qap_estimates(a, c2)
  expect_equal(ac$n_valid, c(x1 = 20L))
  expect_true(all(is.na(ac$abs[, "x2"])))
})
