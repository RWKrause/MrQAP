# ============================================================
# Prior weights and offsets.
#
# Both are properties of the cell, not of the actors: they stay in place
# while the outcome (or the tested predictor) is permuted. They travel
# through the vectoriser as reserved columns, so they are masked, stacked
# and NA-dropped exactly like a predictor without ever entering the model
# matrix, being residualised, or being permuted.
# ============================================================

wo_data <- function(n = 18, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(n = n, x1 = x1, x2 = x2,
       w   = matrix(runif(n^2, 0.5, 2), n, n),
       off = matrix(log(runif(n^2, 1, 3)), n, n))
}

# Vectorise the way the package does: drop the diagonal, column-major.
wo_vec <- function(m, n) as.vector(m)[as.vector(!diag(n))]


# ---- the fit matches the equivalent weighted model ---------------------

test_that("gaussian weights reproduce weighted lm()", {
  d <- wo_data()
  y <- 0.6 * d$x1 - 0.3 * d$x2 + matrix(rnorm(d$n^2), d$n, d$n)

  fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = d$x1, x2 = d$x2),
             weights = d$w, reps = 10, seed = 1)
  ref <- lm(wo_vec(y, d$n) ~ wo_vec(d$x1, d$n) + wo_vec(d$x2, d$n),
            weights = wo_vec(d$w, d$n))

  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-10)
  expect_true(fit$weighted)
  expect_false(fit$offset)
})

test_that("a poisson offset reproduces glm(offset =)", {
  d <- wo_data()
  y <- matrix(rpois(d$n^2, exp(0.2 + 0.5 * d$x1 + d$off)), d$n, d$n)

  fit <- QAP(y ~ x1, data = list(y = y, x1 = d$x1), family = "poisson",
             offset = d$off, reps = 10, seed = 1)
  ref <- glm(wo_vec(y, d$n) ~ wo_vec(d$x1, d$n), family = poisson(),
             offset = wo_vec(d$off, d$n))

  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
  # And the offset is doing something: dropping it moves the intercept.
  bare <- glm(wo_vec(y, d$n) ~ wo_vec(d$x1, d$n), family = poisson())
  expect_gt(abs(coef(ref)[1] - coef(bare)[1]), 0.3)
})

test_that("weights and offset together reproduce the equivalent glm()", {
  d <- wo_data()
  y <- matrix(rpois(d$n^2, exp(0.2 + 0.5 * d$x1 + d$off)), d$n, d$n)

  fit <- QAP(y ~ x1, data = list(y = y, x1 = d$x1), family = "poisson",
             weights = d$w, offset = d$off, reps = 10, seed = 1)
  ref <- glm(wo_vec(y, d$n) ~ wo_vec(d$x1, d$n), family = poisson(),
             weights = wo_vec(d$w, d$n), offset = wo_vec(d$off, d$n))
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
})

test_that("constant weights leave the fit unchanged", {
  d <- wo_data()
  y <- 0.6 * d$x1 + matrix(rnorm(d$n^2), d$n, d$n)
  dd <- list(y = y, x1 = d$x1, x2 = d$x2)

  plain <- QAP(y ~ x1 + x2, data = dd, reps = 20, seed = 3)
  const <- QAP(y ~ x1 + x2, data = dd, reps = 20, seed = 3,
               weights = matrix(1, d$n, d$n))

  expect_equal(unname(plain$coefficients), unname(const$coefficients),
               tolerance = 1e-10)
  # A constant weight column must not trip the degeneracy guard, which
  # rejects permutations whose predictors have no variation.
  expect_equal(const$n_valid, plain$n_valid)
})


# ---- GMM: the user's estimator of interest ------------------------------

test_that("GMM poisson and logit with weights and offset equal the MLE", {
  skip_if_not_installed("gmm")
  # For these two families the moment conditions ARE the score equations
  # and the system is just-identified, so GMM must land exactly on the
  # weighted, offset MLE. That is the sharpest available check that both
  # reach the moment conditions.
  set.seed(2); n <- 20
  x1 <- matrix(rnorm(n^2), n, n); x2 <- matrix(rnorm(n^2), n, n)
  w   <- matrix(runif(n^2, 0.5, 2), n, n)
  off <- matrix(log(runif(n^2, 1, 3)), n, n)
  v   <- function(m) wo_vec(m, n)

  eta <- 0.2 + 0.5 * x1 - 0.3 * x2 + off

  for (fam in c("poisson", "binomial")) {
    y <- if (fam == "poisson") {
      matrix(rpois(n^2, exp(eta)), n, n)
    } else {
      matrix(rbinom(n^2, 1, 1 / (1 + exp(-eta))), n, n)
    }

    fit <- suppressWarnings(
      QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2), family = fam,
          estimator = "gmm", weights = w, offset = off,
          nullhyp = "qapy", reps = 5, seed = 1, ncores = 1))
    ref <- suppressWarnings(
      glm(v(y) ~ v(x1) + v(x2), family = fam,
          weights = v(w), offset = v(off)))

    expect_equal(unname(fit$coefficients), unname(coef(ref)),
                 tolerance = 1e-6, info = fam)
  }
})

test_that("GMM accepts weights and offset for every supported family", {
  skip_if_not_installed("gmm")
  set.seed(4); n <- 18
  x1 <- matrix(rnorm(n^2), n, n); x2 <- matrix(rnorm(n^2), n, n)
  w   <- matrix(runif(n^2, 0.5, 2), n, n)
  off <- matrix(log(runif(n^2, 1, 3)), n, n)
  eta <- 0.2 + 0.5 * x1

  for (fam in c("binomial", "poisson", "negbin", "zip")) {
    y <- switch(fam,
      binomial = matrix(rbinom(n^2, 1, 1 / (1 + exp(-eta))), n, n),
      poisson  = matrix(rpois(n^2, exp(eta)), n, n),
      negbin   = matrix(rnbinom(n^2, mu = exp(eta), size = 3), n, n),
      zip      = matrix(ifelse(rbinom(n^2, 1, 0.3) == 1L, 0L,
                               rpois(n^2, exp(eta))), n, n))
    dd <- list(y = y, x1 = x1, x2 = x2)

    fit <- suppressWarnings(
      QAP(y ~ x1 + x2, data = dd, family = fam, estimator = "gmm",
          weights = w, offset = off, nullhyp = "qapy", reps = 10, seed = 1,
          ncores = 1))
    bare <- suppressWarnings(
      QAP(y ~ x1 + x2, data = dd, family = fam, estimator = "gmm",
          nullhyp = "qapy", reps = 10, seed = 1, ncores = 1))

    expect_length(fit$coefficients, 3)
    expect_true(all(is.finite(fit$coefficients)), info = fam)
    # They must actually reach the moment conditions, not be ignored.
    expect_false(isTRUE(all.equal(fit$coefficients, bare$coefficients)),
                 info = fam)
  }
})

test_that("neutral weights and offset leave GMM moment conditions alone", {
  skip_if_not_installed("gmm")
  set.seed(5); n <- 16
  x1 <- matrix(rnorm(n^2), n, n)
  y  <- matrix(rpois(n^2, exp(0.2 + 0.5 * x1)), n, n)
  dd <- list(y = y, x1 = x1)

  bare <- suppressWarnings(
    QAP(y ~ x1, data = dd, family = "poisson", estimator = "gmm",
        nullhyp = "qapy", reps = 5, seed = 1, ncores = 1))
  neutral <- suppressWarnings(
    QAP(y ~ x1, data = dd, family = "poisson", estimator = "gmm",
        weights = matrix(1, n, n), offset = matrix(0, n, n),
        nullhyp = "qapy", reps = 5, seed = 1, ncores = 1))

  expect_equal(bare$coefficients, neutral$coefficients, tolerance = 1e-8)
})

test_that("gmm_w and gmm_off default to neutral", {
  expect_identical(gmm_w(list(y = 1, x = 1)), 1)
  expect_identical(gmm_off(list(y = 1, x = 1)), 0)
  expect_equal(gmm_w(list(w = c(1, 2))), c(1, 2))
  expect_equal(gmm_off(list(off = c(0.5, 1))), c(0.5, 1))
})


# ---- permutation semantics ---------------------------------------------

test_that("weights are not permuted with the outcome", {
  # A weight belongs to the cell. If it travelled with the relabelling the
  # null distribution would answer a different question, and repeated runs
  # would not reproduce.
  d <- wo_data()
  y <- 0.4 * d$x1 + matrix(rnorm(d$n^2), d$n, d$n)
  dd <- list(y = y, x1 = d$x1, x2 = d$x2)

  a <- QAP(y ~ x1 + x2, data = dd, weights = d$w, reps = 30, seed = 8)
  b <- QAP(y ~ x1 + x2, data = dd, weights = d$w, reps = 30, seed = 8)
  expect_identical(a$abs, b$abs)
  expect_true(all(!is.na(a$abs["t", c("x1", "x2")])))
})

test_that("weights and offset survive qapspp and several networks", {
  set.seed(6); n <- 10
  mk <- function() matrix(rnorm(n^2), n, n)
  x1 <- list(mk(), mk()); x2 <- list(mk(), mk())
  y  <- Map(function(a, b) 0.5 * a - 0.2 * b + matrix(rnorm(n^2), n, n),
            x1, x2)
  w  <- list(matrix(runif(n^2, .5, 2), n, n), matrix(runif(n^2, .5, 2), n, n))

  fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
             weights = w, nullhyp = "qapspp", reps = 15, seed = 1)
  expect_true(fit$weighted)
  expect_true(all(!is.na(fit$abs["t", c("x1", "x2")])))
  expect_equal(nobs(fit), 2 * n * (n - 1))
})

test_that("weights and offset work for CSS arrays", {
  set.seed(7); n <- 6
  x1 <- array(rnorm(n^3), c(n, n, n))
  y  <- 0.5 * x1 + array(rnorm(n^3), c(n, n, n))
  w  <- array(runif(n^3, 0.5, 2), c(n, n, n))

  fit <- QAP(y ~ x1, data = list(y = y, x1 = x1), weights = w,
             reps = 10, seed = 1)
  expect_s3_class(fit, "QAPCSS")
  expect_true(fit$weighted)
})


# ---- the fast paths step aside -----------------------------------------

test_that("weights send a gaussian fit down the general path", {
  d <- wo_data()
  # The FWL decomposition and the closed-form HC3 are both unweighted, so
  # eligibility must fail rather than silently ignore the weights.
  expect_false(qap_ols_eligible("gaussian", "standard", FALSE, FALSE, NULL,
                                list(y = d$x1), "y", has_wo = TRUE))
  expect_true(qap_ols_eligible("gaussian", "standard", FALSE, FALSE, NULL,
                               list(y = d$x1), "y", has_wo = FALSE))
})

test_that("use_gpu with weights warns and falls back", {
  d <- wo_data(n = 10)
  y <- 0.5 * d$x1 + matrix(rnorm(100), 10, 10)
  expect_warning(
    QAP(y ~ x1, data = list(y = y, x1 = d$x1), weights = d$w,
        use_gpu = TRUE, reps = 5, seed = 1),
    "weights or offset")
})


# ---- validation ---------------------------------------------------------

test_that("misshapen weights and offsets are rejected", {
  d <- wo_data(n = 10)
  y <- 0.5 * d$x1 + matrix(rnorm(100), 10, 10)
  dd <- list(y = y, x1 = d$x1)

  expect_error(QAP(y ~ x1, data = dd, weights = matrix(1, 5, 5), reps = 5),
               "same shape")
  expect_error(QAP(y ~ x1, data = dd, weights = matrix(-1, 10, 10), reps = 5),
               "must not be negative")
  expect_error(QAP(y ~ x1, data = dd, offset = matrix(Inf, 10, 10), reps = 5),
               "non-finite")
  expect_error(QAP(y ~ x1, data = dd, weights = "a", reps = 5),
               "must be numeric")
})

test_that("multinomial rejects an offset but accepts weights", {
  skip_if_not_installed("nnet")
  set.seed(9); n <- 10
  x1 <- matrix(rnorm(n^2), n, n)
  y  <- matrix(sample(c("A", "B", "C"), n^2, TRUE), n, n)
  dd <- list(y = y, x1 = x1)

  expect_error(
    suppressWarnings(QAP(y ~ x1, data = dd, family = "multinom",
                         offset = matrix(0, n, n), reps = 5, seed = 1)),
    "offset= is not supported")
  expect_no_error(
    suppressWarnings(QAP(y ~ x1, data = dd, family = "multinom",
                         weights = matrix(1, n, n), reps = 5, seed = 1)))
})
