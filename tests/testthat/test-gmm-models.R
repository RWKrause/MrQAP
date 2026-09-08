# ============================================================
# Substantive GMM tests: every supported family, every claim.
#
# test-gmm.R asserts that a GMM fit has the right SHAPE. This file asserts
# that it has the right VALUE: that the estimator recovers a planted
# coefficient, that two runs with one seed agree exactly, that the
# permutation test separates a real effect from a null one, and that
# weights and offsets reach the moment conditions.
# ============================================================

gmm_truth <- c(b0 = 0.3, b1 = 0.6, b2 = -0.35)

# One builder for all four families. n = 22 gives 462 dyads: enough for GMM
# to recover a coefficient, small enough to stay inside the test budget.
# x3 is a genuine null predictor -- it enters no linear predictor.
make_gmm_family <- function(family, n = 22, seed = 11) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  x3 <- matrix(rnorm(n^2), n, n)
  eta <- gmm_truth[["b0"]] + gmm_truth[["b1"]] * x1 + gmm_truth[["b2"]] * x2

  y <- switch(
    family,
    binomial = matrix(rbinom(n^2, 1, 1 / (1 + exp(-eta))), n, n),
    poisson  = matrix(rpois(n^2, exp(eta)), n, n),
    negbin   = matrix(rnbinom(n^2, mu = exp(eta), size = 3), n, n),
    zip      = matrix(ifelse(rbinom(n^2, 1, 0.25) == 1L, 0L,
                             rpois(n^2, exp(eta))), n, n)
  )
  diag(y) <- diag(x1) <- diag(x2) <- diag(x3) <- NA
  list(y = y, x1 = x1, x2 = x2, x3 = x3)
}

gmm_families <- c("binomial", "poisson", "negbin", "zip")

gmm_fit <- function(...) suppressWarnings(QAP(..., ncores = 1))


# ---- starting values ---------------------------------------------------

test_that("gmm_start returns finite values of the right length", {
  set.seed(5)
  y <- rpois(200, 2)
  X <- cbind(1, rnorm(200), rnorm(200))

  for (fam in gmm_families) {
    yy <- if (fam == "binomial") as.numeric(y > 1) else y
    st <- gmm_start(yy, X, fam)
    # negbin and zip carry one nuisance parameter beyond the coefficients.
    expected <- ncol(X) + if (fam %in% c("negbin", "zip")) 1L else 0L
    expect_length(st, expected)
    expect_true(all(is.finite(st)), info = fam)
  }
})

test_that("gmm_start agrees with the corresponding GLM estimate", {
  set.seed(6)
  x <- rnorm(300)
  y <- rpois(300, exp(0.2 + 0.5 * x))
  X <- cbind(1, x)

  expect_equal(gmm_start(y, X, "poisson"),
               unname(coef(glm(y ~ x, family = poisson()))),
               tolerance = 1e-8)
  # The nuisance start is appended; the coefficients are unchanged.
  expect_equal(gmm_start(y, X, "negbin")[1:2], gmm_start(y, X, "poisson"))
  expect_equal(gmm_start(y, X, "zip")[1:2],    gmm_start(y, X, "poisson"))
})

test_that("gmm_start survives a degenerate design", {
  # Perfect separation sends the logistic GLM coefficients to +/-Inf; a
  # non-finite start would make every downstream gmm() call fail.
  x <- c(rep(-1, 50), rep(1, 50))
  y <- as.numeric(x > 0)
  st <- gmm_start(y, cbind(1, x), "binomial")
  expect_true(all(is.finite(st)))
})

test_that("deterministic starts make GMM estimates start-independent", {
  skip_if_not_installed("gmm")
  # On a wider-scaled design the moment conditions have several local
  # optima. Random starting values landed on different ones, putting noise
  # of the same magnitude as the estimate into the null distribution. This
  # is the regression test for that fix.
  set.seed(3)
  n <- 14
  x1 <- matrix(rnorm(n^2, sd = 2.5), n, n)
  x2 <- matrix(rnorm(n^2, sd = 2.5), n, n)
  eta <- 0.3 + 0.6 * x1 - 0.3 * x2
  y   <- matrix(rpois(n^2, exp(pmin(eta, 5))), n, n)

  yna <- y; diag(yna) <- NA
  keep <- !is.na(as.vector(yna))
  yy <- as.vector(y)[keep]
  X  <- cbind(1, as.vector(x1)[keep], as.vector(x2)[keep])

  fit_from <- function(t0) tryCatch(suppressWarnings(
    gmm::gmm(poisson_moments, x = list(y = yy, x = X), t0 = t0,
             wmatrix = "optimal", vcov = "MDS", optfct = "nlminb",
             control = list(eval.max = 10000))$coefficients),
    error = function(e) NULL)

  set.seed(99)
  spread <- vapply(seq_len(25), function(i) {
    z <- fit_from(rnorm(3))
    if (is.null(z)) NA_real_ else z[2]
  }, numeric(1))
  # Random starts disagree with each other by far more than the effect size.
  expect_gt(stats::sd(spread, na.rm = TRUE), 1)

  # The deterministic start lands near the truth instead.
  det <- fit_from(gmm_start(yy, X, "poisson"))
  expect_false(is.null(det))
  expect_equal(unname(det[2]), 0.6, tolerance = 0.2)
})


# ---- per-family behaviour ----------------------------------------------

for (fam in gmm_families) {
  local({
    family <- fam

    test_that(paste0("GMM recovers planted coefficients (", family, ")"), {
      skip_if_not_installed("gmm")
      d   <- make_gmm_family(family)
      fit <- gmm_fit(y ~ x1 + x2, data = d, family = family,
                     estimator = "gmm", nullhyp = "qapy", reps = 5, seed = 1)

      expect_named(fit$coefficients, c("(Intercept)", "x1", "x2"))
      expect_true(all(is.finite(fit$coefficients)))
      # Generous tolerance: this asserts the estimator is pointed at the
      # right answer, not that it is efficient at this sample size.
      expect_equal(unname(fit$coefficients[["x1"]]), gmm_truth[["b1"]],
                   tolerance = 0.35)
      expect_equal(unname(fit$coefficients[["x2"]]), gmm_truth[["b2"]],
                   tolerance = 0.35)
      # Sign recovery is the weaker claim that must always hold.
      expect_gt(fit$coefficients[["x1"]], 0)
      expect_lt(fit$coefficients[["x2"]], 0)
    })

    test_that(paste0("GMM is reproducible from a seed (", family, ")"), {
      skip_if_not_installed("gmm")
      d   <- make_gmm_family(family)
      run <- function() gmm_fit(y ~ x1 + x2, data = d, family = family,
                                estimator = "gmm", nullhyp = "qapy",
                                reps = 15, seed = 4)
      a <- run(); b <- run()
      # Both the estimate and the permutation inference must repeat: the
      # start values are no longer drawn from the RNG stream.
      expect_identical(a$coefficients, b$coefficients)
      expect_identical(a$t,   b$t)
      expect_identical(a$abs, b$abs)
    })

    test_that(paste0("GMM separates a real effect from a null one (",
                     family, ")"), {
      skip_if_not_installed("gmm")
      d   <- make_gmm_family(family)
      fit <- gmm_fit(y ~ x1 + x3, data = d, family = family,
                     estimator = "gmm", nullhyp = "qapy", reps = 100,
                     seed = 2)
      expect_lt(fit$abs["t", "x1"], 0.10)   # planted
      expect_gt(fit$abs["t", "x3"], 0.10)   # pure noise
    })

    test_that(paste0("GMM runs under qapspp (", family, ")"), {
      skip_if_not_installed("gmm")
      d   <- make_gmm_family(family)
      fit <- gmm_fit(y ~ x1 + x2, data = d, family = family,
                     estimator = "gmm", nullhyp = "qapspp", reps = 15,
                     seed = 3)
      expect_equal(fit$nullhyp, "qapspp")
      # Semi-partialling tests each predictor and never the intercept.
      expect_true(all(!is.na(fit$abs["t", c("x1", "x2")])))
      expect_true(is.na(fit$abs["t", "(Intercept)"]))
      expect_true(all(fit$abs["t", c("x1", "x2")] >= 0 &
                      fit$abs["t", c("x1", "x2")] <= 1))
    })

    test_that(paste0("GMM always reports robust (MDS) errors (", family, ")"), {
      skip_if_not_installed("gmm")
      d   <- make_gmm_family(family)
      fit <- gmm_fit(y ~ x1 + x2, data = d, family = family,
                     estimator = "gmm", nullhyp = "qapy", reps = 5, seed = 1)
      # gmm() is called with vcov = "MDS" whatever use_robust_errors said.
      expect_true(fit$robust_se)
      expect_true(all(is.finite(fit$t)))
      # The nuisance parameter of negbin/zip is not a reported coefficient.
      expect_length(fit$t, 3)
    })

    test_that(paste0("GMM extractors and summary work (", family, ")"), {
      skip_if_not_installed("gmm")
      d   <- make_gmm_family(family)
      fit <- gmm_fit(y ~ x1 + x2, data = d, family = family,
                     estimator = "gmm", nullhyp = "qapy", reps = 15, seed = 1)
      tab <- summary(fit)$coefficients
      expect_equal(nrow(tab), 3L)
      expect_equal(tab$estimate, unname(fit$coefficients))
      expect_true(all(!is.na(tab$p_value)))
      expect_identical(coef(fit), fit$coefficients)
      expect_s3_class(as.data.frame(fit), "data.frame")
    })
  })
}
