# ============================================================
# Statistical correctness tests.
#
# The rest of the suite checks that things run; these check that the
# numbers are right. Each block corresponds to a fixed defect.
# ============================================================

# ---- undirected mode ---------------------------------------------------

test_that("make_qap_data keeps one observation per dyad in graph mode", {
  n <- 6
  y  <- matrix(1, n, n)
  x1 <- matrix(1, n, n)

  directed <- make_qap_data(y = y, x = list(x1 = x1),
                            diag = FALSE, mode = "digraph", net = 1)
  undirected <- make_qap_data(y = y, x = list(x1 = x1),
                              diag = FALSE, mode = "graph", net = 1)

  expect_equal(nrow(directed), n * (n - 1))
  expect_equal(nrow(undirected), n * (n - 1) / 2)

  # Only upper-triangle cells survive.
  expect_true(all(undirected$sv < undirected$rv))
})

test_that("make_qap_data keeps the diagonal once in graph mode", {
  n <- 5
  y <- matrix(1, n, n)
  pred <- make_qap_data(y = y, x = list(x1 = matrix(1, n, n)),
                        diag = TRUE, mode = "graph", net = 1)
  expect_equal(nrow(pred), n * (n + 1) / 2)
})

test_that("make_qap_data rejects an unknown mode", {
  expect_error(
    make_qap_data(y = matrix(1, 3, 3), x = list(x1 = matrix(1, 3, 3)),
                  diag = FALSE, mode = "undirected", net = 1)
  )
})

test_that("QAPglm undirected matches lm() on the upper triangle", {
  set.seed(42)
  n <- 12
  x1 <- matrix(rnorm(n^2), n, n); x1[lower.tri(x1)] <- t(x1)[lower.tri(x1)]
  x2 <- matrix(rnorm(n^2), n, n); x2[lower.tri(x2)] <- t(x2)[lower.tri(x2)]
  y  <- 1 + 2 * x1 - 0.5 * x2 + matrix(rnorm(n^2, sd = 0.3), n, n)
  y[lower.tri(y)] <- t(y)[lower.tri(y)]

  fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
                mode = "undirected", reps = 5, seed = 1)

  ut  <- upper.tri(y)
  ref <- lm(y[ut] ~ x1[ut] + x2[ut])

  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(fit$r.squared, summary(ref)$r.squared, tolerance = 1e-8)
})

test_that("undirected QAPglm uses half the observations of directed", {
  set.seed(7)
  n <- 10
  x1 <- matrix(rnorm(n^2), n, n); x1[lower.tri(x1)] <- t(x1)[lower.tri(x1)]

  nobs <- function(mode) {
    nrow(make_qap_data(y = x1, x = list(x1 = x1), diag = FALSE,
                       mode = mode, net = 1))
  }
  expect_equal(nobs("graph") * 2, nobs("digraph"))
})

test_that("qapspp on undirected data keeps the residualised matrix symmetric", {
  set.seed(3)
  n <- 8
  x1 <- matrix(rnorm(n^2), n, n); x1[lower.tri(x1)] <- t(x1)[lower.tri(x1)]
  pred <- make_qap_data(y = x1, x = list(x1 = x1), diag = FALSE,
                        mode = "graph", net = 1)
  xR  <- rnorm(nrow(pred))
  out <- residuals_to_matrix(xR, x1, pred, large = FALSE, mode = "graph")

  expect_true(isSymmetric(unname(out)))
  # The residuals really did land in the upper triangle.
  expect_equal(unname(out[pred$location]), xR)
})

test_that("directed residual writeback is left untouched", {
  set.seed(4)
  n <- 6
  x1 <- matrix(rnorm(n^2), n, n)
  pred <- make_qap_data(y = x1, x = list(x1 = x1), diag = FALSE,
                        mode = "digraph", net = 1)
  xR  <- rnorm(nrow(pred))
  out <- residuals_to_matrix(xR, x1, pred, large = FALSE, mode = "digraph")
  expect_equal(unname(out[pred$location]), xR)
})


# ---- RMPerm ------------------------------------------------------------

test_that("RMPerm handles permutation groups of size one", {
  set.seed(1)
  n <- 5
  m <- matrix(seq_len(n^2), n, n)
  g <- c(1, 1, 1, 1, 2)   # group 2 has a single member

  for (i in 1:50) {
    p <- RMPerm(m, groups = g)
    expect_equal(dim(p), dim(m))
    expect_setequal(as.vector(p), as.vector(m))
  }
})

test_that("RMPerm never moves a node out of its group", {
  set.seed(2)
  n <- 6
  g <- c(1, 1, 1, 2, 2, 3)
  # Encode group membership in the matrix so it can be read back.
  m <- outer(g, g, function(a, b) a * 10 + b)

  for (i in 1:50) {
    p <- RMPerm(m, groups = g)
    expect_equal(p %/% 10, outer(g, rep(1, n)))
    expect_equal(p %% 10, outer(rep(1, n), g))
  }
})

test_that("RMPerm forwards CSS = TRUE through list recursion", {
  n <- 4
  a <- array(seq_len(n^3), dim = c(n, n, n))

  set.seed(99); direct  <- RMPerm(a, CSS = TRUE)
  set.seed(99); vialist <- RMPerm(list(a), CSS = TRUE)[[1]]
  expect_equal(direct, vialist)

  set.seed(99); css    <- RMPerm(a, CSS = TRUE)
  set.seed(99); noncss <- RMPerm(a, CSS = FALSE)
  expect_false(identical(css, noncss))
})

test_that("RMPerm rejects non-square input and mismatched groups", {
  expect_error(RMPerm(matrix(1, 3, 4)), "square")
  expect_error(RMPerm(matrix(1, 4, 4), groups = c(1, 1)), "length")
})


# ---- HC3 / robust standard errors --------------------------------------

test_that("HC3 matches sandwich::vcovHC for a linear model", {
  skip_if_not_installed("sandwich")
  set.seed(8)
  N <- 200
  x1 <- rnorm(N); x2 <- rnorm(N)
  y  <- 1 + x1 - x2 + rnorm(N, sd = abs(x1))   # heteroskedastic
  ref <- lm(y ~ x1 + x2)

  expect_equal(unname(HC3(cbind(x1 = x1, x2 = x2), residuals(ref))),
               unname(sqrt(diag(sandwich::vcovHC(ref, type = "HC3")))),
               tolerance = 1e-8)
})

test_that("robust SEs for a binomial fit use the GLM sandwich, not HC3", {
  skip_if_not_installed("sandwich")
  set.seed(12)
  n <- 14
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  eta <- -0.5 + x1
  y <- matrix(rbinom(n^2, 1, 1 / (1 + exp(-eta))), n, n)

  fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
                family = "binomial", use_robust_errors = TRUE,
                reps = 5, seed = 1)

  ref  <- fit$simple_fit
  want <- coef(ref) / sqrt(diag(sandwich::vcovHC(ref, type = "HC3")))

  expect_equal(unname(fit$t[names(want)]), unname(want), tolerance = 1e-8)
})

test_that("gaussian robust SEs still use the closed-form HC3", {
  set.seed(9)
  n <- 14
  x1 <- matrix(rnorm(n^2), n, n)
  y  <- x1 + matrix(rnorm(n^2), n, n)

  fit <- QAP(y ~ x1, data = list(y = y, x1 = x1),
                use_robust_errors = TRUE, reps = 5, seed = 1)

  ref  <- fit$simple_fit
  want <- coef(ref) / HC3(as.matrix(model.frame(ref)[, "x1", drop = FALSE]),
                          residuals(ref))
  expect_equal(unname(fit$t), unname(want), tolerance = 1e-8)
})

test_that("robust SEs warn and fall back for mixed models", {
  skip_if_not_installed("lme4")
  set.seed(13)
  n <- 12
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- x1 + matrix(rnorm(n^2), n, n)

  expect_warning(
    QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
           random_intercept_sender = TRUE, use_robust_errors = TRUE,
           reps = 3, seed = 1),
    "not defined for mixed models"
  )
})


# ---- coefficient / t-value alignment -----------------------------------

test_that("aliased predictors keep coefficients and t aligned", {
  set.seed(14)
  n <- 10
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- x1 * 2            # perfectly collinear
  y  <- x1 + matrix(rnorm(n^2), n, n)

  expect_warning(
    fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
                  reps = 5, seed = 1),
    "perfect collinearity"
  )

  expect_equal(length(fit$coefficients), length(fit$t))
  expect_equal(names(fit$coefficients), names(fit$t))
  expect_true(is.na(fit$coefficients[["x2"]]))
  expect_true(is.na(fit$t[["x2"]]))
})

test_that("align_to_coefs fills gaps with NA and preserves order", {
  v <- c(b = 2, a = 1)
  out <- align_to_coefs(v, c("a", "b", "c"))
  expect_equal(names(out), c("a", "b", "c"))
  expect_equal(unname(out), c(1, 2, NA))
})


# ---- input validation --------------------------------------------------

test_that("non-conformable predictors are rejected", {
  n <- 6
  y  <- matrix(rnorm(n^2), n, n)
  x1 <- matrix(rnorm(25), 5, 5)      # wrong size
  expect_error(
    QAP(y ~ x1, data = list(y = y, x1 = x1), reps = 2),
    "dimensions"
  )

  expect_error(
    QAP(y ~ x1, data = list(y = y, x1 = rnorm(n^2)), reps = 2),
    "must be a matrix"
  )
})

test_that("a predictor list of the wrong length is rejected", {
  n <- 5
  ys <- list(matrix(rnorm(n^2), n, n), matrix(rnorm(n^2), n, n))
  xs <- list(matrix(rnorm(n^2), n, n))
  expect_error(
    QAP(y ~ x1, data = list(y = ys, x1 = xs), reps = 2),
    "network"
  )
})


# ---- per-network permutation groups ------------------------------------

test_that("perm_networks applies the grouping belonging to each network", {
  set.seed(15)
  g1 <- c(1, 1, 2, 2)
  g2 <- c(1, 2, 2, 2, 2)
  m1 <- outer(g1, g1, function(a, b) a * 10 + b)
  m2 <- outer(g2, g2, function(a, b) a * 10 + b)

  for (i in 1:25) {
    out <- perm_networks(list(m1, m2), groups = list(g1, g2))
    expect_equal(out[[1]] %/% 10, outer(g1, rep(1, length(g1))))
    expect_equal(out[[2]] %/% 10, outer(g2, rep(1, length(g2))))
  }
})

test_that("perm_networks rejects a group list of the wrong length", {
  m <- matrix(1, 3, 3)
  expect_error(perm_networks(list(m, m), groups = list(c(1, 1, 2))),
               "networks")
})

test_that("QAPglm runs with per-network groups", {
  set.seed(16)
  n <- 6
  mk <- function() matrix(rnorm(n^2), n, n)
  d <- list(y = list(mk(), mk()), x1 = list(mk(), mk()))
  g <- list(rep(1:2, each = 3), rep(1:3, each = 2))

  fit <- QAP(y ~ x1, data = d, groups = g, reps = 5, seed = 1)
  expect_s3_class(fit, "QAPRegression")
  expect_false(is.na(fit$abs[2, "x1"]))
})


# ---- p-value calibration -----------------------------------------------

test_that("permutation p-values are calibrated under the null", {
  skip_on_cran()
  set.seed(20)
  n <- 14
  reps <- 60
  n_sim <- 100

  pvals <- vapply(seq_len(n_sim), function(i) {
    x1 <- matrix(rnorm(n^2), n, n)
    x2 <- matrix(rnorm(n^2), n, n)
    y  <- matrix(rnorm(n^2), n, n)      # independent of x1 and x2
    fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
                  nullhyp = "qapspp", reps = reps)
    fit$abs[2, "x1"]
  }, numeric(1))

  # Uniform(0, 1) implies mean 0.5; allow generous Monte Carlo slack.
  expect_gt(mean(pvals), 0.35)
  expect_lt(mean(pvals), 0.65)
  # Type-I error near nominal.
  expect_lt(mean(pvals <= 0.05), 0.20)
})

test_that("a planted effect is detected", {
  set.seed(21)
  n <- 16
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- 3 * x1 + matrix(rnorm(n^2, sd = 0.4), n, n)

  fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
                nullhyp = "qapspp", reps = 200, seed = 1)

  expect_lt(fit$abs[2, "x1"], 0.05)
  expect_gt(fit$abs[2, "x2"], 0.05)
})


# ---- QAPcss flag ordering / reporting ----------------------------------

test_that("undirected CSS drops sender/receiver random intercepts", {
  skip_if_not_installed("lme4")
  set.seed(22)
  n <- 6
  y  <- array(rnorm(n^3), dim = c(n, n, n))
  x1 <- array(rnorm(n^3), dim = c(n, n, n))

  expect_warning(
    fit <- QAP(y ~ x1, data = list(y = y, x1 = x1),
                  mode = "undirected",
                  random_intercept_sender = TRUE,
                  random_intercept_perceiver = TRUE,
                  reps = 3, seed = 1),
    "sender/receiver random intercepts set to FALSE"
  )
  expect_false(fit$random[["sender"]])
  expect_false(fit$random[["receiver"]])
  expect_true(fit$random[["perceiver"]])
})

test_that("QAPcss only reports grouped permutations when groups were given", {
  set.seed(23)
  n <- 5
  y  <- array(rnorm(n^3), dim = c(n, n, n))
  x1 <- array(rnorm(n^3), dim = c(n, n, n))

  ungrouped <- QAP(y ~ x1, data = list(y = y, x1 = x1), reps = 3, seed = 1)
  expect_null(ungrouped$groups)
  expect_false(any(grepl("within groups",
                         capture.output(print(ungrouped)))))

  grouped <- QAP(y ~ x1, data = list(y = y, x1 = x1),
                    groups = c(1, 1, 2, 2, 2), reps = 3, seed = 1)
  expect_false(is.null(grouped$groups))
  expect_true(any(grepl("within groups", capture.output(print(grouped)))))
})


# ---- zero-inflated robust standard errors ------------------------------

test_that("zip robust SEs use the sandwich estimator, not HC3 on deviance residuals", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("sandwich")
  set.seed(31)
  n <- 20
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  lam <- exp(0.4 + 0.6 * x1)
  y <- matrix(ifelse(rbinom(n^2, 1, 0.3) == 1, 0, rpois(n^2, lam)), n, n)

  fit <- suppressMessages(
    QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
           family = "zip", use_robust_errors = TRUE, reps = 3, seed = 1)
  )

  ref <- fit$simple_fit
  se  <- sqrt(diag(sandwich::sandwich(ref)))
  names(se) <- sub("^count_", "", names(se))
  want <- ref$coefficients$count / se[names(ref$coefficients$count)]

  expect_equal(unname(fit$t), unname(want), tolerance = 1e-8)
  expect_false(anyNA(fit$t))
})

test_that("GMM reports robust SEs and warns that the flag is redundant", {
  skip_if_not_installed("gmm")
  set.seed(32)
  n <- 16
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- matrix(rbinom(n^2, 1, 1 / (1 + exp(-(x1)))), n, n)

  expect_warning(
    fit <- QAP(y ~ x1 + x2, data = list(y = y, x1 = x1, x2 = x2),
                  family = "binomial", estimator = "gmm",
                  use_robust_errors = TRUE, nullhyp = "qapy",
                  reps = 5, seed = 1),
    "redundant"
  )
  # gmm(vcov = "MDS") standard errors are robust, so this is still TRUE.
  expect_true(fit$robust_se)
  expect_false(anyNA(fit$t))
})
