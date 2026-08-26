# ============================================================
# Golden-master case definitions.
#
# Shared by the fixture generator (data-raw/make_golden.R) and
# test-golden.R, so the two can never drift apart.
#
# Every case pins `seed` and states `nullhyp` EXPLICITLY, so the
# fixture does not depend on any default that later refactors change.
# Keep the cases small: the whole grid runs on every test pass.
# ============================================================

# ---- deterministic data builders ---------------------------------------

golden_matrix_data <- function(n = 10, seed = 101, na = FALSE) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- 1 + 1.5 * x1 - 0.7 * x2 + matrix(rnorm(n^2, sd = 0.5), n, n)
  if (na) {
    y[2, 3]  <- NA
    x1[5, 6] <- NA
  }
  list(y = y, x1 = x1, x2 = x2)
}

golden_symmetric_data <- function(n = 10, seed = 102) {
  set.seed(seed)
  sym <- function() {
    m <- matrix(rnorm(n^2), n, n)
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  x1 <- sym(); x2 <- sym()
  # Noise must itself be symmetric, or the two triangles disagree; without
  # any noise the fit is perfect and every t-value is Inf, which makes a
  # useless (and numerically brittle) golden case.
  e  <- sym()
  y  <- 1 + 1.5 * x1 - 0.7 * x2 + 0.5 * e
  list(y = y, x1 = x1, x2 = x2)
}

# Weak signal, so the permutation p-values land strictly between 0 and 1.
# The strongly-signalled cases above saturate every p-value at 0, which makes
# them blind to changes in how the PERMUTATION t-values are computed.
golden_weak_data <- function(n = 12, seed = 108, na = FALSE) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- 0.12 * x1 + matrix(rnorm(n^2), n, n)
  if (na) y[3, 4] <- NA
  list(y = y, x1 = x1, x2 = x2)
}

golden_binary_data <- function(n = 12, seed = 103) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  eta <- -0.4 + 1.2 * x1 - 0.5 * x2
  y <- matrix(rbinom(n^2, 1, 1 / (1 + exp(-eta))), n, n)
  list(y = y, x1 = x1, x2 = x2)
}

golden_multinet_data <- function(n = 8, seed = 104) {
  set.seed(seed)
  mk <- function() matrix(rnorm(n^2), n, n)
  x1 <- list(mk(), mk()); x2 <- list(mk(), mk())
  y  <- Map(function(a, b) 1 + a - 0.5 * b + matrix(rnorm(n^2, sd = .5), n, n),
            x1, x2)
  list(y = y, x1 = x1, x2 = x2)
}

golden_css_data <- function(n = 7, seed = 105) {
  set.seed(seed)
  mk <- function() array(rnorm(n^3), dim = c(n, n, n))
  x1 <- mk(); x2 <- mk()
  y  <- 1 + x1 - 0.5 * x2 + array(rnorm(n^3, sd = 0.5), dim = c(n, n, n))
  list(y = y, x1 = x1, x2 = x2)
}

golden_css_binary_data <- function(n = 7, seed = 106) {
  set.seed(seed)
  mk <- function() array(rnorm(n^3), dim = c(n, n, n))
  x1 <- mk(); x2 <- mk()
  eta <- -0.3 + x1
  y <- array(rbinom(n^3, 1, 1 / (1 + exp(-eta))), dim = c(n, n, n))
  list(y = y, x1 = x1, x2 = x2)
}

golden_symmetric_weak <- function(n = 12, seed = 109) {
  set.seed(seed)
  sym <- function() {
    m <- matrix(rnorm(n^2), n, n)
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  x1 <- sym(); x2 <- sym()
  y  <- 0.12 * x1 + sym()
  list(y = y, x1 = x1, x2 = x2)
}

golden_two_mode_data <- function(nr = 9, nc = 6, seed = 110) {
  set.seed(seed)
  x1 <- matrix(rnorm(nr * nc), nr, nc)
  x2 <- matrix(rnorm(nr * nc), nr, nc)
  list(y = 0.15 * x1 + matrix(rnorm(nr * nc), nr, nc), x1 = x1, x2 = x2)
}

golden_category_data <- function(n = 12, seed = 107) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- matrix(sample(c("FP", "TN", "TP"), n^2, replace = TRUE), n, n)
  list(y = y, x1 = x1, x2 = x2)
}

# ---- the grid ----------------------------------------------------------
#
# Each entry is a zero-argument function returning a fitted object.
# NOTE: when QAP()/QAP() are replaced by QAP(), rewrite the calls
# here only. The stored fixture values must not change.

golden_cases <- function() {
  list(
    gauss_qapy = function()
      QAP(y ~ x1 + x2, data = golden_matrix_data(),
             nullhyp = "qapy", reps = 25, seed = 1),

    gauss_qapspp = function()
      QAP(y ~ x1 + x2, data = golden_matrix_data(),
             nullhyp = "qapspp", reps = 25, seed = 2),

    gauss_undirected = function()
      QAP(y ~ x1 + x2, data = golden_symmetric_data(),
             mode = "undirected", nullhyp = "qapspp", reps = 25, seed = 3),

    gauss_robust = function()
      QAP(y ~ x1 + x2, data = golden_matrix_data(),
             use_robust_errors = TRUE,
             nullhyp = "qapy", reps = 25, seed = 4),

    gauss_diag = function()
      QAP(y ~ x1 + x2, data = golden_matrix_data(),
             diag = TRUE, nullhyp = "qapy", reps = 25, seed = 5),

    gauss_groups = function()
      QAP(y ~ x1 + x2, data = golden_matrix_data(),
             groups = rep(1:2, each = 5),
             nullhyp = "qapspp", reps = 25, seed = 6),

    binom_qapy = function()
      QAP(y ~ x1 + x2, data = golden_binary_data(),
             family = "binomial", nullhyp = "qapy", reps = 25, seed = 7),

    binom_qapspp = function()
      QAP(y ~ x1 + x2, data = golden_binary_data(),
             family = "binomial", nullhyp = "qapspp", reps = 25, seed = 8),

    binom_robust = function()
      QAP(y ~ x1 + x2, data = golden_binary_data(),
             family = "binomial", use_robust_errors = TRUE,
             nullhyp = "qapy", reps = 25, seed = 9),

    multinet_qapy = function()
      QAP(y ~ x1 + x2, data = golden_multinet_data(),
             nullhyp = "qapy", reps = 25, seed = 10),

    multinet_qapspp_groups = function()
      QAP(y ~ x1 + x2, data = golden_multinet_data(),
             groups = list(rep(1:2, each = 4), rep(1:4, each = 2)),
             nullhyp = "qapspp", reps = 25, seed = 11),

    na_qapy = function()
      QAP(y ~ x1 + x2, data = golden_matrix_data(na = TRUE),
             nullhyp = "qapy", reps = 25, seed = 12),

    na_qapspp = function()
      QAP(y ~ x1 + x2, data = golden_matrix_data(na = TRUE),
             nullhyp = "qapspp", reps = 25, seed = 13),

    comparison = function()
      QAP(y ~ x1 + x2, data = golden_category_data(),
             family = "binomial", nullhyp = "qapy", reps = 25, seed = 14,
             comparison = list(commission = c("FP", "TN"))),

    css_qapy = function()
      QAP(y ~ x1 + x2, data = golden_css_data(),
             nullhyp = "qapy", reps = 15, seed = 15),

    css_qapspp = function()
      QAP(y ~ x1 + x2, data = golden_css_data(),
             nullhyp = "qapspp", reps = 15, seed = 16),

    css_binom = function()
      QAP(y ~ x1 + x2, data = golden_css_binary_data(),
             family = "binomial", nullhyp = "qapy", reps = 15, seed = 17),

    # --- weak signal: p-values must be informative, not saturated -----
    weak_qapy = function()
      QAP(y ~ x1 + x2, data = golden_weak_data(),
          nullhyp = "qapy", reps = 40, seed = 20),

    weak_qapspp = function()
      QAP(y ~ x1 + x2, data = golden_weak_data(),
          nullhyp = "qapspp", reps = 40, seed = 21),

    weak_robust_qapy = function()
      QAP(y ~ x1 + x2, data = golden_weak_data(),
          use_robust_errors = TRUE,
          nullhyp = "qapy", reps = 40, seed = 22),

    weak_robust_qapspp = function()
      QAP(y ~ x1 + x2, data = golden_weak_data(),
          use_robust_errors = TRUE,
          nullhyp = "qapspp", reps = 40, seed = 23),

    weak_undirected = function()
      QAP(y ~ x1 + x2, data = golden_symmetric_weak(),
          mode = "undirected", nullhyp = "qapspp", reps = 40, seed = 24),

    weak_na = function()
      QAP(y ~ x1 + x2, data = golden_weak_data(na = TRUE),
          nullhyp = "qapspp", reps = 40, seed = 25),

    weak_two_mode = function()
      QAP(y ~ x1 + x2, data = golden_two_mode_data(),
          nullhyp = "qapspp", reps = 40, seed = 26),

    css_groups = function()
      QAP(y ~ x1 + x2, data = golden_css_data(),
             groups = c(1, 1, 1, 2, 2, 2, 2),
             nullhyp = "qapy", reps = 15, seed = 18)
  )
}

# ---- what we compare ---------------------------------------------------
#
# Reduce a fit to the numbers that define it. Deliberately excludes the
# stored model object, which holds calls and environments that are not
# stable across a refactor and are not results.

golden_snapshot <- function(fit) {
  num <- function(v) {
    if (is.null(v)) return(NULL)
    if (is.list(v)) return(lapply(v, num))
    storage.mode(v) <- "double"
    v
  }
  # For a comparison model fit$base is a LIST of fits, one per comparison,
  # and fit$coefficients is absent; pull each comparison's estimates so the
  # case actually constrains something.
  if (is.null(fit$coefficients) && is.list(fit$base) &&
      is.null(fit$base$coefficients)) {
    base <- list(coefficients = lapply(fit$base, `[[`, "coefficients"),
                 t            = lapply(fit$base, `[[`, "t"))
  } else {
    base <- if (!is.null(fit$coefficients)) fit else fit$base
  }
  list(
    coefficients = num(base$coefficients),
    t            = num(base$t),
    lower        = num(fit$lower),
    larger       = num(fit$larger),
    abs          = num(fit$abs),
    r.squared    = num(base$r.squared),
    theta        = num(base$theta),
    nullhyp      = fit$nullhyp,
    reps         = fit$reps,
    family       = fit$family
  )
}

golden_run_all <- function() {
  cases <- golden_cases()
  out <- vector("list", length(cases))
  names(out) <- names(cases)
  for (nm in names(cases)) {
    out[[nm]] <- golden_snapshot(suppressWarnings(suppressMessages(cases[[nm]]())))
  }
  out
}
