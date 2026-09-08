# ============================================================
# Standard extractor methods, and the permutation draws they rest on.
# ============================================================

md <- function(n = 12, seed = 1, signal = 0.4) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(y = signal * x1 + matrix(rnorm(n^2), n, n), x1 = x1, x2 = x2)
}

cat_data <- function(n = 12, seed = 5) {
  set.seed(seed)
  list(y  = matrix(sample(c("FP", "TN", "TP"), n^2, TRUE), n, n),
       x1 = matrix(rnorm(n^2), n, n),
       x2 = matrix(rnorm(n^2), n, n))
}

# ---- the retained null distribution ------------------------------------

test_that("the stored draws reproduce the reported p-values exactly", {
  for (nh in c("qapy", "qapspp")) {
    f  <- QAP(y ~ x1 + x2, data = md(), nullhyp = nh, reps = 60, seed = 1)
    nd <- f$null_dist
    expect_false(is.null(nd))
    expect_equal(colnames(nd$b), names(f$coefficients))

    for (v in c("x1", "x2")) {
      expect_equal(mean(abs(nd$t[, v]) >= abs(f$t[[v]]), na.rm = TRUE),
                   unname(f$abs["t", v]), tolerance = 1e-12,
                   info = paste(nh, v))
      expect_equal(mean(nd$t[, v] <= f$t[[v]], na.rm = TRUE),
                   unname(f$lower["t", v]), tolerance = 1e-12,
                   info = paste(nh, v))
    }
  }
})

test_that("qapspp leaves the untested intercept column empty", {
  f <- QAP(y ~ x1 + x2, data = md(), nullhyp = "qapspp", reps = 40, seed = 1)
  expect_true(all(is.na(f$null_dist$b[, "(Intercept)"])))
  expect_false(anyNA(f$null_dist$b[, "x1"]))
})

test_that("less_mem suppresses the draws", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1, less_mem = TRUE)
  expect_null(f$null_dist)
  expect_error(confint(f), "less_mem")
})

test_that("draws survive the fast and general paths identically", {
  a <- QAP(y ~ x1 + x2, data = md(), nullhyp = "qapspp", reps = 40, seed = 3)
  old <- options(MrQAP.disable_fast_ols = TRUE); on.exit(options(old))
  b <- QAP(y ~ x1 + x2, data = md(), nullhyp = "qapspp", reps = 40, seed = 3)
  expect_equal(a$null_dist, b$null_dist, tolerance = 1e-12)
})

# ---- coef --------------------------------------------------------------

test_that("coef returns the coefficient vector", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1)
  expect_equal(coef(f), f$coefficients)
  expect_named(coef(f), c("(Intercept)", "x1", "x2"))
})

test_that("coef returns one vector per comparison", {
  f <- QAP(y ~ x1 + x2, data = cat_data(), family = "binomial",
           nullhyp = "qapy", reps = 20, seed = 1,
           comparison = list(commission = c("FP", "TN"),
                             omission   = c("TP", "TN")))
  cf <- coef(f)
  expect_type(cf, "list")
  expect_named(cf, c("commission", "omission"))
  expect_named(cf$commission, c("(Intercept)", "x1", "x2"))
})

# ---- summary -----------------------------------------------------------

test_that("summary returns a usable coefficient table", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 100, seed = 1)
  s <- summary(f)

  expect_s3_class(s, "summary.QAP")
  expect_s3_class(s$coefficients, "data.frame")
  expect_equal(s$coefficients$term, c("(Intercept)", "x1", "x2"))
  expect_equal(s$coefficients$estimate, unname(coef(f)))
  expect_equal(s$reps, 100)
  expect_equal(s$nullhyp, "qapspp")
})

test_that("summary agrees with the fit it came from", {
  f <- QAP(y ~ x1 + x2, data = md(), nullhyp = "qapy", reps = 60, seed = 1)
  s <- summary(f)
  expect_equal(s$coefficients$p_value, unname(f$abs["t", ]))
  expect_equal(s$coefficients$statistic, unname(f$t))
})

test_that("summary prints the key facts", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 100, seed = 1)
  out <- capture.output(print(summary(f)))
  expect_true(any(grepl("MRQAP", out)))
  expect_true(any(grepl("semi-partialling", out)))
  expect_true(any(grepl("Coefficients", out)))
  expect_true(any(grepl("R-squared", out)))
  expect_true(any(grepl("100 permutations", out)))
})

test_that("summary works for CSS and GLM fits", {
  set.seed(2)
  m <- 6
  a1 <- array(rnorm(m^3), dim = c(m, m, m))
  cs <- QAP(y ~ x1, data = list(y = 0.5 * a1 + array(rnorm(m^3), dim = dim(a1)),
                                x1 = a1), reps = 20, seed = 1)
  expect_s3_class(summary(cs), "summary.QAP")
  expect_true(any(grepl("Cognitive Social Structure",
                        capture.output(print(summary(cs))))))

  set.seed(3)
  n <- 12
  x1 <- matrix(rnorm(n^2), n, n)
  yb <- matrix(rbinom(n^2, 1, 1 / (1 + exp(-x1))), n, n)
  g <- QAP(y ~ x1, data = list(y = yb, x1 = x1), family = "binomial",
           reps = 20, seed = 1)
  sg <- summary(g)
  # non-gaussian families get the exponentiated estimate
  expect_true("exp_estimate" %in% names(sg$coefficients))
  expect_equal(sg$coefficients$exp_estimate, exp(sg$coefficients$estimate))
})

test_that("summary handles comparison models", {
  f <- QAP(y ~ x1 + x2, data = cat_data(), family = "binomial",
           nullhyp = "qapy", reps = 20, seed = 1,
           comparison = list(commission = c("FP", "TN")))
  s <- summary(f)
  expect_type(s$coefficients, "list")
  expect_named(s$coefficients, "commission")
  out <- capture.output(print(s))
  expect_true(any(grepl("Comparison: commission", out)))
})

# ---- p-value formatting ------------------------------------------------

test_that("a permutation p-value of 0 is shown as below the resolution", {
  # printCoefmat() would render this as <2e-16, inventing precision that
  # 200 permutations cannot possibly deliver.
  expect_equal(format_perm_p(0, 200), "<0.005")
  expect_equal(format_perm_p(0, 1000), "<0.001")
  expect_equal(format_perm_p(NA_real_, 100), "-")
  expect_equal(format_perm_p(0.25, 100), "0.25")
})

test_that("printed output never claims 2e-16", {
  f <- QAP(y ~ x1 + x2, data = md(signal = 3), nullhyp = "qapy",
           reps = 100, seed = 1)
  out <- capture.output(print(summary(f)))
  expect_false(any(grepl("2e-16", out)))
  expect_true(any(grepl("<0.01", out)))
  expect_true(any(grepl("smallest distinguishable", out)))
})

# ---- confint -----------------------------------------------------------

test_that("confint uses the permutation spread, not the null percentiles", {
  f <- QAP(y ~ x1 + x2, data = md(n = 20, signal = 0.4), nullhyp = "qapy",
           reps = 300, seed = 1)
  ci <- confint(f)

  expect_equal(dim(ci), c(3L, 2L))
  expect_equal(rownames(ci), names(coef(f)))

  # centred on the estimate
  expect_equal(rowMeans(ci), coef(f), tolerance = 1e-8)

  # width is 2 * z * sd(permutation draws)
  s <- apply(f$null_dist$b, 2, sd, na.rm = TRUE)
  expect_equal(unname(ci[, 2] - ci[, 1]),
               unname(2 * qnorm(0.975) * s), tolerance = 1e-8)

  # a percentile interval of the null would NOT contain the estimate
  pct <- quantile(f$null_dist$b[, "x1"], c(.025, .975))
  expect_false(coef(f)[["x1"]] > pct[1] && coef(f)[["x1"]] < pct[2])
})

test_that("confint respects level and parm", {
  f <- QAP(y ~ x1 + x2, data = md(), nullhyp = "qapy", reps = 200, seed = 1)
  w95 <- diff(confint(f)["x1", ])
  w99 <- diff(confint(f, level = 0.99)["x1", ])
  expect_gt(w99, w95)

  one <- confint(f, parm = "x1")
  expect_equal(dim(one), c(1L, 2L))
  expect_equal(rownames(one), "x1")
})

test_that("confint gives NA for the untested intercept under qapspp", {
  f <- QAP(y ~ x1 + x2, data = md(), nullhyp = "qapspp", reps = 60, seed = 1)
  ci <- confint(f)
  expect_true(all(is.na(ci["(Intercept)", ])))
  expect_false(anyNA(ci["x1", ]))
})

test_that("confint returns one matrix per comparison", {
  f <- QAP(y ~ x1 + x2, data = cat_data(), family = "binomial",
           nullhyp = "qapy", reps = 40, seed = 1,
           comparison = list(commission = c("FP", "TN")))
  ci <- confint(f)
  expect_type(ci, "list")
  expect_named(ci, "commission")
  expect_equal(dim(ci$commission), c(3L, 2L))
})

# ---- the rest of the API -----------------------------------------------

test_that("nobs counts the dyads actually used", {
  n <- 12
  f <- QAP(y ~ x1 + x2, data = md(n = n), reps = 20, seed = 1)
  expect_equal(nobs(f), n * (n - 1))          # diagonal dropped

  fd <- QAP(y ~ x1 + x2, data = md(n = n), diag = TRUE, reps = 20, seed = 1)
  expect_equal(nobs(fd), n^2)
})

test_that("fitted, residuals, logLik, AIC and BIC delegate to the baseline", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1)
  expect_length(fitted(f), nobs(f))
  expect_length(residuals(f), nobs(f))
  expect_s3_class(logLik(f), "logLik")
  expect_equal(AIC(f), AIC(f$simple_fit))
  expect_equal(BIC(f), BIC(f$simple_fit))
  expect_equal(unname(fitted(f)), unname(fitted(f$simple_fit)))
})

test_that("extractors needing the model object fail clearly under less_mem", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1, less_mem = TRUE)
  expect_error(fitted(f), "less_mem")
  expect_error(residuals(f), "less_mem")
  expect_error(logLik(f), "less_mem")
  expect_error(vcov(f), "less_mem")
  # nobs() does NOT need the model object: the count is recorded when the
  # data are vectorised, so less_mem no longer costs it.
  full <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1)
  expect_equal(nobs(f), nobs(full))
  expect_false(is.na(nobs(f)))
  # coef still works: it does not need the model object
  expect_named(coef(f), c("(Intercept)", "x1", "x2"))
  # summary() too: it reads coefficients and the p-value matrices.
  expect_equal(nrow(summary(f)$coefficients), 3L)
})

test_that("less_mem actually drops the model object", {
  # It used to skip only the top-level copy while the same object survived
  # in fit$base$base_model, so the fit was no smaller (4.68 of 4.69 MB at
  # n = 120) even as the extractors refused to run.
  f <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1, less_mem = TRUE)
  expect_null(f$simple_fit)
  expect_null(f$base$base_model)
  expect_null(f$null_dist)

  full <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1)
  expect_lt(as.numeric(object.size(f)),
            as.numeric(object.size(full)) / 2)
  # And it is the same model: only the storage differs.
  expect_equal(coef(f), coef(full))
  expect_equal(f$abs, full$abs)
})

test_that("less_mem drops the model for comparison fits too", {
  f <- QAP(y ~ x1 + x2, data = cat_data(), family = "binomial",
           comparison = list(commission = c("FP", "TN")),
           reps = 10, seed = 1, less_mem = TRUE)
  expect_null(f$simple_fits)
  expect_true(all(vapply(f$base, function(b) is.null(b$base_model),
                         logical(1))))
  expect_named(coef(f), "commission")
})

test_that("vcov returns the model-based covariance", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1)
  expect_equal(vcov(f), vcov(f$simple_fit))
})

test_that("as.data.frame gives a tidy table", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 40, seed = 1)
  d <- as.data.frame(f)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 3)
  expect_true(all(c("term", "estimate", "statistic", "p_value") %in% names(d)))
  expect_equal(d$estimate, unname(coef(f)))
})

test_that("as.data.frame stacks comparisons with a label column", {
  f <- QAP(y ~ x1 + x2, data = cat_data(), family = "binomial",
           nullhyp = "qapy", reps = 20, seed = 1,
           comparison = list(commission = c("FP", "TN"),
                             omission   = c("TP", "TN")))
  d <- as.data.frame(f)
  expect_true("comparison" %in% names(d))
  expect_setequal(unique(d$comparison), c("commission", "omission"))
  expect_equal(nrow(d), 6)
})

test_that("the shared parent class does not disturb existing dispatch", {
  f <- QAP(y ~ x1 + x2, data = md(), reps = 20, seed = 1)
  expect_s3_class(f, "QAPRegression")
  expect_s3_class(f, "QAP")
  expect_output(print(f), "OLS Network Model")
})
