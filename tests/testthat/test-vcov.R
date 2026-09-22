# Tests for the vcov= argument, which hands a variance-covariance
# specification straight to fixest::feglm(), and for the fixest dispatch
# around it.

skip_no_fixest <- function() skip_if_not_installed("fixest")

vcov_data <- function(n = 10, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(y  = 0.8 * x1 - 0.4 * x2 + matrix(rnorm(n^2), n, n),
       x1 = x1, x2 = x2)
}

test_that("feglm is reached with no cluster and no fixed effects", {
  skip_no_fixest()
  d <- vcov_data()
  f <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = "iid")
  expect_s3_class(f$simple_fit, "fixest")
  expect_identical(f$vcov_spec, "iid")
})

test_that("feglm is reached with fixed effects and no cluster", {
  skip_no_fixest()
  d <- vcov_data()
  # The | FE syntax alone must route through fixest; cluster = NULL is what
  # feglm() defaults to, so leaving it unset is not a fallback to glm().
  f <- QAP(y ~ x1 + x2 | sv, data = d, reps = 5, seed = 1)
  expect_s3_class(f$simple_fit, "fixest")
  expect_null(f$vcov_spec)
  expect_true(all(is.finite(f$t[c("x1", "x2")])))
})

test_that("a fixest fit without fixed effects keeps a single, real intercept", {
  skip_no_fixest()
  d <- vcov_data()
  f <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = "iid")
  expect_identical(names(f$coefficients), c("(Intercept)", "x1", "x2"))
  expect_false(is.na(f$coefficients[["(Intercept)"]]))
  expect_length(f$t, 3L)
})

test_that("an absorbed intercept is still reported as NA", {
  skip_no_fixest()
  d <- vcov_data()
  f <- QAP(y ~ x1 + x2 | sv, data = d, reps = 5, seed = 1, vcov = "hetero")
  expect_identical(names(f$coefficients), c("(Intercept)", "x1", "x2"))
  expect_true(is.na(f$coefficients[["(Intercept)"]]))
})

test_that("different vcov specifications give different standard errors", {
  skip_no_fixest()
  d <- vcov_data()
  iid <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = "iid")
  het <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = "hetero")
  # Same fit, so the coefficients must agree exactly; only the t-values move.
  expect_equal(iid$coefficients, het$coefficients)
  expect_false(isTRUE(all.equal(iid$t, het$t)))
})

test_that("vcov = ~var is the same request as fixest_se_cluster = 'var'", {
  skip_no_fixest()
  d <- vcov_data()
  by_formula <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = ~sv)
  by_cluster <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1,
                    fixest_se_cluster = "sv")
  expect_equal(unname(by_formula$t), unname(by_cluster$t))
})

test_that("a vcov formula's variables are carried into the model frame", {
  skip_no_fixest()
  d <- vcov_data()
  p <- parse_qap_formula(y ~ x1 + x2, fixest_vcov = ~sv)
  expect_true(p$use_fixest)
  expect_true("sv" %in% p$all_data_vars)
})

test_that("a character vcov does not add data variables", {
  skip_no_fixest()
  p <- parse_qap_formula(y ~ x1 + x2, fixest_vcov = "hetero")
  expect_true(p$use_fixest)
  expect_identical(p$all_data_vars, c("x1", "x2"))
})

test_that("vcov and fixest_se_cluster together are rejected", {
  skip_no_fixest()
  d <- vcov_data()
  expect_error(
    QAP(y ~ x1 + x2, data = d, reps = 5, vcov = "hetero",
        fixest_se_cluster = "sv"),
    "not both", fixed = TRUE)
})

test_that("an unusable vcov is rejected before anything is fitted", {
  skip_no_fixest()
  d <- vcov_data()
  expect_error(QAP(y ~ x1 + x2, data = d, reps = 5, vcov = 42),
               "must be a character name", fixed = TRUE)
  expect_error(QAP(y ~ x1 + x2, data = d, reps = 5,
                   vcov = c("iid", "hetero")),
               "single character string", fixed = TRUE)
})

test_that("vcov takes precedence over use_robust_errors, with a warning", {
  skip_no_fixest()
  d <- vcov_data()
  expect_warning(
    f <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = "hetero",
             use_robust_errors = TRUE),
    "takes precedence", fixed = TRUE)
  expect_false(f$robust_se)
  expect_identical(f$vcov_spec, "hetero")
})

test_that("vcov is dropped, and said to be dropped, alongside random effects", {
  skip_no_fixest()
  skip_if_not_installed("lme4")
  d <- vcov_data()
  # nullhyp = "qapy" because residualising a formula that carries an inline
  # (1|sv) is separately broken; this test is about the vcov, not that.
  expect_warning(
    f <- QAP(y ~ x1 + x2 + (1 | sv), data = d, reps = 5, seed = 1,
             nullhyp = "qapy", vcov = "hetero"),
    "ignored as well", fixed = TRUE)
  expect_null(f$vcov_spec)
})

test_that("summary reports the vcov that was actually used", {
  skip_no_fixest()
  d <- vcov_data()
  f <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = ~sv)
  out <- paste(capture.output(print(summary(f))), collapse = "\n")
  expect_match(out, "fixest vcov = ~sv", fixed = TRUE)
})

test_that("runs with different vcov specifications are not pooled", {
  skip_no_fixest()
  d <- vcov_data()
  iid <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = "iid")
  het <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 2, vcov = "hetero")
  expect_error(combine_qap_estimates(list(iid, het)),
               "same vcov_spec", fixed = TRUE)
})

test_that("two runs asking for the same vcov formula still pool", {
  skip_no_fixest()
  d <- vcov_data()
  # Each formula is built in its own frame, so the two are not identical()
  # even though they say the same thing; pooling must compare the spec.
  one <- local(QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = ~sv))
  two <- local(QAP(y ~ x1 + x2, data = d, reps = 5, seed = 2, vcov = ~sv))
  expect_false(identical(one$vcov_spec, two$vcov_spec))
  expect_s3_class(combine_qap_estimates(list(one, two)), "QAP")
})

test_that("a gaussian fixest fit reports an R-squared rather than NA", {
  skip_no_fixest()
  d <- vcov_data(n = 15)
  fx <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1, vcov = "iid")
  ols <- QAP(y ~ x1 + x2, data = d, reps = 5, seed = 1)
  expect_false(is.na(fx$r.squared))
  # With no absorbed fixed effects fixest's squared correlation is exactly
  # the R-squared the ordinary least-squares path reports.
  expect_equal(fx$r.squared, ols$r.squared, tolerance = 1e-8)
})

test_that("format_vcov_spec labels every accepted kind of specification", {
  expect_null(format_vcov_spec(NULL))
  expect_identical(format_vcov_spec("hetero"), "hetero")
  expect_identical(format_vcov_spec(~sv), "~sv")
  expect_identical(format_vcov_spec(function(x) x), "<function>")
})
