# Random effects written directly into the user's formula.
#
# These used to work under nullhyp = "qapy" and crash under the default
# "qapspp": the residualisation step took its random structure from the
# random_intercept_* flags only, so a formula-supplied (1|g) left
# has_random TRUE -- which routes residualisation through lme4 -- while the
# formula handed to lmer() had no random part at all.

skip_no_lme4 <- function() skip_if_not_installed("lme4")

re_data <- function(n = 10, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(y  = 0.8 * x1 - 0.4 * x2 + matrix(rnorm(n^2), n, n),
       x1 = x1, x2 = x2,
       # A grouping variable the user supplies, as opposed to the structural
       # sv/rv/nv/pv that QAP() generates.
       grp = matrix(rep(seq_len(n), each = n), n, n))
}

re_fit <- function(...) suppressWarnings(suppressMessages(QAP(...)))


# ---- parse_qap_formula --------------------------------------------------

test_that("parse_qap_formula returns the random terms the user wrote", {
  p <- parse_qap_formula(y ~ x1 + x2 + (1 | sv))
  expect_true(p$has_random)
  expect_identical(p$random_terms, "(1 | sv)")
})

test_that("parse_qap_formula carries grouping variables into the data vars", {
  p <- parse_qap_formula(y ~ x1 + x2 + (1 | grp))
  # Without this the grouping column never reaches the model frame and the
  # fit failed with "object 'grp' not found".
  expect_true("grp" %in% p$all_data_vars)
  expect_identical(p$main, c("x1", "x2"))
})

test_that("parse_qap_formula reports no random terms for other shapes", {
  expect_identical(parse_qap_formula(y ~ x1 + x2)$random_terms, character(0))
  # findbars() would read the whole of "x1 + x2 | fe" as a bar, so it must
  # not be called on a fixest formula.
  p <- parse_qap_formula(y ~ x1 + x2 | fe)
  expect_identical(p$random_terms, character(0))
  expect_true(p$use_fixest)
})


# ---- the reported crash -------------------------------------------------

test_that("an inline random effect runs under both null hypotheses", {
  skip_no_lme4()
  d <- re_data()
  spp <- re_fit(y ~ x1 + x2 + (1 | sv), data = d, reps = 3, seed = 1)
  qy  <- re_fit(y ~ x1 + x2 + (1 | sv), data = d, reps = 3, seed = 1,
                nullhyp = "qapy")

  expect_equal(spp$nullhyp, "qapspp")
  expect_equal(qy$nullhyp, "qapy")
  # The baseline fit is the same model either way, so the estimates must
  # agree exactly; only the null distribution differs.
  expect_equal(spp$coefficients, qy$coefficients)
  expect_true(all(is.finite(spp$t[c("x1", "x2")])))
})

test_that("several inline random effects run under qapspp", {
  skip_no_lme4()
  d <- re_data()
  f <- re_fit(y ~ x1 + x2 + (1 | sv) + (1 | rv), data = d, reps = 3, seed = 1)
  expect_true(all(is.finite(f$t[c("x1", "x2")])))
})

test_that("a user-named grouping variable works under both null hypotheses", {
  skip_no_lme4()
  d <- re_data()
  spp <- re_fit(y ~ x1 + x2 + (1 | grp), data = d, reps = 3, seed = 1)
  qy  <- re_fit(y ~ x1 + x2 + (1 | grp), data = d, reps = 3, seed = 1,
                nullhyp = "qapy")
  expect_true(all(is.finite(spp$t[c("x1", "x2")])))
  expect_equal(spp$coefficients, qy$coefficients)
})


# ---- both sources of random effects at once -----------------------------

test_that("an inline random effect combines with a flag-generated one", {
  skip_no_lme4()
  d <- re_data()
  both <- re_fit(y ~ x1 + x2 + (1 | sv), data = d, reps = 3, seed = 1,
                 random_intercept_receiver = TRUE)
  # Same model as writing both terms out by hand.
  written <- re_fit(y ~ x1 + x2 + (1 | sv) + (1 | rv), data = d, reps = 3,
                    seed = 1)
  expect_equal(both$coefficients, written$coefficients)
  expect_true(both$random[["receiver"]])
})

test_that("a flag does not duplicate a term the formula already has", {
  skip_no_lme4()
  d <- re_data()
  dup  <- re_fit(y ~ x1 + x2 + (1 | sv), data = d, reps = 3, seed = 1,
                 random_intercept_sender = TRUE)
  once <- re_fit(y ~ x1 + x2 + (1 | sv), data = d, reps = 3, seed = 1)
  # lme4 accepts (1|sv) twice but then fits two variance components for one
  # grouping, which is neither what the flag nor the formula asked for.
  expect_equal(dup$coefficients, once$coefficients)
  expect_equal(dup$t, once$t)
})

test_that("build_internal_formula drops an already-present term", {
  expect_identical(
    deparse(build_internal_formula(y ~ x1 + (1 | sv), ris = TRUE,
                                   existing_terms = "(1 | sv)")),
    deparse(y ~ x1 + (1 | sv)))
  expect_identical(
    deparse(build_internal_formula(y ~ x1 + (1 | sv), rir = TRUE,
                                   existing_terms = "(1 | sv)")),
    deparse(y ~ x1 + (1 | sv) + (1 | rv)))
})


# ---- has_random is about random effects, not about parentheses ----------

test_that("a transformed predictor is rejected, and clearly", {
  d <- re_data()
  d$x1 <- abs(d$x1) + 1
  # This was read as a mixed model because the formula contains "(", and
  # failed inside lme4 complaining about missing random effects.
  expect_error(QAP(y ~ log(x1) + x2, data = d, reps = 3),
               "must be plain variables", fixed = TRUE)
  expect_error(QAP(y ~ I(x1^2) + x2, data = d, reps = 3),
               "must be plain variables", fixed = TRUE)
})

test_that("an interaction of plain variables is still allowed", {
  d <- re_data()
  # x1 * x2 puts x1:x2 in the model matrix while still permuting x1 and x2
  # themselves, so it is not a transformed predictor.
  f <- re_fit(y ~ x1 * x2, data = d, reps = 3, seed = 1)
  expect_true("x1:x2" %in% names(f$coefficients))
})

test_that("a plain model is not treated as having random effects", {
  d <- re_data()
  f <- re_fit(y ~ x1 + x2, data = d, reps = 3, seed = 1)
  expect_false(any(f$random))
  expect_s3_class(f$simple_fit, "lm")
})
