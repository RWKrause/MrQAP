# ============================================================
# Tests for utility functions in utils.R
# ============================================================

test_that("parse_qap_formula parses simple formula", {
  p <- parse_qap_formula(y ~ x1 + x2)
  expect_equal(p$dependent, "y")
  expect_equal(p$main, c("x1", "x2"))
  expect_null(p$fixed_effects)
  expect_false(p$has_random)
  expect_false(p$use_fixest)
  expect_equal(p$all_data_vars, c("x1", "x2"))
})

test_that("parse_qap_formula detects fixest FE", {
  p <- parse_qap_formula(y ~ x1 + x2 | fe1)
  expect_equal(p$main, c("x1", "x2"))
  expect_equal(p$fixed_effects, "fe1")
  expect_true(p$use_fixest)
  expect_false(p$has_random)
  expect_equal(sort(p$all_data_vars), c("fe1", "x1", "x2"))
})

test_that("parse_qap_formula detects fixest_se_cluster", {
  p <- parse_qap_formula(y ~ x1, fixest_se_cluster = "cl")
  expect_true(p$use_fixest)
  expect_true("cl" %in% p$all_data_vars)
})

test_that("build_internal_formula adds random intercepts", {
  f <- build_internal_formula(y ~ x1, ris = TRUE, rir = TRUE)
  fs <- paste(deparse(f, width.cutoff = 500), collapse = " ")
  expect_true(grepl("1\\s*\\|\\s*sv", fs))
  expect_true(grepl("1\\s*\\|\\s*rv", fs))
})

test_that("build_internal_formula returns original when no extras", {
  f <- build_internal_formula(y ~ x1 + x2)
  expect_equal(deparse(f), deparse(y ~ x1 + x2))
})

test_that("validate_qap_input catches missing variables", {
  data <- list(y = matrix(0, 3, 3))
  parsed <- list(dependent = "y", all_data_vars = c("x1"))
  expect_error(validate_qap_input(data, parsed),
               "Predictor 'x1' not found")
})

test_that("validate_qap_input catches non-matrix", {
  data <- list(y = 1:9, x1 = matrix(0, 3, 3))
  parsed <- list(dependent = "y", all_data_vars = "x1")
  expect_error(validate_qap_input(data, parsed, css = FALSE),
               "must be a matrix")
})

test_that("compare_perm_to_baseline works with xi", {
  base <- list(coefficients = c(a = 1, b = 2),
               t            = c(a = 3, b = 4))
  res <- compare_perm_to_baseline(c(a = 0.5, b = 3),
                                   c(a = 2,   b = 5),
                                   base, xi = "b")
  expect_length(res$lower, 2)   # 2 rows (coef + t)
  expect_length(res$larger, 2)
})

test_that("compare_perm_to_baseline works without xi", {
  base <- list(coefficients = c(a = 1, b = 2),
               t            = c(a = 3, b = 4))
  res <- compare_perm_to_baseline(c(a = 0.5, b = 3),
                                   c(a = 2,   b = 5),
                                   base, xi = NULL)
  expect_true(is.matrix(res$lower))
  expect_equal(dim(res$lower), c(2, 2))
})

test_that("aggregate_perm_results computes proportions", {
  r1 <- list(lower = c(TRUE, FALSE), larger = c(FALSE, TRUE),
             abs = c(TRUE, TRUE))
  r2 <- list(lower = c(FALSE, FALSE), larger = c(TRUE, TRUE),
             abs = c(FALSE, TRUE))
  agg <- aggregate_perm_results(list(r1, r2), reps = 2)
  expect_equal(agg$lower, c(0.5, 0.0))
  expect_equal(agg$larger, c(0.5, 1.0))
  expect_equal(agg$abs, c(0.5, 1.0))
})
