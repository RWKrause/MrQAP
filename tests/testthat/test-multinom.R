# ============================================================
# Multinomial fits.
#
# A multinomial fit's coefficients are a (ncat - 1) x k MATRIX, so the
# permutation result matrices stack a block of b-comparisons above a block
# of t-comparisons rather than a single row of each. Every consumer has to
# unpack the block, not index row 2.
# ============================================================

make_mn <- function(n = 12, seed = 42) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  y  <- matrix(sample(c("A", "B", "C"), n^2, replace = TRUE), n, n)
  list(y = y, x1 = x1, x2 = x2)
}

mn_fit <- function(...) suppressWarnings(QAP(...))


test_that("multinomial runs under qapspp, the default nullhyp", {
  skip_if_not_installed("nnet")
  # Regression test: qap_init_pmats() sized the result matrices with
  # length(coefficients) -- (ncat - 1) * k -- instead of ncol(), and took
  # names() of a matrix, which is NULL. The per-predictor assignment
  # out$lower[, xi] then failed, so QAP(family = "multinom") errored
  # outright on its own default settings.
  d <- make_mn()
  expect_no_error(
    fit <- mn_fit(y ~ x1 + x2, data = d, family = "multinom",
                  nullhyp = "qapspp", reps = 10, seed = 1)
  )
  expect_equal(fit$nullhyp, "qapspp")
  expect_equal(dim(fit$abs), c(4L, 3L))   # 2 categories x (b, t) by 3 terms
})

test_that("the permutation matrices are labelled by block and category", {
  skip_if_not_installed("nnet")
  d <- make_mn()
  for (nh in c("qapy", "qapspp")) {
    fit <- mn_fit(y ~ x1 + x2, data = d, family = "multinom",
                  nullhyp = nh, reps = 10, seed = 1)
    cats <- rownames(fit$coefficients)
    expect_equal(rownames(fit$abs),
                 c(paste0("b:", cats), paste0("t:", cats)), info = nh)
    expect_equal(colnames(fit$abs), colnames(fit$coefficients), info = nh)
  }
})

test_that("summary() reports one p-value per category-by-term", {
  skip_if_not_installed("nnet")
  d   <- make_mn()
  fit <- mn_fit(y ~ x1 + x2, data = d, family = "multinom",
                nullhyp = "qapy", reps = 20, seed = 1)
  tab <- summary(fit)$coefficients

  # (ncat - 1) categories x k terms, flattened.
  expect_equal(nrow(tab), prod(dim(fit$coefficients)))
  expect_setequal(tab$term,
                  as.vector(outer(rownames(fit$coefficients),
                                  colnames(fit$coefficients),
                                  paste, sep = ":")))
  # Every one of them is populated: these used to be uniformly NA.
  expect_true(all(!is.na(tab$p_value)))
  expect_true(all(tab$p_value >= 0 & tab$p_value <= 1))
})

test_that("each reported p-value comes from its own category's t-block", {
  skip_if_not_installed("nnet")
  d   <- make_mn()
  fit <- mn_fit(y ~ x1 + x2, data = d, family = "multinom",
                nullhyp = "qapy", reps = 20, seed = 1)
  tab <- summary(fit)$coefficients

  # The claim the old code got wrong: "B:x1" must read fit$abs["t:B", "x1"],
  # not row 2 of the matrix (which is category C's *coefficient* row).
  for (cat in rownames(fit$coefficients)) {
    for (term in colnames(fit$coefficients)) {
      row <- tab[tab$term == paste(cat, term, sep = ":"), ]
      expect_equal(row$p_value,  unname(fit$abs[paste0("t:", cat), term]))
      expect_equal(row$p_lower,  unname(fit$lower[paste0("t:", cat), term]))
      expect_equal(row$p_upper,  unname(fit$larger[paste0("t:", cat), term]))
      # And the estimate it sits beside is the matching matrix cell.
      expect_equal(row$estimate, unname(fit$coefficients[cat, term]))
    }
  }
})

test_that("qapspp blanks the intercept for every category", {
  skip_if_not_installed("nnet")
  d   <- make_mn()
  fit <- mn_fit(y ~ x1 + x2, data = d, family = "multinom",
                nullhyp = "qapspp", reps = 10, seed = 1)
  tab <- summary(fit)$coefficients

  ic <- grepl(":\\(Intercept\\)$", tab$term)
  expect_equal(sum(ic), nrow(fit$coefficients))   # one per category
  expect_true(all(is.na(tab$p_value[ic])))
  # The tested predictors are not blanked.
  expect_true(all(!is.na(tab$p_value[!ic])))
})

test_that("print() shows one row per category-by-term", {
  skip_if_not_installed("nnet")
  d   <- make_mn()
  fit <- mn_fit(y ~ x1 + x2, data = d, family = "multinom",
                nullhyp = "qapy", reps = 10, seed = 1)
  out <- capture.output(print(fit))
  for (nm in as.vector(outer(rownames(fit$coefficients),
                             colnames(fit$coefficients), paste, sep = ":")))
    expect_true(any(grepl(nm, out, fixed = TRUE)), info = nm)
})

test_that("as.data.frame() matches summary() for a multinomial fit", {
  skip_if_not_installed("nnet")
  d   <- make_mn()
  fit <- mn_fit(y ~ x1 + x2, data = d, family = "multinom",
                nullhyp = "qapy", reps = 10, seed = 1)
  expect_equal(as.data.frame(fit), summary(fit)$coefficients)
})

test_that("confint() explains itself instead of blaming less_mem", {
  skip_if_not_installed("nnet")
  d   <- make_mn()
  fit <- mn_fit(y ~ x1 + x2, data = d, family = "multinom",
                nullhyp = "qapy", reps = 10, seed = 1)
  # The draws are a matrix per replication and are never stacked, so
  # null_dist is absent by construction -- refitting cannot help.
  expect_error(confint(fit), "not available for .*multinomial")
})

test_that("the reference category is respected", {
  skip_if_not_installed("nnet")
  d   <- make_mn()
  fit <- mn_fit(y ~ x1, data = d, family = "multinom", reference = "B",
                nullhyp = "qapy", reps = 10, seed = 1)
  # B is the baseline, so only A and C get coefficient rows.
  expect_setequal(rownames(fit$coefficients), c("A", "C"))
  expect_true(all(!is.na(summary(fit)$coefficients$p_value)))
})
