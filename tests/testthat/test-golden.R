# ============================================================
# Golden master.
#
# Locks the numeric output of an 18-case grid so the API unification,
# multi-mode generalisation and hot-path optimisation cannot silently
# change a result.
#
# A FAILURE HERE MEANS A REFACTOR CHANGED RESULTS. Diagnose it. Do not
# regenerate the fixture unless the change is intended and understood
# (data-raw/make_golden.R).
# ============================================================

golden_path <- test_path("fixtures", "golden.rds")

test_that("the golden fixture is present", {
  expect_true(file.exists(golden_path))
})

test_that("all golden cases reproduce exactly", {
  skip_if_not(file.exists(golden_path))
  expected <- readRDS(golden_path)
  actual   <- golden_run_all()

  expect_named(actual, names(expected))

  for (nm in names(expected)) {
    expect_equal(actual[[nm]], expected[[nm]],
                 tolerance = 1e-10,
                 info = paste0("golden case '", nm, "' changed"))
  }
})

test_that("the grid still covers every path it is supposed to", {
  # Guards against a case being quietly dropped from helper-golden.R.
  expected <- readRDS(golden_path)
  nms <- names(expected)

  expect_true(all(c("gauss_qapy", "gauss_qapspp") %in% nms))
  expect_true("gauss_undirected" %in% nms)
  expect_true(any(grepl("^binom", nms)))
  expect_true(any(grepl("^multinet", nms)))
  expect_true(any(grepl("^na_", nms)))
  expect_true("comparison" %in% nms)
  expect_true(any(grepl("^css_", nms)))
  expect_true(any(grepl("groups", nms)))

  # Non-degenerate: coefficients and t must be finite, or the case
  # constrains nothing.
  for (nm in nms) {
    cf <- expected[[nm]]$coefficients
    tv <- expected[[nm]]$t
    flat <- function(v) if (is.list(v)) unlist(v) else v
    expect_true(all(is.finite(flat(cf))),
                info = paste(nm, "has non-finite coefficients"))
    expect_true(all(is.finite(flat(tv))),
                info = paste(nm, "has non-finite t-values"))
  }
})
