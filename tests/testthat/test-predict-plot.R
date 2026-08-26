# ============================================================
# predict(), plot(), and the broom methods.
# ============================================================

pp_data <- function(n = 12, seed = 1) {
  set.seed(seed)
  x1 <- matrix(rnorm(n^2), n, n)
  x2 <- matrix(rnorm(n^2), n, n)
  list(n = n, y = 0.6 * x1 - 0.3 * x2 + matrix(rnorm(n^2), n, n),
       x1 = x1, x2 = x2)
}

pp_fit <- function(d, ...) {
  QAP(y ~ x1 + x2, data = list(y = d$y, x1 = d$x1, x2 = d$x2),
      reps = 60, seed = 1, ...)
}


# ---- predict ------------------------------------------------------------

test_that("predict returns one value per modelled dyad", {
  d <- pp_data(); f <- pp_fit(d)
  expect_length(predict(f), d$n * (d$n - 1))
  expect_equal(predict(f), unname(fitted(f$simple_fit)))
})

test_that("type = 'matrix' restores the network shape", {
  d <- pp_data(); f <- pp_fit(d)
  m <- predict(f, type = "matrix")

  expect_equal(dim(m), c(d$n, d$n))
  # The diagonal never entered the model, so it comes back NA.
  expect_true(all(is.na(diag(m))))
  expect_equal(sum(!is.na(m)), d$n * (d$n - 1))
  # Off-diagonal values are the vector, in the package's own order.
  expect_equal(as.vector(m)[as.vector(!diag(d$n))], predict(f))
})

test_that("type = 'link' works for a linear model", {
  # predict.lm() has no "link" type, so asking for one used to error even
  # though a linear model's link is the identity.
  d <- pp_data(); f <- pp_fit(d)
  expect_equal(predict(f, type = "link"), predict(f, type = "response"))
})

test_that("type = 'link' differs from the response scale for a GLM", {
  set.seed(2); n <- 12
  x1 <- matrix(rnorm(n^2), n, n)
  f  <- QAP(y ~ x1, data = list(y = matrix(rbinom(n^2, 1, 0.4), n, n),
                                x1 = x1),
            family = "binomial", reps = 20, seed = 1)
  lk <- predict(f, type = "link")
  rp <- predict(f, type = "response")
  expect_false(isTRUE(all.equal(lk, rp)))
  expect_equal(1 / (1 + exp(-lk)), rp, tolerance = 1e-10)
  expect_true(all(rp >= 0 & rp <= 1))
})

test_that("newdata predicts for a different network", {
  d <- pp_data(); f <- pp_fit(d)
  set.seed(9)
  nd <- list(x1 = matrix(rnorm(d$n^2), d$n, d$n),
             x2 = matrix(rnorm(d$n^2), d$n, d$n))

  p <- predict(f, newdata = nd)
  expect_length(p, d$n * (d$n - 1))
  expect_false(isTRUE(all.equal(p, predict(f))))

  # Feeding the original data back reproduces the fitted values.
  same <- predict(f, newdata = list(x1 = d$x1, x2 = d$x2))
  expect_equal(same, predict(f), tolerance = 1e-10)

  expect_equal(dim(predict(f, newdata = nd, type = "matrix")),
               c(d$n, d$n))
})

test_that("predict handles undirected and CSS shapes", {
  set.seed(3); n <- 10
  sym <- function() { m <- matrix(rnorm(n^2), n, n)
                      m[lower.tri(m)] <- t(m)[lower.tri(m)]; m }
  x1 <- sym(); y <- 0.5 * x1 + sym()
  u <- QAP(y ~ x1, data = list(y = y, x1 = x1), mode = "undirected",
           reps = 20, seed = 1)
  mu <- predict(u, type = "matrix")
  expect_equal(dim(mu), c(n, n))
  # Only the upper triangle was modelled.
  expect_true(all(is.na(mu[lower.tri(mu, diag = TRUE)])))
  expect_equal(sum(!is.na(mu)), n * (n - 1) / 2)

  a1 <- array(rnorm(6^3), c(6, 6, 6))
  cs <- QAP(y ~ x1, data = list(y = 0.5 * a1 + array(rnorm(6^3), c(6, 6, 6)),
                                x1 = a1), reps = 15, seed = 1)
  mc <- predict(cs, type = "matrix")
  expect_equal(dim(mc), c(6, 6, 6))
})

test_that("predict returns one matrix per network for multi-network fits", {
  set.seed(4); n <- 8
  mk <- function() matrix(rnorm(n^2), n, n)
  x1 <- list(mk(), mk())
  y  <- lapply(x1, function(a) 0.5 * a + matrix(rnorm(n^2), n, n))

  f <- QAP(y ~ x1, data = list(y = y, x1 = x1), reps = 20, seed = 1)
  m <- predict(f, type = "matrix")
  expect_length(m, 2L)
  expect_equal(dim(m[[1]]), c(n, n))
  expect_equal(sum(vapply(m, function(z) sum(!is.na(z)), integer(1))),
               2 * n * (n - 1))
})

test_that("predict refuses what it cannot do", {
  d <- pp_data()
  lm_less <- pp_fit(d, less_mem = TRUE)
  expect_error(predict(lm_less), "less_mem")

  f <- pp_fit(d)
  expect_error(predict(f, newdata = list(nope = d$x1)),
               "none of the model's predictors")
  expect_error(predict(f, newdata = "x"), "named list")

  set.seed(5); n <- 10
  cd <- list(y = matrix(sample(c("FP", "TN", "TP"), n^2, TRUE), n, n),
             x1 = matrix(rnorm(n^2), n, n))
  cf <- QAP(y ~ x1, data = cd, family = "binomial", reps = 10, seed = 1,
            comparison = list(commission = c("FP", "TN")))
  expect_error(predict(cf), "not defined for a model fitted with comparison")
})


# ---- plot ---------------------------------------------------------------

test_that("plot draws one panel per tested coefficient", {
  d <- pp_data(); f <- pp_fit(d)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)
  expect_invisible(plot(f))
  expect_no_error(plot(f, which = "x1"))
  expect_no_error(plot(f, statistic = "b"))
})

test_that("plot refuses the intercept under qapspp", {
  d <- pp_data(); f <- pp_fit(d, nullhyp = "qapspp")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(plot(f, which = "(Intercept)"),
               "not tested under semi-partialling")
})

test_that("plot needs the retained draws", {
  d <- pp_data(); f <- pp_fit(d, less_mem = TRUE)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(plot(f), "less_mem")
})


# ---- broom --------------------------------------------------------------

test_that("tidy follows broom's column convention", {
  skip_if_not_installed("broom")
  d <- pp_data(); f <- pp_fit(d)
  tb <- broom::tidy(f)

  expect_s3_class(tb, "data.frame")
  expect_true(all(c("term", "estimate", "statistic", "p.value") %in%
                    names(tb)))
  expect_equal(tb$estimate, unname(coef(f)))
  # std.error is deliberately absent: no single standard error corresponds
  # to the permutation p-value.
  expect_false("std.error" %in% names(tb))
})

test_that("tidy can attach permutation intervals", {
  skip_if_not_installed("broom")
  d <- pp_data(); f <- pp_fit(d)
  tb <- broom::tidy(f, conf.int = TRUE)
  expect_true(all(c("conf.low", "conf.high") %in% names(tb)))

  ci <- confint(f)
  row <- tb[tb$term == "x1", ]
  expect_equal(row$conf.low,  unname(ci["x1", 1]))
  expect_equal(row$conf.high, unname(ci["x1", 2]))
})

test_that("glance gives a one-row summary", {
  skip_if_not_installed("broom")
  d <- pp_data(); f <- pp_fit(d)
  g <- broom::glance(f)

  expect_equal(nrow(g), 1L)
  expect_equal(g$nobs, nobs(f))
  expect_equal(g$reps, f$reps)
  expect_equal(g$family, "gaussian")
  expect_equal(g$r.squared, f$r.squared)
})
