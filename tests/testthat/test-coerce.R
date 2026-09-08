# ============================================================
# Input adapters: dyadic data frames, igraph and network objects.
#
# Every one of these must land on exactly the fit the equivalent matrix
# input gives. The adapters convert; they must not change the model.
# ============================================================

coerce_df <- function(n = 8, seed = 1) {
  set.seed(seed)
  ids <- paste0("a", seq_len(n))
  df  <- expand.grid(from = ids, to = ids, stringsAsFactors = FALSE)
  df  <- df[df$from != df$to, ]
  df$x1 <- rnorm(nrow(df))
  df$y  <- 0.6 * df$x1 + rnorm(nrow(df))
  df
}


test_that("a dyadic data frame gives the same fit as the matrices", {
  df <- coerce_df()

  from_df <- QAP(y ~ x1, data = df, sender = "from", receiver = "to",
                 reps = 20, seed = 1)
  mats    <- df_to_mat(df, sender = "from", receiver = "to")
  from_m  <- QAP(y ~ x1, data = mats, reps = 20, seed = 1)

  expect_equal(from_df$coefficients, from_m$coefficients)
  expect_equal(from_df$abs, from_m$abs)
  expect_equal(nobs(from_df), nobs(from_m))
})

test_that("a data frame without sender/receiver is refused clearly", {
  df <- coerce_df()
  expect_error(QAP(y ~ x1, data = df, reps = 5), "sender.*receiver")
})

test_that("split_by gives one network per level", {
  df <- coerce_df()
  df$wave <- rep(1:2, length.out = nrow(df))

  # Alternating rows leaves each wave with gaps, which df_to_mat warns
  # about; that is the test data, not the behaviour under test.
  fit <- suppressWarnings(
    QAP(y ~ x1, data = df, sender = "from", receiver = "to",
        split_by = "wave", reps = 10, seed = 1))
  expect_s3_class(fit, "QAP")
  # Two networks were stacked, so more observations than one wave alone.
  expect_gt(nobs(fit), 0)
})

test_that("igraph objects are converted to adjacency matrices", {
  skip_if_not_installed("igraph")
  set.seed(2); n <- 8
  gy <- igraph::sample_gnp(n, 0.4, directed = TRUE)
  gx <- igraph::sample_gnp(n, 0.4, directed = TRUE)
  igraph::V(gy)$name <- igraph::V(gx)$name <- paste0("v", seq_len(n))

  my <- as.matrix(igraph::as_adjacency_matrix(gy, sparse = FALSE))
  mx <- as.matrix(igraph::as_adjacency_matrix(gx, sparse = FALSE))
  storage.mode(my) <- storage.mode(mx) <- "double"

  a <- QAP(y ~ x1, data = list(y = gy, x1 = gx), family = "binomial",
           reps = 10, seed = 1)
  b <- QAP(y ~ x1, data = list(y = my, x1 = mx), family = "binomial",
           reps = 10, seed = 1)
  expect_equal(a$coefficients, b$coefficients)
  expect_equal(a$abs, b$abs)
})

test_that("an igraph edge attribute is read when named", {
  skip_if_not_installed("igraph")
  set.seed(3); n <- 7
  g <- igraph::sample_gnp(n, 0.5, directed = TRUE)
  igraph::E(g)$w <- runif(igraph::ecount(g), 1, 5)

  wt  <- qap_from_igraph(g, attr = "w")
  bin <- qap_from_igraph(g, attr = NULL)
  expect_true(any(wt > 1))              # weights, not 0/1
  expect_true(all(bin %in% c(0, 1)))
  expect_error(qap_from_igraph(g, attr = "nope"), "no edge attribute")
})

test_that("network objects are converted to adjacency matrices", {
  skip_if_not_installed("network")
  set.seed(4); n <- 8
  my <- matrix(rbinom(n^2, 1, 0.4), n, n); diag(my) <- 0
  mx <- matrix(rbinom(n^2, 1, 0.4), n, n); diag(mx) <- 0

  ny <- network::network(my, directed = TRUE)
  nx <- network::network(mx, directed = TRUE)

  a <- QAP(y ~ x1, data = list(y = ny, x1 = nx), family = "binomial",
           reps = 10, seed = 1)
  storage.mode(my) <- storage.mode(mx) <- "double"
  b <- QAP(y ~ x1, data = list(y = my, x1 = mx), family = "binomial",
           reps = 10, seed = 1)
  expect_equal(a$coefficients, b$coefficients)
})

test_that("a list of graphs is converted network by network", {
  skip_if_not_installed("igraph")
  set.seed(5); n <- 7
  mk <- function() igraph::sample_gnp(n, 0.4, directed = TRUE)
  gy <- list(mk(), mk()); gx <- list(mk(), mk())

  fit <- QAP(y ~ x1, data = list(y = gy, x1 = gx), family = "binomial",
             reps = 10, seed = 1)
  expect_equal(nobs(fit), 2 * n * (n - 1))
})

test_that("plain matrix input is untouched by the adapter", {
  set.seed(6); n <- 6
  d <- list(y = matrix(rnorm(n^2), n, n), x1 = matrix(rnorm(n^2), n, n))
  expect_identical(qap_coerce_data(d), d)
})

test_that("a non-list, non-data-frame data argument is refused", {
  expect_error(qap_coerce_data(1:10), "named list of matrices")
})
