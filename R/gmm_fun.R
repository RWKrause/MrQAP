# ============================================================
# Moment conditions for the GMM estimator.
#
# Each takes gmm()'s data list, which carries y and x and -- when the model
# was fitted with them -- the prior weights w and the offset off. Weights
# scale each observation's contribution to the moment; the offset enters
# the linear predictor with a coefficient fixed at one. Both default to
# neutral, so a fit without them is unchanged.
# ============================================================

#' Prior weights from a gmm() data list, defaulting to 1
#' @param data The list passed to \code{gmm::gmm()}.
#' @return Numeric vector, or the scalar 1.
#' @keywords internal
gmm_w <- function(data) if (is.null(data$w)) 1 else as.numeric(data$w)

#' Offset from a gmm() data list, defaulting to 0
#' @param data The list passed to \code{gmm::gmm()}.
#' @return Numeric vector, or the scalar 0.
#' @keywords internal
gmm_off <- function(data) if (is.null(data$off)) 0 else as.numeric(data$off)


#' Auxiliary functions for GMM estimation - Poisson
#'
#' @param theta numeric; Vector of parameter values, see \code{?gmm::gmm}.
#' @param data matrix; Data matrix for all X, see \code{?gmm::gmm}.
#'
#' @returns Returns moment the moment conditions
#'

poisson_moments <- function(theta, data) {
  Y <- as.numeric(data$y)
  X <- data.matrix(data$x)
  lambda_hat <- exp(X %*% theta + gmm_off(data))
  residuals <- as.vector(Y - lambda_hat)
  g <- (gmm_w(data) * residuals) * X
  return(g)
}

#' Auxiliary functions for GMM estimation - Logistic
#'
#' @param theta numeric; Vector of parameter values, see \code{?gmm::gmm}.
#' @param data matrix; Data matrix for all X, see \code{?gmm::gmm}.
#'
#' @returns Returns moment the moment conditions
#'
logit_moments <- function(theta, data) {
  Y <- data$y
  X <- data$x
  prob <- 1 / (1 + exp(-1 * (X %*% theta + gmm_off(data))))
  residuals <- as.vector(Y - prob)
  g <- (gmm_w(data) * residuals) * X
  return(g)
}



#' Getting residuals for gmm logit
#'
#' @param gmmo gmm() output
#'
#' @returns residuals of the gmm()
logit_resid <- function(gmmo) {
  Y <- gmmo$dat$y
  X <- gmmo$dat$x
  prob <- 1 / (1 + exp(-1 * (X %*% gmmo$coefficients + gmm_off(gmmo$dat))))
  residuals <- as.vector(Y - prob)
  return(residuals)
}

#' Getting residuals for gmm poisson
#'
#' @param gmmo gmm() output
#'
#' @returns residuals of the gmm()
poisson_resid <- function(gmmo) {
  Y <- gmmo$dat$y
  X <- gmmo$dat$x
  lambda_hat <- exp(X %*% gmmo$coefficients + gmm_off(gmmo$dat))
  residuals <- as.vector(Y - lambda_hat)
  return(residuals)
}


#' GMM moment conditions for negative binomial regression
#'
#' The parameter vector \code{theta} contains regression coefficients
#' (length \code{ncol(X)}) followed by \code{log(alpha)} where alpha is
#' the dispersion parameter.  The mean is \eqn{\mu = \exp(X \beta)} and
#' the variance is \eqn{\mu + \alpha \mu^2}.
#'
#' @param theta numeric; parameter vector (regression coefs + log(alpha)).
#' @param data list with \code{y} and \code{x}.
#'
#' @return Matrix of moment conditions.
#' @keywords internal
negbin_moments <- function(theta, data) {
  Y <- as.numeric(data$y)
  X <- data.matrix(data$x)
  p <- ncol(X)
  beta   <- theta[1:p]
  alpha  <- exp(theta[p + 1])  # ensure alpha > 0

  mu   <- as.vector(exp(X %*% beta + gmm_off(data)))
  resid <- Y - mu
  V    <- mu + alpha * mu^2
  w    <- gmm_w(data)

  # Moment 1: E[(Y - mu) * X / V] = 0  (score-type)
  g1 <- (w * resid / V) * X
  # Moment 2: E[(Y - mu)^2 / V - 1] = 0  (variance moment)
  g2 <- w * ((resid^2 / V) - 1)

  cbind(g1, g2)
}


#' Residuals from a GMM negative binomial fit
#'
#' @param gmmo gmm() output from negbin_moments.
#'
#' @return Numeric vector of residuals.
#' @keywords internal
negbin_resid <- function(gmmo) {
  Y <- as.numeric(gmmo$dat$y)
  X <- data.matrix(gmmo$dat$x)
  p <- ncol(X)
  beta <- gmmo$coefficients[1:p]
  mu   <- as.vector(exp(X %*% beta + gmm_off(gmmo$dat)))
  as.vector(Y - mu)
}


#' GMM moment conditions for zero-inflated Poisson regression
#'
#' The parameter vector \code{theta} contains regression coefficients
#' (length \code{ncol(X)}) followed by \code{logit(pi)} where pi is the
#' zero-inflation probability.  The model is:
#' \deqn{P(Y=0) = \pi + (1-\pi) \exp(-\lambda)}
#' \deqn{P(Y=k) = (1-\pi) \lambda^k \exp(-\lambda) / k!,  k \geq 1}
#' with \eqn{\lambda = \exp(X \beta)}.
#'
#' @param theta numeric; parameter vector (regression coefs + logit(pi)).
#' @param data list with \code{y} and \code{x}.
#'
#' @return Matrix of moment conditions.
#' @keywords internal
zip_moments <- function(theta, data) {
  Y <- as.numeric(data$y)
  X <- data.matrix(data$x)
  p <- ncol(X)
  beta <- theta[1:p]
  pi_z <- 1 / (1 + exp(-theta[p + 1]))  # zero-inflation prob

  lambda <- as.vector(exp(X %*% beta + gmm_off(data)))
  # E[Y] under ZIP = (1 - pi) * lambda
  mu    <- (1 - pi_z) * lambda
  resid <- Y - mu
  w     <- gmm_w(data)

  # Moment 1: score w.r.t. beta: E[(Y - mu) * X] = 0
  g1 <- (w * resid) * X

  # Moment 2: zero indicator moment
  # P(Y=0) = pi + (1-pi)*exp(-lambda)
  p0 <- pi_z + (1 - pi_z) * exp(-lambda)
  is_zero <- as.numeric(Y == 0)
  g2 <- w * (is_zero - p0)

  cbind(g1, g2)
}


#' Residuals from a GMM zero-inflated Poisson fit
#'
#' @param gmmo gmm() output from zip_moments.
#'
#' @return Numeric vector of residuals.
#' @keywords internal
zip_resid <- function(gmmo) {
  Y <- as.numeric(gmmo$dat$y)
  X <- data.matrix(gmmo$dat$x)
  p <- ncol(X)
  beta <- gmmo$coefficients[1:p]
  pi_z <- 1 / (1 + exp(-gmmo$coefficients[p + 1]))
  lambda <- as.vector(exp(X %*% beta + gmm_off(gmmo$dat)))
  mu <- (1 - pi_z) * lambda
  as.vector(Y - mu)
}


#' Deterministic starting values for a GMM fit
#'
#' \code{gmm()} was previously started from \code{rnorm()}. That is a poor
#' start and a different one on every permutation: the Poisson, negative
#' binomial and ZIP moment conditions evaluate \code{exp(X \%*\% theta)}, which
#' overflows for a start drawn without reference to the data, so fits fail
#' routinely. A failed permutation is dropped and only warned about, and the
#' failures are not independent of the permuted data, so the surviving set is
#' not a random subsample. Surviving fits that land on different local optima
#' inject noise straight into the null distribution.
#'
#' The corresponding GLM estimate is a consistent, cheap and deterministic
#' start for all four supported families.
#'
#' @param y Numeric vector; the response.
#' @param X Numeric matrix; the design, \strong{including} the intercept
#'   column, exactly as \code{fit_qap_model()} builds it.
#' @param family Character; \code{"binomial"}, \code{"poisson"},
#'   \code{"negbin"} or \code{"zip"}.
#' @param w Numeric vector of prior weights, or NULL.
#' @param off Numeric vector of offsets, or NULL.
#'
#' @return Numeric vector of starting values: one per column of \code{X}, plus
#'   a nuisance start for \code{negbin} (\code{log(alpha)}) and \code{zip}
#'   (\code{logit(pi)}).
#' @keywords internal

gmm_start <- function(y, X, family, w = NULL, off = NULL) {
  fam <- if (family == "binomial") stats::binomial() else stats::poisson()

  b <- tryCatch(
    suppressWarnings(stats::glm.fit(
      X, y, family = fam,
      weights = if (is.null(w)) rep(1, length(y)) else w,
      offset  = if (is.null(off)) rep(0, length(y)) else off)$coefficients),
    error = function(e) rep(0, ncol(X))
  )
  # A separated or rank-deficient design yields NA/Inf coefficients; zero is
  # a usable start, an infinite one is not.
  b[!is.finite(b)] <- 0
  b <- unname(b)

  if (family == "negbin") {
    # log(alpha); alpha = 1 is the neutral start.
    return(c(b, 0))
  }
  if (family == "zip") {
    # logit(pi) from the observed zero share, clamped away from the
    # boundaries so the logit is finite.
    p0 <- mean(y == 0, na.rm = TRUE)
    p0 <- min(max(p0, 0.01), 0.99)
    return(c(b, log(p0 / (1 - p0))))
  }
  b
}
