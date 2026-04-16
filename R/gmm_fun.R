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
  lambda_hat <- exp(X %*% theta)
  residuals <- as.vector(Y - lambda_hat)
  g <- residuals * X
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
  prob <- 1 / (1 + exp(-1 * (X %*% theta)))
  residuals <- as.vector(Y - prob)
  g <- residuals * X
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
  prob <- 1 / (1 + exp(-1 * (X %*% gmmo$coefficients)))
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
  lambda_hat <- exp(X %*% gmmo$coefficients)
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

  mu   <- as.vector(exp(X %*% beta))
  resid <- Y - mu
  V    <- mu + alpha * mu^2

  # Moment 1: E[(Y - mu) * X / V] = 0  (score-type)
  g1 <- (resid / V) * X
  # Moment 2: E[(Y - mu)^2 / V - 1] = 0  (variance moment)
  g2 <- (resid^2 / V) - 1

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
  mu   <- as.vector(exp(X %*% beta))
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

  lambda <- as.vector(exp(X %*% beta))
  # E[Y] under ZIP = (1 - pi) * lambda
  mu    <- (1 - pi_z) * lambda
  resid <- Y - mu

  # Moment 1: score w.r.t. beta: E[(Y - mu) * X] = 0
  g1 <- resid * X

  # Moment 2: zero indicator moment
  # P(Y=0) = pi + (1-pi)*exp(-lambda)
  p0 <- pi_z + (1 - pi_z) * exp(-lambda)
  is_zero <- as.numeric(Y == 0)
  g2 <- is_zero - p0

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
  lambda <- as.vector(exp(X %*% beta))
  mu <- (1 - pi_z) * lambda
  as.vector(Y - mu)
}
