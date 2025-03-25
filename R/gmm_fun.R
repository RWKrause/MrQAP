#' Auxiliary functions for estimation gmm
#'
#' @param theta numeric; Vetor of parameter values, see \code{?gmm::gmm}.
#' @param data matrix; Data matrix for all X, see \code{?gmm::gmm}.
#'
#' @returns Returns moment the moment conditions
#'

logit_moments <- function(theta, data) {
  Y <- data$y
  X <- data$x
  prob <- 1 / (1 + exp(-X %*% theta))
  residuals <- as.vector(Y - prob)
  g <- residuals * X
  return(g)
}

poisson_moments <- function(theta, data) {
  Y <- data$y
  X <- data$x
  lambda_hat <- exp(X %*% theta)
  residuals <- as.vector(y - lambda_hat)
  g <- residuals * X
  return(g)
}

logit_resid <- function(gmmo) {
  Y <- gmmo$dat$y
  X <- gmmo$dat$X
  prob <- 1 / (1 + exp(-X %*% gmmo$coefficients))
  residuals <- as.vector(Y - prob)
  return(residuals)
}

poisson_resid <- function(gmmo) {
  Y <- gmmo$dat$y
  X <- gmmo$dat$X
  lambda_hat <- exp(X %*% gmmo$coefficients)
  residuals <- as.vector(y - lambda_hat)
  return(residuals)
}
