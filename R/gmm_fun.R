#' Auxiliary functions for GMM estimation - Poisson
#'
#' @param theta numeric; Vetor of parameter values, see \code{?gmm::gmm}.
#' @param data matrix; Data matrix for all X, see \code{?gmm::gmm}.
#'
#' @returns Returns moment the moment conditions
#'

poisson_moments <- function(theta, data) {
  Y <- data$y
  X <- data$x
  lambda_hat <- exp(X %*% theta)
  residuals <- as.vector(Y - lambda_hat)
  g <- residuals * X
  return(g)
}

#' Auxiliary functions for GMM estimation - Logistic
#'
#' @param theta numeric; Vetor of parameter values, see \code{?gmm::gmm}.
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
