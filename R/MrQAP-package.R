#' @keywords internal
"_PACKAGE"

# The S3 generics these methods extend must be in the namespace, or
# registerS3methods() fails when the installed package is loaded.
#' @importFrom stats AIC BIC as.formula coef coefficients confint cor fitted
#' @importFrom stats glm lm logLik na.omit nobs poisson qnorm quantile rbinom
#' @importFrom stats relevel residuals rnorm sd setNames vcov
#' @importFrom stats predict binomial glm.fit formula
#' @importFrom graphics hist abline box
#' @importFrom utils head
NULL
