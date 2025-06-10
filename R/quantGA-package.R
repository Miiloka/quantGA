#' quantGA: Fits a linear quantile regression model using genetic algorithm
#'
#' Provides tools to fit linear quantile regression models using a genetic
#' algorithm to estimate the coefficients.
#'
#' @docType package
#' @name quantGA
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline legend lines
#' @importFrom stats IQR dnorm model.frame model.matrix model.response pnorm printCoefmat rnorm runif sd var
#' @useDynLib quantGA, .registration = TRUE
## usethis namespace: end
NULL

