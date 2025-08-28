#' Right truncated negative binomial distribution
#'
#' Density, distribution function, and quantile function for the right truncated negative binomial distribution with parameters size, prob, and truncation point.
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param mu mean of distribution
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param truncation_point right censoring point
#'
#' @name truncnbinom
NULL
#> NULL


#' @rdname truncnbinom
#' @importFrom stats dnbinom pnbinom
#' @export
dtruncnbinom <- function(x, mu, size, truncation_point){
  dnbinom(x = x, size = size, mu = mu) / (pnbinom(truncation_point, size = size, mu = mu))
}

#' @rdname truncnbinom
#' @importFrom stats pnbinom
#' @export
ptruncnbinom <- function(q, mu, size, truncation_point){
  pnbinom(q = q, size = size, mu = mu) / (pnbinom(truncation_point, size = size, mu = mu))
}

#' @rdname truncnbinom
#' @importFrom stats qnbinom
#' @export
qtruncnbinom <- function(p, mu, size, truncation_point){
  qnbinom(p = p*(pnbinom(truncation_point, size=size, mu=mu)), size = size, mu = mu)
}
