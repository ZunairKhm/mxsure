#' Right truncated negative binomial distribution
#'
#' Density, distribution function, and quantile function for the right truncated negative binomial distribution with parameters size, prob, and truncation point.
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param mu mean of distribution
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param right_truncation right censoring point
#'
#' @name truncnbinom
NULL
#> NULL


#' @rdname truncnbinom
#' @importFrom stats dnbinom pnbinom
#' @export
dtruncnbinom <- function(x, mu, size, right_truncation){
  dnbinom(x = x, size = size, mu = mu) / (pnbinom(right_truncation, size = size, mu = mu))
}

#' @rdname truncnbinom
#' @importFrom stats pnbinom
#' @export
ptruncnbinom <- function(q, mu, size, right_truncation){
  pnbinom(q = q, size = size, mu = mu) / (pnbinom(right_truncation, size = size, mu = mu))
}

#' @rdname truncnbinom
#' @importFrom stats qnbinom
#' @export
qtruncnbinom <- function(p, mu, size, right_truncation){
  qnbinom(p = p*(pnbinom(right_truncation, size=size, mu=mu)), size = size, mu = mu)
}
