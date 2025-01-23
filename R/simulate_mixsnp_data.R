#' Simulate Mixed related SNP distance dataset
#'
#' @param lambda mutation rate input (SNPs per day)
#' @param k proportion of mixed dataset that is linked
#' @param error_param chance that each SNP is missed
#' @param n size of dataset to simulate
#' @param unrelmu mean of normal distribution for SNP distances for unrelated pairs in years
#' @param unrelsd sd of normal distribution for SNP distances for unrelated pairs in years
#' @param reltimelimit maximum time allowed between samples (calculates a rate for an exponential distribution to model this time difference as)
#'
#' @return SNP distance dataset with a mixture of related and unrelated pairs with time differences
#' @export
#'
#' @examples
simulate_mixsnp_data <- function(lambda, k,unrelmu=100, unrelsd=NA, error_param=NA, n=100, reltimelimit=c(1)){
  mix_snp_dist <- map_dfr(1:n, ~{

    tt <- rexp(1, rate = -log(1-0.99)/(reltimelimit*365.25)) #time distribution for sample time different
    td <- abs(rnorm(1, mean=unrelmu*365.25, sd = ifelse(anyNA(unrelsd),unrelmu*365.25, unrelsd*365.25))) #time distribution for unrelated evo time
    if (runif(1) <k) {
      dd <- rpois(n = 1, lambda = tt*lambda)
      rr <- "Related"
    } else {
      dd <- rpois(1, lambda = td*lambda)
      rr <- "Unrelated"
    }
    if (!is.na(error_param)){
      dd <- rbinom(1, dd, 1-error_param) #error param is chance each SNP difference is missed
    }
    return(
      tibble(snp_dist=dd, time_dist=tt, relation=rr)
    )

  })
}
