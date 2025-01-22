#' Simulate Mixed related SNP distance dataset
#'
#' @param lambda mutation rate input (SNPs per day)
#' @param k proportion of mixed dataset that is linked
#' @param nbmu mean of the negative binomial distribution underlying the unrelated SNP distances
#' @param error_param chance that each SNP is missed
#' @param n size of dataset to simulate
#'
#' @return SNP distance dataset with a mixture of related and unrelated pairs with time differences
#' @export
#'
#' @examples
simulate_mixsnp_data <- function(lambda, k,unrelmu=100*365.25, unrelsd=30*365.25, error_param=NA, n=1000){
  mix_snp_dist <- map_dfr(1:n, ~{

    tt <- rexp(1, rate = 0.02) #time distribution for sample time different
    td <- abs(rnorm(1, mean=unrelmu, sd=unrelsd)) #time distribution for unrelated evo time
    if (runif(1)<k){
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
