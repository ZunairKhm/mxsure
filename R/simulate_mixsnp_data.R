#' Simulate Mixed related SNP distance dataset
#'
#' @param lambda mutation rate input (SNPs per day)
#' @param k proportion of mixed dataset that is linked, set to 0 to produce a distant dataset
#' @param error_param chance that each SNP is missed
#' @param n size of dataset to simulate
#' @param unrel_rate rate of gamma distribution for time of evolution unrelated pairs in years
#' @param unrel_shape shape of gamma distribution for time of evolution for unrelated pairs in years
#' @param rel_timelimit maximum time allowed between samples in years (calculates time difference from uniform dist.)
#' @param rel_timemin minimum time allowed between samples in years
#' @param truncation_point maximum allowed SNP distances, if set too low can cause all SNP distances to be 0
#'
#' @return SNP distance dataset with a mixture of related and unrelated pairs with time differences
#' @export
#'
#' @examples
simulate_mixsnp_data <- function(lambda, k, unrel_shape=100, unrel_rate=1, error_param=NA, n=100, rel_timelimit=1, rel_timemin=0, truncation_point=NA){
  mix_snp_dist <- map_dfr(1:n, ~{

    if(is.na(truncation_point)){
      truncation_point <- Inf
    }

    tt <- runif(1, rel_timemin*365.25, rel_timelimit*365.25) #rexp(1, rate = -log(1-0.99)/(reltimelimit*365.25)) #time distribution for sample time different
    td <- rgamma(1, unrel_shape,rate= unrel_rate)*365.25 #abs(rnorm(1, mean=unrelmu*365.25, sd = ifelse(anyNA(unrelsd),unrelmu*365.25, unrelsd*365.25))) #time distribution for unrelated evo time
    if (runif(1) <k) {
      dd <- qpois(runif(1)*ppois(truncation_point, tt*lambda), tt*lambda) #rpois(n = 1, lambda = tt*lambda)
      rr <- "Related"
    } else {
      dd <- qpois(runif(1)*ppois(truncation_point, td*lambda), td*lambda) #rpois(1, lambda = td*lambda)
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



