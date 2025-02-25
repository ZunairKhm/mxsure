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
    tt <- runif(1, rel_timemin*365.25, rel_timelimit*365.25)
    td <- rgamma(1, unrel_shape,rate= unrel_rate)*365.25
    if (runif(1) <k) {
      truncation_correction <- ifelse(ppois(truncation_point, tt*lambda)!=0, ppois(truncation_point, tt*lambda), NA)
      dd <- qpois(runif(1)*truncation_correction, tt*lambda)
      rr <- "Related"
    } else {
      truncation_correction <- ifelse(ppois(truncation_point, td*lambda)!=0, ppois(truncation_point, td*lambda), NA)
      dd <- qpois(runif(1)*truncation_correction, td*lambda)
      rr <- "Unrelated"
    }
    if (!is.na(error_param)){
      dd <- rbinom(1, dd, 1-error_param) #error param is chance each SNP difference is missed
    }
    return(tibble(snp_dist=dd, time_dist=tt, relation=rr))
  })

  mix_snp_dist <-  na.omit(mix_snp_dist)
    if(nrow(mix_snp_dist)==0){
      warning("No data could be simulated due to the cumulative density at the truncation point for all simulated evolution times being 0")
    } else if(nrow(mix_snp_dist)<n){
      warning(paste0(n-nrow(mix_snp_dist) ," data point(s) could not be simulated due to the cumulative density at the truncation point for the corresponding simulated evolution times being 0"))
    }
    return(mix_snp_dist)
}



