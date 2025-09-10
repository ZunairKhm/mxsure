#' Simulate mixed related SNP distance dataset
#'
#' @param lambda mutation rate input (SNPs per year)
#' @param k proportion of mixed dataset that is linked, set to 0 to produce a distant dataset
#' @param error_param chance that each SNP is missed. For simulating systematic error.
#' @param n size of dataset to simulate
#' @param rel_timemax maximum time allowed between samples in years (calculates time difference from uniform dist.)
#' @param rel_timemin minimum time allowed between samples in years
#' @param right_truncation maximum allowed SNP distances, if set too low can cause all SNP distances to be 0
#' @param unrel_mean mean of gamma distribution for which the evolutionary time for the unrelated datapoints are chosen
#' @param unrel_sd sd of gamma distribution for which the evolutionary time for the unrelated datapoints are chosen
#'
#' @importFrom purrr map_dfr
#' @importFrom stats na.omit
#'
#' @return SNP distance dataset with a mixture of related and unrelated pairs with time differences
#' @export
#'
#' @examples
#' simulate_mixsnp_data(1, 0.8)
#'
simulate_mixsnp_data <- function(lambda, k, unrel_mean=25, unrel_sd=12.5, error_param=NA, n=100, rel_timemax=1, rel_timemin=0, right_truncation=NA){

  # n <- n/dpois(right_truncation, lambda*(unrel_shape/unrel_rate))
  unrel_shape <- (unrel_mean/unrel_sd)^2
  unrel_rate <- unrel_mean/(unrel_sd^2)

  mix_snp_dist <- map_dfr(1:n, ~{

    if(is.na(right_truncation)){
      right_truncation <- Inf
    }
    tt <- runif(1, rel_timemin*365.25, rel_timemax*365.25)
    td <- rgamma(1, unrel_shape,rate= unrel_rate)*365.25
    if (runif(1) <k) {
      truncation_correction <- ifelse(ppois(right_truncation, tt*lambda/365.25)>(1/n), ppois(right_truncation, tt*lambda/365.25), NA)
      dd <- qpois(runif(1)*truncation_correction, tt*lambda/365.25)
      rr <- "Related"
    } else {
      truncation_correction <- ifelse(ppois(right_truncation, td*lambda/365.25)>(1/n), ppois(right_truncation, td*lambda/365.25), NA)
      dd <- qpois(runif(1)*truncation_correction, td*lambda/365.25)
      rr <- "Unrelated"
    }
    if (!is.na(error_param)){
      dd <- rbinom(1, dd, 1-error_param) #error param is chance each SNP difference is missed
    }
    return(tibble(snp_dist=dd, time_dist=tt, relation=rr))
  })

  # mix_snp_dist <- filter(mix_snp_dist, snp_dist<=right_truncation)

  mix_snp_dist <-  na.omit(mix_snp_dist)
    if(nrow(mix_snp_dist)==0){
      warning("No data could be simulated due to the cumulative density at the truncation point for all simulated evolution times being lower than 1/n")
    } else if(nrow(mix_snp_dist)<n){
      warning(paste0(n-nrow(mix_snp_dist) ," data point(s) could not be simulated due to the cumulative density at the truncation point for the corresponding simulated evolution times being lower than 1/n"))
    }
    return(mix_snp_dist)
}



