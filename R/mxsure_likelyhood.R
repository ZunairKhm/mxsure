#' Reports Likelyhoods for related and unrelated distributions for each datapoint provided
#'
#'
#'
#' @param mixed_snp_dist vector of SNP distances from mixed dataset
#' @param unrelated_snp_dist vector of SNP distances from unrelated dataset
#' @param mixed_time_dist vector of time differences for each SNP distacne in the mixed dataset
#' @param mixed_sites vector of sites considered for each SNP distance in the mixed dataset
#' @param truncation_point maximum limit of SNP distances to consider
#'
#' @importFrom dplyr distinct mutate
#' @importFrom tibble tibble
#' @importFrom stats dpois dnbinom
#'
#' @return a dataframe with SNP distances, time differences, sites considered and the likeyhoods of related and unrelated models fitting for each datapoint
#' @export
mxsure_likelyhood <- function(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, truncation_point=2000){

  unrelated_snp_dist <- unrelated_snp_dist[unrelated_snp_dist<truncation_point]
  mix_res <- suppressWarnings(mxsure_estimate(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, truncation_point = 2000))

  if(is.na(mean(mixed_sites, na.rm=TRUE))){
    mixed_sites <- 1
  }

  LH <-  distinct(tibble(snp_dist=mixed_snp_dist, time_dist=mixed_time_dist))
  LH <- LH |>
    mutate(rel_loglh = (dpois(snp_dist, (mix_res$lambda*(time_dist/365.25)*mean(mixed_sites)+mix_res$intercept), log=TRUE) #/ ppois(truncation_point, (mix_res$lambda*(time_dist/365.25)*mean(mixed_sites)+mix_res$intercept))
    ),
    unrel_loglh = (dnbinom(snp_dist, mu = mix_res$nb_mu, size = mix_res$nb_size, log = TRUE) - pnbinom(truncation_point, mu = mix_res$nb_mu, size = mix_res$nb_size, log.p=TRUE)))|>
    mutate(logLHR = (rel_loglh-unrel_loglh))|>
    mutate(rel_lh=exp(rel_loglh),
           unrel_lh=exp(unrel_loglh),
           LHR=exp(logLHR))


  return(LH)

}
