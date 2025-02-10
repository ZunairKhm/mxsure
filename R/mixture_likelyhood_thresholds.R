#' Mixture Likelyhood Thresholds
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param year_range range of years to calculate thresholds over
#'
#' @return Threshold corresponding to the highest SNP distance where it is more likely to be related than unrelated
#' @export
#'
#' @examples
mixture_likelyhood_thresholds <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites, year_range=c(seq(0.5, 10, 0.5))){
  if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){
  nb_fit <- suppressWarnings(MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial"))
  mix_res <- suppressWarnings(mixture_snp_cutoff(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites))


    threshold <- map_dbl(year_range, function(year){
    x <- 1:500
    rel_prob <- mix_res$k * dpois(x, (mix_res$lambda*year*365.25*mean(trans_sites))/1000000)
    unrel_prob <- (1-mix_res$k) * dnbinom(x, mu = nb_fit$estimate['mu'], size = nb_fit$estimate['size'])
    thresh <- x[max(which(rel_prob > unrel_prob))]



 })
    result <- tibble(year_range, threshold, estimated_fp=NA, prop_pos=NA)
      result$estimated_fp <- modify(result$threshold, ~{sum(unrelated_snp_dist<=.x)/length(unrelated_snp_dist)})
      result$prop_pos <- modify(result$threshold,  ~{sum(trans_snp_dist<=.x)/length(trans_snp_dist)})


  return(result)


  }else {
    warning("Insufficient data points to fit distributions!")
    return(tibble(NA))
  }
}
