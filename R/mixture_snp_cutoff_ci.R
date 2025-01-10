#' Bootstrapped confidence intervals for mixture SNP cutoffs
#'
#' Creates confidence intervals for outputs from mixture_snp_cutoff using bootstrapping. Utilises the furrr package to allow for parallel computation, to enable use the plan() function from the future package.
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_size size of each bootstrap sample
#' @param sample_n number of bootstrap sampling to conduct
#'
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom diptest dip.test
#'
#' @return Confidence intervals
#'
#' @export
mixture_snp_cutoff_ci <- function(trans_snp_dist,unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                                  sample_size=length(trans_snp_dist), sample_n=1000, confidence_level=0.95){
  if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){

  mix_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)

  bootstrapresults <- furrr::future_map_dfr(1:sample_n, ~{
    x <- slice_sample(mix_data, n= sample_size, replace = TRUE)
    y <- plyr::try_default(
      suppressWarnings(
        mixture_snp_cutoff(
          x$snp_dist,unrelated_snp_dist, x$time_dist, x$sites
          ), classes = "warning"),
                           data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA))
   y[1:4]
  },.progress=TRUE,.options = furrr::furrr_options(seed = TRUE))

  bootstrapresults <- bootstrapresults[complete.cases(bootstrapresults), ]

  if(!anyNA(bootstrapresults)){
    p <- diptest::dip.test(bootstrapresults$lambda)
    dip_pvalue <- p$p.value
    if(p$p.value<0.05){
      warning("Bootstrap samples indicate multimodal distribution: estimates may be unreliable")
    }}
  else {dip_pvalue <- NA}

  lowerres <- bootstrapresults|>
    summarise(across(everything(),  ~quantile(.x, 1-confidence_level, na.rm=TRUE)))
  upperres <- bootstrapresults|>
    summarise(across(everything(), ~quantile(.x, confidence_level, na.rm=TRUE)))

  ci <- bind_rows(lowerres, upperres)
  res <- list(confidence_intervals=ci, raw_results=bootstrapresults, dip_pvalue=dip_pvalue)
  # return(res)
  } else {
    warning("Insufficient data points to fit distributions!")

    lowerres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)
    upperres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)

    ci <- bind_rows(lowerres, upperres)
    res <- list(confidence_intervals=ci, raw_results=NA, dip_pvalue=NA)
  }
}
