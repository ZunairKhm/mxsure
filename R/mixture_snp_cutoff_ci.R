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
#' @param start_params initial parametrs for optim, if NA (as default) will try 3 different start parameters and produce the highest likelyhood result. Specifying the start parameters minimises computing time.
#' @param truncation_point a SNP distance limit for the data, if set to NA will estimate as if there is no limit
#' @param confidence_level confidence level to produce confidence intervals
#'
#' @importFrom furrr future_map_dfr furrr_options
#'
#' @return Confidence intervals
#'
#' @export
mixture_snp_cutoff_ci <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                                  sample_size=length(trans_snp_dist),truncation_point=NA, sample_n=1000, confidence_level=0.95, start_params=NA){

  if(is.na(truncation_point)){
    truncation_point <- Inf
  }

  if (sample_n==0){
    lowerres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)
    upperres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)

    ci <- bind_rows(lowerres, upperres)
    res <- list(confidence_intervals=ci, raw_results=NA)
  }
  else if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){

  mix_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)

  if (anyNA(start_params)){
  test_result <- suppressWarnings(
      mixture_snp_cutoff(
        mix_data$snp_dist,unrelated_snp_dist, mix_data$time_dist, mix_data$sites,truncation_point=truncation_point, start_params = NA
      ), classes = "warning")
  start_params <- c(test_result[3], test_result[2])
  }

  mix_data <- filter(mix_data, snp_dist<truncation_point)
  unrelated_snp_dist <- unrelated_snp_dist[unrelated_snp_dist<truncation_point]

  #bootstrapping both close and distant data sets allowing for parallelisationg
  bootstrapresults <- furrr::future_map_dfr(1:sample_n, ~{
    x <- slice_sample(mix_data, n= sample_size, replace = TRUE)
    y <- sample(unrelated_snp_dist, size=length(unrelated_snp_dist), replace=TRUE)
    z <- plyr::try_default( #suppresses warnings and errors from mixture_snp_cutoffs
      suppressWarnings(
        mixture_snp_cutoff(
          x$snp_dist,y, x$time_dist, x$sites, truncation_point=truncation_point, start_params = start_params
          ), classes = "warning"),
                           data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA))
   z[1:4]
  },.progress=TRUE, .options = furrr::furrr_options(seed = TRUE))

  bootstrapresults <- bootstrapresults[complete.cases(bootstrapresults), ] #removes failed results from optim failures

  lowerres <- bootstrapresults|> #finds quantiles for confidence intervals
    summarise(across(everything(),  ~quantile(.x, 1-confidence_level, na.rm=TRUE)))
  upperres <- bootstrapresults|>
    summarise(across(everything(), ~quantile(.x, confidence_level, na.rm=TRUE)))

  ci <- bind_rows(lowerres, upperres)
  res <- list(confidence_intervals=ci, raw_results=bootstrapresults)
  } else {
    warning("Insufficient data points to fit distributions!")

    lowerres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)
    upperres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)

    ci <- bind_rows(lowerres, upperres)
    res <- list(confidence_intervals=ci, raw_results=NA)
  }
}
