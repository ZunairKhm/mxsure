#' Bootstrapped confidence intervals for mixture SNP cutoffs
#'
#' Creates confidence intervals for outputs from mxsure_estimate using bootstrapping. Utilises the furrr package to allow for parallel computation, to enable use the plan() function from the future package.
#'
#' @param mixed_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param mixed_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param mixed_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_size size of each bootstrap sample
#' @param sample_n number of bootstrap sampling to conduct
#' @param start_params initial parametrs for optim, if NA (as default) will try 3 different start parameters and produce the highest likelyhood result. Specifying the start parameters minimises computing time.
#' @param truncation_point a SNP distance limit for the data, if set to NA will estimate as if there is no limit
#' @param confidence_level confidence level to produce confidence intervals
#' @param lambda_bounds bounds of rate estimation in SNPs/year/site if given time and site data
#' @param k_bounds bounds of related proportion estimation
#' @param intercept_bounds bounds of intercept estimation
#'
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom dplyr bind_rows summarise across everything
#' @importFrom tibble tibble
#' @importFrom stats complete.cases
#'
#' @return Confidence intervals for all estimates produced by mxsure_estimate.
#'
#' @export
mxsure_ci <- function(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist=NA, mixed_sites=NA,
                      sample_size=length(mixed_snp_dist),truncation_point=NA, sample_n=100, confidence_level=0.95, start_params="Efficient",
                      tree=NA, sampleA=NA, sampleB=NA,
                      lambda_bounds = c(0, 1), k_bounds=c(0,1), intercept_bounds=c(-Inf, Inf), single_branch_lambda_bounds = c(0, Inf), branch_intercept=NA,
                      quiet=FALSE){

  snp_dist <-NULL

  if(is.na(truncation_point)){
    truncation_point <- Inf
  }

  if (sample_n==0){
    lowerres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)
    upperres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)

    ci <- bind_rows(lowerres, upperres)
    res <- list(confidence_intervals=ci, raw_results=NA)
  }
  else if ((length(mixed_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){

  mix_data <- tibble(snp_dist=mixed_snp_dist, time_dist=mixed_time_dist, sites=mixed_sites)

  mix_data$sampleA <- sampleA
  mix_data$sampleB <- sampleB

  mix_data <- filter(mix_data, snp_dist<truncation_point)
  unrelated_snp_dist <- unrelated_snp_dist[unrelated_snp_dist<truncation_point]



  if (anyNA(start_params)){
    start_params <- NA
    } else if(all(start_params=="Efficient")){
  test_result <- #suppressWarnings(
      mxsure_estimate(
        mix_data$snp_dist,unrelated_snp_dist, mix_data$time_dist, mix_data$sites, truncation_point=truncation_point, start_params = NA,
        tree=tree, sampleA=mix_data$sampleA, sampleB=mix_data$sampleB, branch_intercept=branch_intercept,
        lambda_bounds = lambda_bounds, k_bounds=k_bounds,  intercept_bounds=intercept_bounds,
        single_branch_lambda_bounds = single_branch_lambda_bounds)
  #, classes = "warning")
  start_params <- as.numeric(c(test_result[3], test_result[2], test_result[4], test_result[7], test_result[8]))
    }


  raw_data <- list()
  #bootstrapping both close and distant data sets allowing for parallelisation
  raw_bootstrapresults <- furrr::future_map_dfr(1:sample_n, ~{
    x <- slice_sample(mix_data, n= nrow(mix_data), replace = TRUE)
    y <- sample(unrelated_snp_dist, size=length(unrelated_snp_dist), replace=TRUE)
    z <- plyr::try_default(
  suppressWarnings(
    mxsure_estimate(
      x$snp_dist, y, x$time_dist, x$sites, truncation_point=truncation_point,
      tree=tree, sampleA=x$sampleA, sampleB=x$sampleB, branch_intercept=branch_intercept,
      start_params = start_params, lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds, single_branch_lambda_bounds = single_branch_lambda_bounds
    )
  , classes = "warning")
  ,data.frame(snp_threshold=NA, lambda=NA, k=NA,intercept=NA, estimated_fp=NA))

  return(z)

  },.progress=!quiet, .options = furrr::furrr_options(seed = TRUE))


  bootstrapresults <- raw_bootstrapresults[complete.cases(raw_bootstrapresults), ] #removes failed results from optimisation failures
  bootstrapresults <- bootstrapresults[1:5]

  lowerres <- bootstrapresults|> #finds quantiles for confidence intervals
    summarise(across(everything(),  ~quantile(.x, 1-confidence_level, na.rm=TRUE)))
  upperres <- bootstrapresults|>
    summarise(across(everything(), ~quantile(.x, confidence_level, na.rm=TRUE)))

  ci <- bind_rows(lowerres, upperres)
  res <- list(confidence_intervals=ci, raw_results=raw_bootstrapresults)

  } else {
    warning("Insufficient data points to fit distributions!")

    lowerres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,intercept=NA,estimated_fp=NA)
    upperres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,intercept=NA,estimated_fp=NA)

    ci <- bind_rows(lowerres, upperres)
    res <- list(confidence_intervals=ci, raw_results=NA)
  }
}
