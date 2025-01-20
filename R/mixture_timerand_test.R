#' Time Randomisation Test for mixture distribution mutation rate estimation
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_size size of each bootstrap sample
#' @param sample_n number of bootstrap sampling to conduct
#' @param start_params initial parameters for optim, if NA (as default) will try 3 different start parameters and produce the highest likelyhood result. Specifying the start parameters minimises computing time.
#' @param ci_data optional input for previously calculated CI data (mixture_snp_cutoff_ci) for computational efficiency
#' @param raw option to output raw results instead of ggplot
#' @param confidence_level confidence level for CIs
#' @param title titel for ggplot
#'
#' @return ggplot comparing point estimates and confidence levels between normal data and time randomised data
#' @export
#'
#' @examples
mixture_timerand_test <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                  sample_size=length(trans_snp_dist), sample_n=500, confidence_level=0.95, start_params=NA, ci_data=NA, raw=FALSE, title=NULL){

  normal_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
  timerand_data <- tibble(snp_dist=trans_snp_dist, time_dist=sample(trans_time_dist, length(trans_time_dist)), sites=trans_sites)


  normal_result <- mixture_snp_cutoff(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites)

  timerand_result <- mixture_snp_cutoff(timerand_data$snp_dist,unrelated_snp_dist, timerand_data$time_dist,timerand_data$sites)

  if(anyNA(ci_data)){
    normal_ci <- mixture_snp_cutoff_ci(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites,
                                      sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level,
                                      start_params = c(normal_result[3], normal_result[2]))
  } else{
    normal_ci <- ci_data
  }

  timerand_ci <- mixture_snp_cutoff_ci(timerand_data$snp_dist,unrelated_snp_dist, timerand_data$time_dist,timerand_data$sites,
                                       sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level,
                                       start_params = c(normal_result[3], normal_result[2]))
  result <- tibble(
    "method"=c("Normal", "Time Randomised"),
    "5%"=c("normal" = normal_ci$confidence_intervals$lambda[1], "time_randomised" = timerand_ci$confidence_intervals$lambda[1], use.names=FALSE),
    "point_est"=c("normal" = normal_result$lambda, "time_randomised" = timerand_result$lambda, use.names=FALSE),
    "95%"=c("normal" = normal_ci$confidence_intervals$lambda[2], "time_randomised" = timerand_ci$confidence_intervals$lambda[2], use.names=FALSE)
  )
  xres_timerand <- list(
    result=result,
    normal_result=normal_result,
    timerand_result= timerand_result,
    normal_ci=normal_ci,
    timerand_ci=timerand_ci,
    normal_data=normal_data,
    timerand_data=timerand_data
  )

  if(raw==TRUE){
  return(xres_timerand)
  } else{
    ggplot(xres_timerand$result, aes(x = method, y = point_est)) +
      geom_point(data=tibble("point_est"=xres_timerand$normal_ci$raw_results$lambda, "method"="Normal"),color="grey60",size=1, alpha=0.5)+
      geom_point(data=tibble("point_est"=xres_timerand$timerand_ci$raw_results$lambda, "method"="Time Randomised"),color="grey60",size=1, alpha=0.5)+
      geom_point(size = 2, color = "red3") +  # Point estimate
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
      labs(title=title,
        x = "Method",
           y = "Rate") +
      theme_minimal()
  }
}
