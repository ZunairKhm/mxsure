#' Time Randomisation Test for mixture distribution mutation rate estimation
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_size size of each bootstrap sample
#' @param sample_n number of bootstrap sampling to conduct
#' @param start_params initial parameters for optim, if NA (as default) will try 3 different start parameters and produce the highest likelihood result. Specifying the start parameters minimises computing time.
#' @param ci_data optional input for previously calculated CI data (mixture_snp_cutoff_ci) for computational efficiency
#' @param confidence_level confidence level for CIs
#' @param title title for ggplot
#'
#' @return ggplot comparing point estimates and confidence levels between normal data and time randomised data
#' @export
#'
#' @examples
mixture_timerand_ci <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                                  sample_size=length(trans_snp_dist), sample_n=500, permutations=3, confidence_level=0.95,
                                  start_params=NA, ci_data=NA,  raw=FALSE, title=NULL){
  #unadjusted result
  normal_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
  normal_result <- mixture_snp_cutoff(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites)
  if(anyNA(ci_data)){
    normal_ci <- mixture_snp_cutoff_ci(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites,
                                       sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level,
                                       start_params = c(normal_result[3], normal_result[2]))
  } else{
    normal_ci <- ci_data
  }

  result <- tibble(
    "method"="Normal",
    "5%"=normal_ci$confidence_intervals$lambda[1],
    "point_est"=normal_result$lambda,
    "95%"=normal_ci$confidence_intervals$lambda[2]
  )

  # time rand loop
  rawtimerand <- normal_ci$raw_results |>
    mutate(method = "Normal")
  for(i in 1:permutations){
  timerand_data <- tibble(snp_dist=trans_snp_dist, time_dist=sample(trans_time_dist, length(trans_time_dist)), sites=trans_sites)
  timerand_result <- mixture_snp_cutoff(timerand_data$snp_dist,unrelated_snp_dist, timerand_data$time_dist,timerand_data$sites)
  timerand_ci <- mixture_snp_cutoff_ci(timerand_data$snp_dist,unrelated_snp_dist, timerand_data$time_dist,timerand_data$sites,
                                       sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level,
                                       start_params = c(normal_result[3], normal_result[2]))
  # Append results to 'result'
  result <- bind_rows(result, tibble(
    method = paste0("Time Randomisation ", i),
    `5%` = timerand_ci$confidence_intervals$lambda[1],
    point_est = timerand_result$lambda,
    `95%` = timerand_ci$confidence_intervals$lambda[2]
  ))

  # Append raw results with a method column
  rawtimerand <- bind_rows(rawtimerand,
                           timerand_ci$raw_results %>% mutate(method = paste0("Time Randomisation ", i)))
  }



    plot <- ggplot(result, aes(x = method, y = point_est)) +
      geom_hline(yintercept = result$point_est[1], color="grey60")+
      geom_point(data=rawtimerand,aes(x=method, y=lambda),color="grey50",size=0.8, alpha=0.3)+
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
      geom_point(size = 2, color = "red3") +  # Point estimate
      scale_y_continuous(#transform = "log10"
                         )+
      labs(title=title,
           x = "Method",
           y = "Rate") +
      theme_minimal()


  xres_timerand <- list(
    result=result,
    rawtimerand=rawtimerand,
    plot=plot
  )
}
