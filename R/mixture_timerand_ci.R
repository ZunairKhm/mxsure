#' Time Randomisation Test for mixture distribution mutation rate estimation
#
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#'
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_size size of each bootstrap sample
#' @param sample_n number of bootstrap sampling to conduct
#' @param start_params initial parameters for optim, if NA (as default) will try 3 different start parameters and produce the highest likelihood result. Specifying the start parameters minimises computing time.
#' @param ci_data optional input for previously calculated CI data (mixture_snp_cutoff_ci) for computational efficiency
#' @param confidence_level confidence level for CIs
#' @param title title for ggplot
#' @param truncation_point SNP distances to truncate at
#' @param permutations number of time permutation to run
#' @param quiet if false will print progress bar for each permutation
#' @param within_individual permute time data within each individuals (requires indiviudal ID codes that must be the exact same)
#' @param subjectA_id needed for within individual permutation
#' @param subjectB_id needed for within individual permutation (must be the same as subject A)
#' @param clustered permute between different time distances
#' @param p_value_type either "above_estimate", "above_low_ci", or "within_ci"
#'
#' @return ggplot comparing point estimates and confidence levels between normal data and time randomised data
#' @export
#'
#' @examples
mixture_timerand_ci <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA, truncation_point=NA,
                                  sample_size=length(trans_snp_dist), sample_n=500, permutations=3, quiet=FALSE, confidence_level=0.95,
                                  within_individual=FALSE, subjectA_id, subjectB_id,
                                  clustered=FALSE, p_value_type="above_estimate",
                                  start_params=NA, ci_data=NA, title=NULL){
  #unadjusted result
  normal_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
  normal_result <- mixture_snp_cutoff(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites, truncation_point=truncation_point)
  if(anyNA(ci_data)){
    normal_ci <- mixture_snp_cutoff_ci(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites, truncation_point=truncation_point,
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
    if(!quiet){print(paste0("Processing Permutation: ", i))}
    if(within_individual){
      if(all(subjectA_id==subjectB_id)){
        timerand_data <- tibble(subject_id= subjectA_id, snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
        timerand_data <- timerand_data|>
          group_by(subject_id)|>
          mutate(time_dist=sample(time_dist, length(time_dist)))|>
          ungroup()
      } else{ warning("subject ID's do not match")}
    }else if (clustered){
      distinct_time_dist <- tibble(time_dist=trans_time_dist)
      distinct_time_dist$permuted <- sample(distinct_time_dist$time_dist)
      timerand_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
      timerand_data$time_dist <- distinct_time_dist$permuted[match(timerand_data$time_dist, distinct_time_dist$time_dist)]
    }else {
      timerand_data <- tibble(snp_dist=trans_snp_dist, time_dist=sample(trans_time_dist, length(trans_time_dist)), sites=trans_sites)
    }

  timerand_result <- mixture_snp_cutoff(timerand_data$snp_dist,unrelated_snp_dist, timerand_data$time_dist, timerand_data$sites, truncation_point=truncation_point)
  timerand_ci <- mixture_snp_cutoff_ci(timerand_data$snp_dist,unrelated_snp_dist, timerand_data$time_dist,timerand_data$sites,
                                       sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level, truncation_point=truncation_point,
                                       start_params = c(normal_result[3], normal_result[2]))
  if(p_value_type=="above_estimate"){
  p_value_n <- sum(timerand_ci$raw_results$lambda>=normal_result$lambda)
  }else if(p_value_type=="within_ci"){
    p_value_n <- sum(timerand_ci$raw_results$lambda<=result$`95%`[result$method=="Normal"]&timerand_ci$raw_results$lambda>=result$`5%`[result$method=="Normal"])
  }else if (p_value_type=="above_low_ci"){
    p_value_n <- sum(timerand_ci$raw_results$lambda>=result$`5%`[result$method=="Normal"])
  }
    p_value_t <- length(timerand_ci$raw_results$lambda)

  # Append results to 'result'
  result <- bind_rows(result, tibble(
    method = paste0("TR ", i),
    `5%` = timerand_ci$confidence_intervals$lambda[1],
    point_est = timerand_result$lambda,
    `95%` = timerand_ci$confidence_intervals$lambda[2],
    p_value_n= p_value_n,
    p_value_t=p_value_t
  ))

  # Append raw results with a method column
  rawtimerand <- bind_rows(rawtimerand,
                           timerand_ci$raw_results %>% mutate(method = paste0("TR ", i)))
  }

  result$method <- factor(result$method, levels = result$method)
  rawtimerand$method <- factor(rawtimerand$method, levels = result$method)





  outcome <- tibble(
    n_permutations=permutations,
    n_overlapping_est=sum(result$`95%`[result$method!="Normal"]>=result$point_est[1]),
    perc_overlapping_est=sum(result$`95%`[result$method!="Normal"]>=result$point_est[1])/length(result$`95%`[result$method!="Normal"]),
    n_overlapping_lowci=sum(result$`95%`[result$method!="Normal"]>=result$`5%`[1]),
    perc_overlapping_lowci=sum(result$`95%`[result$method!="Normal"]>=result$`5%`[1])/length(result$`95%`[result$method!="Normal"]),
    p_value=(sum(result$p_value_n, na.rm=TRUE)+1)/(sum(result$p_value_t, na.rm=TRUE)+1)
  )

    plot <- ggplot(result, aes(x = method, y = point_est)) +
      geom_hline(yintercept = result$point_est[1], color="grey60")+
      geom_point(data=rawtimerand,aes(x=method, y=lambda),color="grey50",size=0.8, alpha=0.3)+
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
      geom_point(size = 2, color = "red3") +  # Point estimate
      scale_y_continuous(#transform = "log10"
                         )+
      annotate("label", label=paste0("p=",format(round(outcome$p_value, 4), nsmall = 4)), x=Inf, y=Inf, vjust=1, hjust=1)+
      labs(title=title,
           x = NULL,
           y = "Rate") +
      theme_minimal()




  xres_timerand <- list(
    result=result,
    raw_results=rawtimerand,
    outcome=outcome,
    plot=plot
  )
}
