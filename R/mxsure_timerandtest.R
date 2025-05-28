#' Time Randomisation Test for mixture distribution mutation rate estimation
#
#' @param mixed_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param mixed_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param mixed_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_size size of each bootstrap sample
#' @param sample_n number of bootstrap sampling to conduct for each time permutation
#' @param start_params initial parameters for optim, if NA (as default) will try 3 different start parameters and produce the highest likelihood result. Specifying the start parameters minimises computing time.
#' @param ci_data optional input for previously calculated CI data (mxsure_ci) for computational efficiency
#' @param confidence_level confidence level for CIs
#' @param title title for ggplot
#' @param truncation_point SNP distances to truncate at
#' @param permutations number of time permutation to run
#' @param quiet if false will print progress bar for each permutation
#' @param within_individual permute time data within each individuals (requires indiviudal ID codes that must be the exact same)
#' @param subjectA_id needed for within individual permutation
#' @param subjectB_id needed for within individual permutation (must be the same as subject A)
#' @param clustered permute between different time distances
#' @param prop_type proportion of bootstraps either "above_estimate", "above_low_ci", or "within_ci"
#' @param lambda_bounds bounds of rate estimation in SNPs/year/site if given time and site data
#' @param k_bounds bounds of related proportion estimation
#' @param intercept_bounds bounds of intercept estimation
#'
#' @importFrom tibble tibble
#' @importFrom dplyr mutate group_by ungroup bind_rows
#' @importFrom tidyr %>%
#' @importFrom ggplot2 ggplot aes scale_color_manual geom_hline geom_errorbar geom_point theme_bw theme
#'
#' @return list of overall outcomes, results from each time randomisation, raw results, and a plot comparing point estimates and confidence levels between original data and time randomised data
#' @export
#'
mxsure_timerandtest <- function(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist=NA, mixed_sites=NA, truncation_point=NA,
                                  sample_size=length(mixed_snp_dist), sample_n=100, permutations=5, quiet=FALSE, confidence_level=0.95,
                                start_params=NA, ci_data=NA, title=NULL, lambda_bounds = c(0, 1), k_bounds=c(0,1), intercept_bounds=c(-Inf, Inf),
                                  within_individual=FALSE, subjectA_id, subjectB_id,
                                  clustered=FALSE, prop_type="above_estimate"
                                ){

  subject_id<- time_dist<- method<- point_est<- overlapping_est<- lambda<- "5%" <- "95%" <- NULL
  #unadjusted result
  original_data <- tibble(snp_dist=mixed_snp_dist, time_dist=mixed_time_dist, sites=mixed_sites)
  original_result <- mxsure_estimate(original_data$snp_dist,unrelated_snp_dist, original_data$time_dist,original_data$sites, truncation_point=truncation_point, lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds)

  if (sample_n==0){
    return(list(
      result=tibble(NA),
      raw_results=tibble(NA),
      outcome=tibble(
        n_permutations=NA,
        any_overlapping_est=NA,
        n_overlapping_est=NA,
        perc_overlapping_est=NA,
        any_overlapping_lowci=NA,
        n_overlapping_lowci=NA,
        perc_overlapping_lowci=NA,
        prop=NA
      ),
      plot=tibble(NA)
    ))
  }


  if(anyNA(ci_data)){
    original_ci <- mxsure_ci(original_data$snp_dist,unrelated_snp_dist, original_data$time_dist,original_data$sites, truncation_point=truncation_point, quiet=quiet,
                                       sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level,
                                       start_params = c(original_result[3], original_result[2], original_result[4]), lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds)
  } else{
    original_ci <- ci_data
  }

  result <- tibble(
    "method"="Original",
    "5%"=original_ci$confidence_intervals$lambda[1],
    "point_est"=original_result$lambda,
    "95%"=original_ci$confidence_intervals$lambda[2],
    overlapping_est=FALSE,
    overlapping_low_ci=FALSE
  )

  # time rand loop
  rawtimerand <- original_ci$raw_results |>
    mutate(method = "Original")

  for(i in 1:permutations){
    if(!quiet){print(paste0("Processing Permutation: ", i))}
    if(within_individual){
      if(all(subjectA_id==subjectB_id)){
        timerand_data <- tibble(subject_id= subjectA_id, snp_dist=mixed_snp_dist, time_dist=mixed_time_dist, sites=mixed_sites)
        timerand_data <- timerand_data|>
          group_by(subject_id)|>
          mutate(time_dist=sample(time_dist, length(time_dist)))|>
          ungroup()
      } else{ warning("subject ID's do not match")}
    }else if (clustered){
      distinct_time_dist <- tibble(time_dist=mixed_time_dist)
      distinct_time_dist$permuted <- sample(distinct_time_dist$time_dist)
      timerand_data <- tibble(snp_dist=mixed_snp_dist, time_dist=mixed_time_dist, sites=mixed_sites)
      timerand_data$time_dist <- distinct_time_dist$permuted[match(timerand_data$time_dist, distinct_time_dist$time_dist)]
    }else {
      timerand_data <- tibble(snp_dist=mixed_snp_dist, time_dist=sample(mixed_time_dist, length(mixed_time_dist)), sites=mixed_sites)
    }

    if(!anyNA(start_params)){
      if (any(start_params=="Efficient")){
        start_params_timerand <- NA
    } else {
      start_params_timerand <- start_params
    }
     } else {
        start_params_timerand <- start_params
      }

  timerand_result <- mxsure_estimate(timerand_data$snp_dist,unrelated_snp_dist, timerand_data$time_dist, timerand_data$sites, truncation_point=truncation_point,
                                        lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds, start_params = start_params_timerand)

  if(!anyNA(start_params)){
  if (any(start_params=="Efficient")){
    start_params_timerand <-as.numeric(c(timerand_result[3], timerand_result[2], timerand_result[4]))
  }}

  timerand_ci <- mxsure_ci(timerand_data$snp_dist, unrelated_snp_dist, timerand_data$time_dist, timerand_data$sites,
                                       sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level, truncation_point=truncation_point,
                                       lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds
                                       ,start_params = start_params_timerand
                                       )
  if(prop_type=="above_estimate"){
  p_value_n <- sum(timerand_ci$raw_results$lambda>=original_result$lambda)
  }else if(prop_type=="within_ci"){
    p_value_n <- sum(timerand_ci$raw_results$lambda <= result$`95%`[result$method=="Original"] & timerand_ci$raw_results$lambda >= result$`5%`[result$method=="Original"])
  }else if (prop_type=="above_low_ci"){
    p_value_n <- sum(timerand_ci$raw_results$lambda >= result$`5%`[result$method=="Original"])
  }
    p_value_t <- length(timerand_ci$raw_results$lambda)

  # Append results to 'result'
  result <- bind_rows(result, tibble(
    method = paste0("TR ", i),
    `5%` = timerand_ci$confidence_intervals$lambda[1],
    point_est = timerand_result$lambda,
    `95%` = timerand_ci$confidence_intervals$lambda[2],
    p_value_n= p_value_n,
    p_value_t=p_value_t,
    overlapping_est=timerand_ci$confidence_intervals$lambda[2]>original_result$lambda,
    overlapping_low_ci=timerand_ci$confidence_intervals$lambda[2]>original_ci$confidence_intervals$lambda[1]
  ))

  # Append raw results with a method column
  rawtimerand <- bind_rows(rawtimerand,
                           timerand_ci$raw_results %>% mutate(method = paste0("TR ", i)))
  }

  result$method <- factor(result$method, levels = result$method)
  rawtimerand$method <- factor(rawtimerand$method, levels = result$method)


  outcome <- tibble(
    n_permutations=permutations,
    any_overlapping_est=sum(result$`95%`[result$method!="Original"]>=result$point_est[1])>0,
    n_overlapping_est=sum(result$`95%`[result$method!="Original"]>=result$point_est[1]),
    perc_overlapping_est=sum(result$`95%`[result$method!="Original"]>=result$point_est[1])/length(result$`95%`[result$method!="Original"]),
    any_overlapping_lowci=sum(result$`95%`[result$method!="Original"]>=result$`5%`[1])>0,
    n_overlapping_lowci=sum(result$`95%`[result$method!="Original"]>=result$`5%`[1]),
    perc_overlapping_lowci=sum(result$`95%`[result$method!="Original"]>=result$`5%`[1])/length(result$`95%`[result$method!="Original"]),
    prop=(sum(result$p_value_n, na.rm=TRUE)+1)/(sum(result$p_value_t, na.rm=TRUE)+1)
  )

    plot <- ggplot(result, aes(x = method, y = point_est, colour = as.factor(overlapping_est))) +
      scale_color_manual(values=c("FALSE"="black", "TRUE"="red3"))+
      geom_hline(yintercept = result$point_est[1], color="grey60")+
      ggbeeswarm::geom_quasirandom(data=rawtimerand,aes(x=method, y=lambda),color="grey50",size=1.1, alpha=0.4, stroke = 0)+
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.3) +
      geom_point(size = 2, color = "red3") +  # Point estimate
      scale_y_continuous(expand = c(0,0), transform = scales::pseudo_log_trans(sigma=1e-10,base=10), breaks=c(0, 10^(-10:-1)))+
      labs(title=title,
           x = NULL,
           y = paste0("Rate (",original_result$lambda_units, ")" ),
           subtitle = paste0("prop.=",format(round(outcome$prop, 4)))) +
      theme_bw()+
      theme(legend.position="none")

  xres_timerand <- list(
    result=result,
    raw_results=rawtimerand,
    outcome=outcome,
    plot=plot
  )
}
