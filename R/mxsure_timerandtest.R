#' Time Randomisation Test for mixture distribution substituion rate estimation
#
#' @param mixed_snp_dist list of SNP distances from a mixed transmission data set
#'
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param mixed_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param mixed_sites list of sites considered for each SNP distance in mixed data set
#' @param bootstraps number of bootstrap sampling to conduct for each time permutation
#' @param start_params initial parameters for bootstrap estimates, if NA will try various different start parameters and produce the highest likelihood result. Specifying the start parameters minimises computing time. If "Efficient" (as default) will use original result as starting paramters.
#' @param ci_data optional input for previously calculated CI data (mxsure_ci) for computational efficiency
#' @param confidence_level confidence level for CIs
#' @param title title for ggplot
#' @param right_truncation SNP distances to truncate at
#' @param permutations number of time permutation to run
#' @param quiet if false will print progress bar for each permutation
#' @param failure_criterion proportion of bootstraps either "above_estimate", "above_low_ci", or "within_ci"
#' @param lambda_bounds bounds of rate estimation in SNPs/year/site if given time and site data
#' @param k_bounds bounds of related proportion estimation
#' @param intercept_bounds bounds of intercept estimation
#' @param tree SNP-scaled phylogenetic tree or list of trees with pairs of tips labelled with sampleA and sampleB. If this is supplied a different model for related SNP distances will be fit that takes into account branch length differences to the MRCA for any given pair of samples provided in sampleA and sampleB.
#' @param sampleA tip labels for sampleA; must be in the correct order with respect to sampleB such that the time distance is calculated as SampleA date-SampleB date (even if this allows for negative numbers)
#' @param sampleB tip labels for sampleB;see above
#' @param branch_offset overide to branch offset to each skellam parameter for divergence correction model
#' @param tree_fulldist_param_bounds bounds of single branch fitting used in tree models
#' @param original_result if supplied the function will not re-estimate the original estimate
#'
#' @importFrom tidyr tibble
#' @importFrom dplyr mutate group_by ungroup bind_rows slice_sample
#' @importFrom tidyr %>%
#' @importFrom ggplot2 ggplot aes scale_color_manual geom_hline geom_errorbar geom_point theme_bw theme
#'
#' @return list of overall outcomes, results from each time randomisation, raw results, and a plot comparing point estimates and confidence levels between original data and time randomised data. Failure bootstraps and failure percentages indicate the number and percentage of bootstrap estimates from time randomised datasets above the original estimate (if prop type is not adjusted)
#'
#'#' @examples
#' mixed_distances <- simulate_mixsnp_data(lambda=5, k=0.8, n=100)
#' distant_distances <- simulate_mixsnp_data(lambda=5, k=0, n=1000)
#' result <- mxsure_estimate(mixed_snp_dist = mixed_distances$snp_dist,
#' unrelated_snp_dist = distant_distances$snp_dist,
#' mixed_time_dist = mixed_distances$time_dist)
#'
#' ci <- mxsure_ci(mixed_snp_dist = mixed_distances$snp_dist,
#' unrelated_snp_dist = distant_distances$snp_dist,
#' mixed_time_dist = mixed_distances$time_dist,
#' bootstraps=5)
#'
#' mxsure_timerandtest(mixed_snp_dist = mixed_distances$snp_dist,
#' unrelated_snp_dist = distant_distances$snp_dist,
#' mixed_time_dist = mixed_distances$time_dist,
#' permutations=1,
#' bootstraps=5,
#' original_result= result,
#' ci_data = ci) #practical useage would require more bootstraps and permutations
#'
#' @export
mxsure_timerandtest <- function(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist=NA, mixed_sites=NA, right_truncation=2000,
                                bootstraps=100, permutations=5, quiet=FALSE, confidence_level=0.95,
                                tree=NA, sampleA=NA, sampleB=NA,
                                start_params='Efficient', original_result=NA, ci_data=NA, title=NULL,
                                lambda_bounds = c(0, 1), k_bounds=c(0,1), intercept_bounds=c(-Inf, Inf),tree_fulldist_param_bounds = c(0, Inf), branch_offset=NA,
                                failure_criterion="above_estimate"
                                ){

  subject_id<- time_dist<- method<- point_est<- overlapping_est<- lambda<- "5%" <- "95%" <- NULL

  if (bootstraps==0){
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


  #unadjusted result
  original_data <- tibble(snp_dist=mixed_snp_dist, time_dist=mixed_time_dist, sites=mixed_sites)

  if(anyNA(original_result)){
    cat("Fitting base result")
  original_result <- mxsure_estimate(original_data$snp_dist,unrelated_snp_dist, original_data$time_dist,original_data$sites, right_truncation=right_truncation,
                                     tree=tree, sampleA=sampleA, sampleB=sampleB, branch_offset=branch_offset, start_params=ifelse(all(start_params=="Efficient"), NA, start_params),
                                     lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds, tree_fulldist_param_bounds = tree_fulldist_param_bounds)
  }



  if(anyNA(ci_data)){
    original_ci <- mxsure_ci(original_data$snp_dist,unrelated_snp_dist, original_data$time_dist,original_data$sites, right_truncation=right_truncation, quiet=quiet,
                             tree=tree, sampleA=sampleA, sampleB=sampleB,branch_offset=branch_offset,
                                        bootstraps=bootstraps, confidence_level=confidence_level,
                                       start_params = c(original_result[3], original_result[2], original_result[4], original_result[7], original_result[8]),
                             lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds)
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

      timerand_data <- tibble(snp_dist=mixed_snp_dist, time_dist=sample(mixed_time_dist, length(mixed_time_dist)), sites=mixed_sites)

    if(!anyNA(start_params)){
      if (any(start_params=="Efficient")){
        start_params_timerand <- NA
    } else {
      start_params_timerand <- start_params
    }
     } else {
        start_params_timerand <- start_params
      }

  timerand_result <- mxsure_estimate(timerand_data$snp_dist,unrelated_snp_dist, timerand_data$time_dist, timerand_data$sites, right_truncation=right_truncation,
                                     tree=tree, sampleA=sampleA, sampleB=sampleB,branch_offset=branch_offset,
                                        lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds, start_params = start_params_timerand,  tree_fulldist_param_bounds = tree_fulldist_param_bounds)

  if(!anyNA(start_params)){
  if (any(start_params=="Efficient")){
    start_params_timerand <-as.numeric(c(timerand_result[3], timerand_result[2], timerand_result[4], timerand_result[7], timerand_result[8]))
  }}

  timerand_ci <- mxsure_ci(timerand_data$snp_dist, unrelated_snp_dist, timerand_data$time_dist, timerand_data$sites,
                                        bootstraps=bootstraps, confidence_level=confidence_level, right_truncation=right_truncation,
                           tree=tree, sampleA=sampleA, sampleB=sampleB,branch_offset=branch_offset,
                                       lambda_bounds = lambda_bounds, k_bounds=k_bounds, intercept_bounds=intercept_bounds
                                       ,start_params = start_params_timerand
                                       )
  if(failure_criterion=="above_estimate"){
  p_value_n <- sum(timerand_ci$raw_results$lambda>=original_result$lambda)
  }else if(failure_criterion=="within_ci"){
    p_value_n <- sum(timerand_ci$raw_results$lambda <= result$`95%`[result$method=="Original"] & timerand_ci$raw_results$lambda >= result$`5%`[result$method=="Original"])
  }else if (failure_criterion=="above_low_ci"){
    p_value_n <- sum(timerand_ci$raw_results$lambda >= result$`5%`[result$method=="Original"])
  }
    p_value_t <- length(timerand_ci$raw_results$lambda)

  # Append results to 'result'
  result <- bind_rows(result, tibble(
    method = paste0("TR ", i),
    `5%` = timerand_ci$confidence_intervals$lambda[1],
    point_est = timerand_result$lambda,
    `95%` = timerand_ci$confidence_intervals$lambda[2],
    failed_bootstraps= p_value_n,
    total_bootstraps=p_value_t,
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
    failure_perc=(sum(result$failed_bootstraps, na.rm=TRUE))/(sum(result$total_bootstraps, na.rm=TRUE))
  )

  if(original_result$lambda_units=="SNPs per year per site"){
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
           subtitle = paste0("fail %=",format(round(outcome$failure_perc, 4)))) +
      theme_bw()+
      theme(legend.position="none")
  }else{
    plot <- ggplot(result, aes(x = method, y = point_est, colour = as.factor(overlapping_est))) +
      scale_color_manual(values=c("FALSE"="black", "TRUE"="red3"))+
      geom_hline(yintercept = result$point_est[1], color="grey60")+
      ggbeeswarm::geom_quasirandom(data=rawtimerand,aes(x=method, y=lambda),color="grey50",size=1.1, alpha=0.4, stroke = 0)+
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.3) +
      geom_point(size = 2, color = "red3") +  # Point estimate
      #scale_y_continuous(expand = c(0,0), transform = scales::pseudo_log_trans(sigma=1e-10,base=10), breaks=c(0, 10^(-10:-1)))+
      labs(title=title,
           x = NULL,
           y = paste0("Rate (",original_result$lambda_units, ")" ),
           subtitle = paste0("fail %=",format(round(outcome$failure_perc, 4)))) +
      theme_bw()+
      theme(legend.position="none")
    }

  xres_timerand <- list(
    result=result,
    raw_results=rawtimerand,
    outcome=outcome,
    plot=plot
  )
}
