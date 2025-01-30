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
#' @param confidence_level confidence level for CIs
#' @param title title for ggplot
#' @param within_individual whether to permute time data within each individual
#' @param subjectA_id subject id for within individual permutation
#' @param subjectB_id subject id for within individual permutation
#'
#' @return ggplot comparing point estimates and confidence levels between normal data and time randomised data
#' @export
#'
#' @examples
mixture_timerand_test <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                                  within_individual=FALSE, subjectA_id=NA, subjectB_id=NA,
                                  sample_size=length(trans_snp_dist), sample_n=500, confidence_level=0.95,
                                  start_params=NA, ci_data=NA, title=NULL){

  if ((length(trans_snp_dist) < 30) | (length(unrelated_snp_dist) < 30)){
    warning("Insufficient data points to fit distributions!")
    return(
      list(
        comparison=c(NA),
        p_value=c(NA),
        normal_point_est=c(NA),
        normal_ci=c(NA),
        timerand_ci=c(NA),
        plot=c(NA)
      )
    )
  }
  normal_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
  normal_point_est <-suppressWarnings(mixture_snp_cutoff(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites))

  if(anyNA(ci_data)){
    normal_ci <- mixture_snp_cutoff_ci(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites,
                                      sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level,
                                      start_params = c(normal_point_est[3], normal_point_est[2]))
  } else{
    normal_ci <- ci_data
  }

  timerand_point_ests <- furrr::future_map_dfr(1:sample_n, ~{
    if(within_individual){
      if(all(subjectA_id==subjectB_id)){
        x <- tibble(subject_id= subjectA_id, snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
        x <- x|>
          group_by(subject_id)|>
          mutate(time_dist=sample(time_dist, length(time_dist)))|>
          ungroup()
      } else{ warning("subject ID's do not match")}
    }else {
      x <- tibble(snp_dist=trans_snp_dist, time_dist=sample(trans_time_dist, length(trans_time_dist)), sites=trans_sites)
    }

    y <- unrelated_snp_dist
    z <- plyr::try_default( #suppresses warnings and errors from mixture_snp_cutoffs
      suppressWarnings(
        mixture_snp_cutoff(
          x$snp_dist,y, x$time_dist, x$sites, start_params = c(normal_point_est[3], normal_point_est[2])
        ), classes = "warning"),
      data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA))
    z[1:4]
  },.progress=TRUE, .options = furrr::furrr_options(seed = TRUE))

  timerand_point_ests <- timerand_point_ests[complete.cases(timerand_point_ests), ] #removes failed results from optim failures

  lowerres <- timerand_point_ests|> #finds quantiles for confidence intervals
    summarise(across(everything(),  ~quantile(.x, 1-confidence_level, na.rm=TRUE)))
  upperres <- timerand_point_ests|>
    summarise(across(everything(), ~quantile(.x, confidence_level, na.rm=TRUE)))

  timerand_ci <- bind_rows(lowerres, upperres)
  timerand_ci <- list(confidence_intervals=timerand_ci, raw_results=timerand_point_ests)

  comparison <- tibble(
  method = c("Normal", "Time Randomised"),
  `5%` = c(normal_ci$confidence_intervals$lambda[1], timerand_ci$confidence_intervals$lambda[1]),
  point_est = c(normal_point_est$lambda, NA),
  `95%` = c(normal_ci$confidence_intervals$lambda[2], timerand_ci$confidence_intervals$lambda[2])
)

  p_value <- (sum(timerand_point_ests$lambda>=normal_point_est$lambda)+1)/
              (length(timerand_point_ests$lambda)+1)

  final_result <- list(
    comparison=comparison,
    p_value=p_value,
    normal_point_est=normal_point_est,
    normal_ci=normal_ci,
    timerand_ci=timerand_ci
  )

  plot <-  ggplot(final_result$comparison, aes(x = method, y = point_est)) +
    geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
    geom_point(data=tibble("point_est"=final_result$normal_ci$raw_results$lambda, "method"="Normal"),color="grey50",size=1, alpha=0.3)+
    geom_point(data=tibble("point_est"=final_result$timerand_ci$raw_results$lambda, "method"="Time Randomised"),color="grey50",size=1, alpha=0.3)+
    geom_point(y=final_result$normal_point_est$lambda, x="Normal", size = 2, color = "red3") +  # Point estimate
    annotate("label", label=paste0("p=",format(round(final_result$p_value, 4), nsmall = 4)), x=Inf, y=Inf, vjust=1, hjust=1)+
    labs(title=title,
         x = "Method",
         y = "Rate") +
    theme_minimal()

  final_result <- list(
    comparison=comparison,
    p_value=p_value,
    normal_point_est=normal_point_est,
    normal_ci=normal_ci,
    timerand_ci=timerand_ci,
    plot=plot
  )

  return(final_result)




}
