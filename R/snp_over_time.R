#' SNP dist/time plot
#'
#' Creates a plot of snp distance over time
#'
#' @param title title for graph passed to ggplot
#' @param jitter whether to jitter data (keeps data above 0 SNPs and 0 time)
#' @param mixed_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param mixed_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param mixed_sites list of sites considered for each SNP distance in mixed data set
#' @param truncation_point a SNP distance limit for the data, if set to NA will estimate as if there is no limit
#' @param max_time the time (in days) utilised to calculate SNP thresholds, only applicable when time differences are provided, if not provided will utilise maximum of supplied times
#' @param ci_data optional input for previously calculated CI data (mxsure_ci) to display
#' @param time_limits x axis limits passed to ggplot
#' @param under_threshold whether to only plot points below the calculated SNP threshold
#' @param p_value optional input to display time randomised p value from mxsure_timerand_test
#'
#' @importFrom ggplot2 ggplot aes scale_y_continuous scale_x_continuous scale_color_manual geom_point geom_step geom_hline labs theme_minimal guides guide_legend
#' @importFrom dplyr distinct mutate recode
#' @importFrom tidyr pivot_longer
#' @importFrom tibble tibble
#' @importFrom stats dpois dnbinom setNames
#' @importFrom viridis viridis
#'
#' @return a plot of SNP distance over time using ggplot
#'
#' @export
snp_over_time <- function(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites=NA, truncation_point=2000, original_result=NA, start_params=NA,
                          max_time=NA, title="SNP Distance Over Time", jitter=TRUE, p_value=NA, ci_data=NA, time_limits=c(0,NA), under_threshold=FALSE,
                          tree=NA, sampleA=NA, sampleB=NA,branch_offset=NA,
                          lambda_bounds=c(0, 1), k_bounds=c(0,1), intercept_bounds=c(-Inf, Inf)){


  snp_dist <- time_dist <-rel_lh <- unrel_lh <-LHR <-LHR_bin <- result <- estimate <- low_ci <- high_ci <- rel_loglh <- unrel_logLH <- logLH <- NULL

  if(is.na(truncation_point)){
    truncation_point <- Inf
  }

  data <- mxsure_likelyhood(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, truncation_point =
                              truncation_point, original_result = original_result, start_params = start_params, tree =
                              tree, sampleA = sampleA, sampleB = sampleB, branch_offset = branch_offset)

  data$time_dist <- abs(data$time_dist)


if(anyNA(original_result)){
  mix_res <- suppressWarnings(mxsure_estimate(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, truncation_point = 2000, start_params=start_params,
                                              tree=tree, sampleA=sampleA, sampleB=sampleB, branch_offset=branch_offset))
}else{
  mix_res <- original_result
}
  if(is.na(mean(mixed_sites, na.rm=TRUE))){
    mixed_sites <- 1
  }

  if(under_threshold){
  data <- filter(data, snp_dist<=mix_res$snp_threshold+1)
  } else {
    data <- filter(data, snp_dist<truncation_point)
  }
  #data$LHR <- LH$LHR[match(paste(data$snp_dist, data$time_dist), paste(LH$snp_dist, LH$time_dist))]
  lhr_levels <- c("LHR < 0.01",
                "0.01 \u2264 LHR < 0.1",
                "0.1 \u2264 LHR < 1",
                "1 \u2264 LHR < 10",
                "10 \u2264 LHR < 100",
                "100 \u2264 LHR")

data <- data |>
  mutate(LHR_bin = cut(logLHR,
                       breaks = c(-Inf, -2, -1, 0, 1, 2, Inf),
                       labels = lhr_levels,
                       right = FALSE)) |>
  mutate(LHR_bin = factor(LHR_bin, levels = lhr_levels))



  if(jitter==TRUE){
  data$snp_dist <- abs(jitter(data$snp_dist))
  data$time_dist <- abs(jitter(data$time_dist))
  }
    lambda <- mix_res$lambda*mean(mixed_sites) #convert snp/year/site to snp/year/genome


    predictive_intervals <- tibble(time_dist=1:max(c(max(data$time_dist), time_limits[2]), na.rm=TRUE))
    predictive_intervals <- predictive_intervals|>
      mutate(estimate=qpois(0.5, (time_dist/365.25)*lambda+mix_res$intercept +ifelse(is.null(mix_res$single_branch_lambda), 0, 2*mix_res$single_branch_lambda) ))

  if(!anyNA(ci_data)){
    ci <- c(ci_data$confidence_intervals$lambda[1]*mean(mixed_sites),
     ci_data$confidence_intervals$lambda[2]*mean(mixed_sites))
    }

    predictive_intervals <- predictive_intervals|>
      mutate(low_ci=qpois(0.025, (time_dist/365.25)*lambda+mix_res$intercept +ifelse(is.null(mix_res$single_branch_lambda), 0, 2*mix_res$single_branch_lambda)),
             high_ci=qpois(0.975, (time_dist/365.25)*lambda+mix_res$intercept +ifelse(is.null(mix_res$single_branch_lambda), 0, 2*mix_res$single_branch_lambda)))


  if(under_threshold){
  ggplot(data, aes(x=time_dist, y=snp_dist, color=LHR_bin))+
    scale_y_continuous(limits = c(0, mix_res$snp_threshold+1), expand = c(0.01,0.01))+
    scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
      scale_color_manual(
        values = setNames(viridis::viridis(8, option = "C", direction = 1)[1:6], lhr_levels),
        drop = FALSE
      ) +
    geom_point(show.legend = TRUE)+
    # geom_abline(intercept=0, slope = lambda)
    # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[1], NA), linetype="dotted")+
    # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[2], NA), linetype="dotted")+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=estimate), color="black")+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=low_ci), color="black", linetype="dotted")+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=high_ci), color="black", linetype="dotted")+
      geom_hline(yintercept=mix_res$snp_threshold, linetype="dashed", alpha=0.75, color="red3")+
    labs(title=title,
         y="SNP Distance",
         x="Time (Days)",
         color="LHR",
         subtitle=paste0(signif(lambda, digits = 3), " SNPs per year",
                         (ifelse(anyNA(ci_data),"", paste0("; 95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                         (ifelse(anyNA(p_value),"", paste0("; p-value=", format(round(p_value, 4)))))
         ))+
    theme_minimal()
  } else {
    ggplot(data, aes(x=time_dist, y=snp_dist, color=LHR_bin))+
      scale_y_continuous(limits = c(0, NA), expand = c(0.01,0.01), transform = scales::pseudo_log_trans(sigma=1, base=10),
                         breaks = c(0, signif(exp(seq(0, log(min(c(truncation_point, max(data$snp_dist)))), length.out=10)), 1)))+
      scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
      scale_color_manual(
        values = setNames(viridis::viridis(8, option = "C", direction = 1)[1:6], lhr_levels),
        drop = FALSE
      ) +
      geom_point(show.legend = TRUE)+
      # geom_abline(intercept=0, slope = lambda)
      # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[1], NA), linetype="dotted")+
      # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[2], NA), linetype="dotted")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=estimate), color="black")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=low_ci), color="black", linetype="dotted")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=high_ci), color="black", linetype="dotted")+
      geom_hline(yintercept=mix_res$snp_threshold, linetype="dashed", alpha=0.75, color="red3")+
      guides(color = guide_legend(override.aes = list(size = 3))) +
      labs(title=title,
           y="SNP Distance",
           x="Time (Days)",
           color="LHR",
           subtitle=paste0(signif(lambda, digits = 3), " SNPs per year",
                           (ifelse(anyNA(ci_data),"", paste0("; 95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                           (ifelse(anyNA(p_value),"", paste0("; p-value=", format(round(p_value, 4)))))
           ))+
      theme_minimal()
  }
    }
