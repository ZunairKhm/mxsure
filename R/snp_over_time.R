#' SNP distance/time plot
#'
#' Creates a plot of snp distance over time
#'
#' @param title title for graph passed to ggplot
#' @param jitter whether to jitter data (keeps data above 0 SNPs and 0 time)
#' @param mixed_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param mixed_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param mixed_sites list of sites considered for each SNP distance in mixed data set
#' @param right_truncation a SNP distance limit for the data, if set to NA will estimate as if there is no limit
#' @param threshold_time the time (in days) utilised to calculate SNP thresholds, only applicable when time differences are provided, if not provided will utilise maximum of supplied times
#' @param ci_data optional input for previously calculated CI data (mxsure_ci) to display on subtitle
#' @param time_limits x axis limits passed to ggplot
#' @param under_threshold whether to only plot points below the calculated SNP threshold
#' @param failure_perc optional input to display time randomised failure percentage from mxsure_timerand_test
#' @param original_result optional input to provide original mxsure_likelyhood() result instead of computing it again
#' @param start_params initial parameters for lambda, k, and intercept parameters for optimisation. If NA (as default) will try a range of different start parameters and produce the highest likelihood result
#' @param tree SNP-scaled phylogenetic tree or list of trees with pairs of tips labelled with sampleA and sampleB. If this is supplied a different model for related SNP distances will be fit that takes into account branch length differences to the MRCA for any given pair of samples provided in sampleA and sampleB.
#' @param sampleA tip labels for sampleA; must be in the correct order with respect to sampleB such that the time distance is calculated as SampleA date-SampleB date (even if this allows for negative numbers)
#' @param sampleB tip labels for sampleB; see above
#' @param branch_offset overide to branch offset to each skellam parameter for divergence correction model
#' @param lambda_bounds bounds of rate estimation in SNPs/year/site if given time and site data
#' @param k_bounds bounds of related proportion estimation
#' @param intercept_bounds bounds of intercept estimation
#'
#' @importFrom ggplot2 ggplot aes scale_y_continuous scale_x_continuous scale_color_manual geom_point geom_step geom_hline labs theme_minimal guides guide_legend scale_linetype_manual
#' @importFrom dplyr distinct mutate recode rowwise
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr tibble
#' @importFrom stats dpois dnbinom setNames
#' @importFrom viridis viridis
#'
#' @return a plot of SNP distance over time using ggplot
#'
#' @examples
#' mixed_distances <- simulate_mixsnp_data(lambda=5, k=0.8, n=100)
#' distant_distances <- simulate_mixsnp_data(lambda=5, k=0, n=1000)
#' result <- mxsure_estimate(mixed_snp_dist = mixed_distances$snp_dist,
#' unrelated_snp_dist = distant_distances$snp_dist,
#' mixed_time_dist = mixed_distances$time_dist)
#'
#' snp_over_time(mixed_snp_dist = mixed_distances$snp_dist,
#' unrelated_snp_dist = distant_distances$snp_dist,
#' mixed_time_dist = mixed_distances$time_dist,
#' original_result = result)
#'
#' @export
snp_over_time <- function(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites=NA, right_truncation=2000, original_result=NA, start_params=NA,
                          threshold_time=NA, title="SNP Distance Over Time", jitter=TRUE, failure_perc=NA, ci_data=NA, time_limits=c(0,NA), under_threshold=FALSE,
                          tree=NA, sampleA=NA, sampleB=NA,branch_offset=NA,
                          lambda_bounds=c(0, 1), k_bounds=c(0,1), intercept_bounds=c(-Inf, Inf)){


  snp_dist <- time_dist <-rel_lh <- unrel_lh <-LHR <-LHR_bin <- result <- estimate <- low_ci <- high_ci <- rel_loglh <- unrel_logLH <- logLH <- NULL

  if(is.na(right_truncation)){
    right_truncation <- Inf
  }

  data <- mxsure_likelihood(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, right_truncation =
                              right_truncation, original_result = original_result, start_params = start_params, tree =
                              tree, sampleA = sampleA, sampleB = sampleB, branch_offset = branch_offset)

  data$time_dist <- abs(data$time_dist)


if(anyNA(original_result)){
  mix_res <- suppressWarnings(mxsure_estimate(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, right_truncation = 2000, start_params=start_params,
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
    data <- filter(data, snp_dist<right_truncation)
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
      rowwise()|>
      mutate(low_ci=ifelse((time_dist/365.25)*lambda+mix_res$intercept<0, -1, qpois(0.025, (time_dist/365.25)*lambda+mix_res$intercept )),
             high_ci=ifelse((time_dist/365.25)*lambda+mix_res$intercept<0, -1, qpois(0.975, (time_dist/365.25)*lambda+mix_res$intercept )),
             estimate= ifelse((time_dist/365.25)*lambda+mix_res$intercept<0, -1, qpois(0.5, (time_dist/365.25)*lambda+mix_res$intercept )))



  if(!anyNA(ci_data)){
    ci <- c(ci_data$confidence_intervals$lambda[1]*mean(mixed_sites),
     ci_data$confidence_intervals$lambda[2]*mean(mixed_sites))
    }



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
      geom_step(data = predictive_intervals, aes(x=time_dist, y=estimate, linetype="Expected SNPs"), color="black")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=low_ci, linetype="95% Predictive Interval"), color="black")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=high_ci, linetype="95% Predictive Interval"), color="black")+
      # SNP threshold
      geom_hline(aes(yintercept=mix_res$snp_threshold, linetype="SNP Threshold"), color="red3")+
      scale_linetype_manual(
        values = c("Expected SNPs"="solid", "95% Predictive Interval"="dotted", "SNP Threshold"="dashed"),
        breaks = c("Expected SNPs", "95% Predictive Interval", "SNP Threshold"),
        name = NULL
      ) +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 3)),
        linetype = guide_legend(
          override.aes = list(
            shape = NA,
            colour = c("black", "black", "red3"),
            size = c(0.8, 0.8, 0.8)
          )
        )
      ) +
    labs(title=title,
         y="SNP Distance",
         x="Time (Days)",
         color="LHR",
         subtitle=paste0(signif(lambda, digits = 3), " SNPs per year",
                         (ifelse(anyNA(ci_data),"", paste0("; 95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                         (ifelse(anyNA(failure_perc),"", paste0("; failure %=", format(round(failure_perc, 4)))))
         ))+
    theme_minimal()
  } else {
    ggplot(data, aes(x=time_dist, y=snp_dist, color=LHR_bin))+
      scale_y_continuous(limits = c(0, NA), expand = c(0.01,0.01), transform = scales::pseudo_log_trans(sigma=1, base=10),
                         breaks = c(0, signif(exp(seq(0, log(min(c(right_truncation, max(data$snp_dist)))), length.out=10)), 1)))+
      scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01), oob=scales::censor)+
      scale_color_manual(
        values = setNames(viridis::viridis(8, option = "C", direction = 1)[1:6], lhr_levels),
        drop = FALSE
      ) +
      geom_point(show.legend = TRUE)+
      # Predictive intervals
      geom_step(data = predictive_intervals, aes(x=time_dist, y=estimate, linetype="Expected SNPs"), color="black")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=low_ci, linetype="95% Predictive Interval"), color="black")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=high_ci, linetype="95% Predictive Interval"), color="black")+
      # SNP threshold
      geom_hline(aes(yintercept=mix_res$snp_threshold, linetype="SNP Threshold"), color="red3")+
      scale_linetype_manual(
        values = c("Expected SNPs"="solid", "95% Predictive Interval"="dotted", "SNP Threshold"="dashed"),
        breaks = c("Expected SNPs", "95% Predictive Interval", "SNP Threshold"),
        name = NULL
      ) +
      guides(
        color = guide_legend(override.aes = list(shape = 16, size = 3)),
        linetype = guide_legend(
          override.aes = list(
            shape = NA,
            colour = c("black", "black", "red3"),
            size = c(0.8, 0.8, 0.8)
          )
        )
      ) +
      labs(title=title,
           y="SNP Distance",
           x="Time (Days)",
           color="LHR",
           subtitle=paste0(signif(lambda, digits = 3), " SNPs per year",
                           (ifelse(anyNA(ci_data),"", paste0("; 95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                           (ifelse(anyNA(failure_perc),"", paste0("; failure %=", format(round(failure_perc, 4)))))
           ))+
      theme_minimal()
  }
    }
