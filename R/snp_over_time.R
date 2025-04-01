library(ggplot2)

#' SNP dist/time plot
#'
#' Creates a plot of snp distance over time
#'
#' @param SNPs list of SNP distances from a mixed transmission data set
#' @param Time list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param lambda mutation rate to plot in the graph (calculated with mixture_snp_cutoff)
#' @param sites list of sites considered for each SNP distance in mixed data set; adjusts lambda to be in SNPs/day if provided
#' @param snp_threshold a threshold to apply to the mixed data set for considering related data (calculated with mixture_snp_cutoff)
#' @param title title for graph
#' @param jitter jitter data (keeps data above 0 SNPs and 0 time)
#' @param p_value optional input to display time randomised p value from mixture_timerand_test
#' @param ci optional input to display CIs of mutation rate
#'
#'
#' @return a plot of SNP distance over time using ggplot
#'
#' @export
snp_over_time <- function(SNPs, Time, lambda, sites=NA, snp_threshold, title="SNPs over Time", jitter=TRUE, p_value=NA, ci=NA, time_limits=c(0,NA)){

  data <- data.frame(SNPs, Time)
  if(jitter==TRUE){
  data$SNPs <- abs(jitter(SNPs))
  data$Time <- abs(jitter(Time))
  }

  if(!is.na(mean(sites, na.rm=TRUE))){
    lambda <- (lambda*mean(sites))/1000000}
    predictive_intervals <- tibble(Time=0:max(c(max(Time), time_limits[2]), na.rm=TRUE))
    predictive_intervals <- predictive_intervals|>
      mutate(estimate=qpois(0.5, Time*lambda))


  if(!is.na(mean(sites, na.rm=TRUE)) & !anyNA(ci)){
    ci[1] <- (ci[1]*mean(sites))/1000000
    ci[2] <- (ci[2]*mean(sites))/1000000}

    predictive_intervals <- predictive_intervals|>
      mutate(low_ci=qpois(0.025, Time*lambda),
             high_ci=qpois(0.975, Time*lambda))





  ggplot(filter(data,SNPs<=snp_threshold), aes(x=Time, y=SNPs))+
    scale_y_continuous(limits = c(0, snp_threshold+1), expand = c(0.01,0.01))+
    scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
    geom_point()+
    # geom_abline(intercept=0, slope = lambda)
    # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[1], NA), linetype="dotted")+
    # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[2], NA), linetype="dotted")+
    geom_step(data = predictive_intervals, aes(x=Time, y=estimate))+
    geom_step(data = predictive_intervals, aes(x=Time, y=low_ci), linetype="dotted")+
    geom_step(data = predictive_intervals, aes(x=Time, y=high_ci), linetype="dotted")+
    geom_abline(intercept=snp_threshold, slope=0, linetype="dashed", alpha=0.75)+
    geom_label(x=Inf, y=Inf, vjust=1, hjust=1,
               label = paste0(signif(lambda*365.25,digits = 3), " SNPs per year",
                              (ifelse(anyNA(ci),"", paste0("\n95% CI: ", signif(ci[1]*365.25, 3), ", ", signif(ci[2]*365.25, 3)))),
                              (ifelse(anyNA(p_value),"", paste0("\np-value=", format(round(p_value, 4)))))
                              ),
               size=5)+
    labs(title=title)+
    theme_minimal()
    }
