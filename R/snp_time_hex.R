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
#'
#' @return a plot of SNP distance over time using ggplot
#'
#' @export
snp_time_hex <- function(SNPs, Time, lambda, sites=NA, snp_threshold, title="SNP-Time Hexagonal Heatmap", jitter=TRUE){
  data <- data.frame(SNPs, Time)
  if(jitter==TRUE){
    data$SNPs <- abs(jitter(SNPs))
    data$Time <- abs(jitter(Time))
  }
  if(!is.na(mean(sites, na.rm=TRUE))){lambda <- (lambda*mean(sites))/1000000}
  labely <- 0.8*max(data$SNPs[data$SNPs<=snp_threshold])
  labelx <- 0.7*max(data$Time[data$SNPs<=snp_threshold])
  ggplot(data, aes(x=Time, y=SNPs))+
    scale_y_continuous(limits = c(0, snp_threshold+1), expand = c(0.01,0.01))+
    scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
    geom_hex(#binwidth=c(1,1)
               )+
    scale_fill_continuous(#limits=c(0,25)
                          )+
    geom_abline(intercept=0, slope = lambda)+
    geom_abline(intercept=snp_threshold, slope=0, linetype="dashed", alpha=0.75)+
    geom_label(x=labelx, y=labely,
               label = paste0(round(lambda*365.25,digits = 4), " SNPs per year"), size=5)+
    labs(title=title)+
    theme_minimal()
}
