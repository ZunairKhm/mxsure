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
#'
#' @return a plot of SNP distance over time using ggplot
#'
#' @export
snp_over_time <- function(SNPs, Time, lambda, sites=NA, snp_threshold){
  data <- data.frame(SNPs, Time)
  data$SNPs <- abs(jitter(SNPs))
  data$Time <- abs(jitter(Time))
  if(!is.na(mean(sites, na.rm=TRUE))){lambda <- (lambda*mean(sites))/1000000}
  labely <- 0.8*max(data$SNPs[data$SNPs<=10])
  labelx <- 0.9*max(data$Time[data$SNPs<=10])
  ggplot(filter(data,SNPs<=10), aes(x=Time, y=SNPs))+
    scale_y_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
    scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
    geom_point()+
    geom_abline(intercept=0, slope = lambda)+
    geom_abline(intercept=snp_threshold, slope=0, linetype="dashed", alpha=0.75)+
    geom_label(x=labelx, y=labely,
               label = paste0(round(lambda*365.25,digits = 4), " SNPs per year"), size=5)
}
