library(ggplot2)

#' SNP distace histogram
#'
#' Creates a histogram of snp distances with or without the unrelated dataset
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param snp_threshold a threshold to overlay on the graph (calculated with mixture_snp_cutoff)
#' @param limits snp distance/x axis limits
#'
#' @return a histogram plot of snp distances using ggplot
#'
#' @export
snp_hist <- function(trans_snp_dist, unrelated_snp_dist=NULL, snp_threshold=NULL, limits=c(NA,100))
  ggplot(data=data.frame(x=append(trans_snp_dist, unrelated_snp_dist)), aes(x))+
  scale_x_continuous(limits = limits)+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=snp_threshold+0.5, color="red")
