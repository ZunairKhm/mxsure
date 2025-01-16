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
snp_hist <- function(trans_snp_dist, unrelated_snp_dist=NULL, snp_threshold=NULL, limits=c(NA,100), title="SNP Distance Histogram", separate_colour=FALSE){
  if (separate_colour==TRUE){
    # Combine the datasets into a data frame
    data <- data.frame(
      SNP_Distance = c(unrelated_snp_dist, trans_snp_dist),
      Dataset = c(rep("Unrelated", length(unrelated_snp_dist)),rep("Mixed", length(trans_snp_dist)))

    )

    # Create the plot
    ggplot(data, aes(x = SNP_Distance, fill = Dataset)) +
      scale_x_continuous(limits = limits) +
      geom_histogram(data = subset(data, Dataset == "Unrelated"),
                     binwidth = 1, position = "identity", alpha = 0.6) +
      geom_histogram(data = subset(data, Dataset == "Mixed"),
                     binwidth = 1, position = "identity", alpha = 0.8) +
      scale_fill_manual(values = c("Mixed" = "darkseagreen4", "Unrelated" = "red3")) +
      geom_vline(xintercept = snp_threshold + 0.5, color = "red3",alpha=0.7, linetype = "dashed") +
      labs(
        title = title,
        x = "SNP Distances",
        y = "Count",
        fill = "Dataset"
      ) +
      theme_minimal()
    } else{
  ggplot(data=data.frame(x=append(trans_snp_dist, unrelated_snp_dist)), aes(x))+
  scale_x_continuous(limits = limits)+
  geom_histogram(binwidth = 1)+
  geom_vline(xintercept=snp_threshold+0.5, color="red3",alpha=0.7, linetype="dashed")+
  labs(title = title, x="SNP Distances", y="Count")+
  theme_minimal()
  }
}
