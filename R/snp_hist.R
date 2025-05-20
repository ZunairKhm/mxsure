#' SNP distace histogram
#'
#' Creates a histogram of snp distances with or without the unrelated dataset
#'
#' @param mixed_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param snp_threshold a threshold to overlay on the graph (calculated with mixture_snp_cutoff)
#' @param limits snp distance/x axis limits
#' @param title title to pass to ggplot
#' @param scales scales to pass to facet_wrap
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous geom_histogram geom_vline facet_wrap labs theme_minimal
#'
#' @return a histogram plot of snp distances using ggplot
#'
#' @export
snp_hist <- function(mixed_snp_dist, unrelated_snp_dist=NULL, snp_threshold=NULL, limits=c(NA,100), title="SNP Distance Histogram", scales="fixed"){
    # Combine the datasets into a data frame
    data <- data.frame(
      SNP_Distance = c(unrelated_snp_dist, mixed_snp_dist),
      Dataset = c(rep("Unrelated", length(unrelated_snp_dist)),rep("Mixed", length(mixed_snp_dist)))
    )
    data$Dataset <- as.factor(data$Dataset)

    # Create the plot
    ggplot(data, aes(x = SNP_Distance)) +
      scale_x_continuous(limits = limits) +
      geom_histogram(binwidth = 1, position = "identity") +
      geom_vline(xintercept = snp_threshold + 0.5, color = "red3", alpha = 0.7, linetype = "dashed") +
      facet_wrap(~Dataset, ncol=1, scales=scales
                 ) +
      labs(
        title = title,
        x = "SNP Distances",
        y = "Count"
      ) +
      theme_minimal()
    }
