#' Threshold Comparison Plot
#'
#' @param mixed_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param snp_thresholds list of SNP thresholds to consider
#' @param method_names list of names associated with each SNP threshold
#' @param identity whether to include an identity threshold as well
#' @param percent_dist_mixed percent distances from a mixed transmission data set, needed if identity thresholds are considered
#' @param percent_dist_distant percent distances from a distant transmission data set, needed if identity thresholds are considered
#' @param identity_thresh identity threshold to use (as percentage)
#' @param k_est proportion related line to plot
#' @param title title of plot
#' @param labels whether to include labels on the plot
#' @param dataset_names names of the mixed/distant datasets
#' @param show_related whether to show related dataset comparisons
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual scale_x_continuous geom_vline geom_text position_dodge labs theme_bw theme
#' @importFrom purrr map2_dbl
#'
#' @return plot comparing proportion of data under each threshold
#' @export
#'
threshold_comparison_plot <- function(mixed_snp_dist, unrelated_snp_dist, snp_thresholds, method_names, show_related=FALSE,
                                 identity=FALSE, percent_dist_mixed, percent_dist_distant, identity_thresh=99.99,
                                 k_est=NULL, title=NULL,labels=TRUE, dataset_names=c("Mixed", "Distant")) {
  data <- tibble(method = method_names, snp_threshold = snp_thresholds)

  data <- data |>
    mutate(
      n_mixed=map2_dbl(snp_threshold, method, ~sum(mixed_snp_dist <= .x)),
      tot_mixed=length(mixed_snp_dist),
      n_distant=map2_dbl(snp_threshold, method, ~sum(unrelated_snp_dist <= .x)),
      tot_distant = length(unrelated_snp_dist)
    )|>
    mutate(
      p_mixed=n_mixed/tot_mixed,
      p_distant=n_distant/tot_distant
    )

  if(identity==TRUE){
   id_data <- tibble(
     method=paste0(identity_thresh, "% Identity"),
     snp_threshold=NA,
     n_mixed=sum(percent_dist_mixed <= (1-(identity_thresh/100))),
     tot_mixed=length(percent_dist_mixed),
     n_distant=sum(percent_dist_distant <= (1-(identity_thresh/100))),
     tot_distant=length(percent_dist_distant)
   )
   id_data <- id_data|>
     mutate(
       p_mixed=n_mixed/tot_mixed,
       p_distant=n_distant/tot_distant
     )
   data <- rbind(data, id_data)
  }

  data <- data|>
  pivot_longer(
    cols = c(n_mixed, tot_mixed, p_mixed, n_distant, tot_distant, p_distant),
    names_to = c(".value", "dataset"),
    names_pattern = "(n|tot|p)_(mixed|distant)"
  ) |>
    mutate(
      dataset = recode(dataset, mixed = dataset_names[1], distant = dataset_names[2]),

    )|>
    mutate(dataset=factor(dataset, levels=dataset_names))|>
    mutate(label = paste0(n, "/", tot, " (", round(p*100, 2), "%)"))

 #return(data)

  if(show_related){

    ggplot(data, aes(x = p, y = method, fill = dataset)) +
      geom_bar(stat = "identity", position = "dodge", width=0.9) +
      scale_fill_manual(values = c( "darkseagreen", "red3")) +
      scale_x_continuous(limits=c(0,1), expand = c(0,0))+
      geom_vline(xintercept=k_est, linetype="solid", colour="grey20", alpha=0.6)+
      geom_text(aes(label = if (labels) label else NA, x = 0.5),
                position = position_dodge(width = 0.9),
                color = "black") +
      labs(
        title = title,
        x = "Shared strain proportion",
        y = "Threshold Method",
        fill = "Relationship"
      ) +
      theme_bw()+
      theme(legend.position = "bottom")

  }else{
    ggplot(filter(data,dataset==dataset_names[2]), aes(x = p, y = method, fill = method)) +
      geom_bar(stat = "identity", position = "dodge", width=0.9) +
      scale_x_continuous(limits=c(0,1), expand = c(0,0))+
      geom_vline(xintercept=k_est, linetype="solid", colour="grey20", alpha=0.6)+
      geom_text(aes(label = if (labels) label else NA, x = 0.5),
                position = position_dodge(width = 0.9),
                color = "black") +
      labs(
        title = title,
        x = "Estimated False Positive Rate",
        y = "Threshold Method",
        fill = NULL
      ) +
      theme_bw()+
      theme(legend.position = "none")
  }


}
