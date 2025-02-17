#' Time Randomisation Test for mixture distribution mutation rate estimation
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_n number of bootstrap sampling to conduct
#' @param start_params initial parameters for optim, if NA (as default) will try 3 different start parameters and produce the highest likelihood result. Specifying the start parameters minimises computing time.
#' @param ci_data optional input for previously calculated CI data (mixture_snp_cutoff_ci) for computational efficiency
#' @param confidence_level confidence level for CIs
#' @param title title for ggplot
#' @param unrelated_time_dist time differences for unrelated sample
#' @param unrelated_sites sites considered for unrelated sample
#' @param percent_permute proportion of mixed sample to permute
#' @param permutations number of permutation to run
#'
#' @return ggplot comparing point estimates and confidence levels between normal data and time randomised data
#' @export
#'
#' @examples
mixture_permutation_ci <- function(trans_snp_dist, trans_time_dist=NA, trans_sites=NA,
                                      unrelated_snp_dist, unrelated_time_dist=NA, unrelated_sites=NA,
                                      percent_permute=1, permutations=3, sample_n=500, confidence_level=0.95, start_params=NA, ci_data=NA, title=NULL){
  normal_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
  unrel_data <- tibble(snp_dist=unrelated_snp_dist, time_dist=unrelated_time_dist, sites=unrelated_sites)
  unrel_data <- unrel_data[complete.cases(unrel_data$time_dist), ]

  #unadjusted result
  normal_result <- mixture_snp_cutoff(normal_data$snp_dist, unrel_data$snp_dist, normal_data$time_dist,normal_data$sites)
  if(anyNA(ci_data)){
    normal_ci <- mixture_snp_cutoff_ci(normal_data$snp_dist,unrel_data$snp_dist, normal_data$time_dist,normal_data$sites,
                                        sample_n=sample_n, confidence_level=confidence_level,
                                       start_params = c(normal_result[3], normal_result[2]))
  } else{
    normal_ci <- ci_data
  }

  result <- tibble(
    "method"="Normal",
    "5%"=normal_ci$confidence_intervals$lambda[1],
    "point_est"=normal_result$lambda,
    "95%"=normal_ci$confidence_intervals$lambda[2]
  )

  result_k <- tibble(
    "method"="Normal",
    "5%"=normal_ci$confidence_intervals$k[1],
    "point_est"=normal_result$k,
    "95%"=normal_ci$confidence_intervals$k[2]
  )

  # mixture permutation loop
  rawperm <- normal_ci$raw_results |>
    mutate(method = "Normal")
  for(i in 1:permutations){
    num_swap_rel <- round(nrow(normal_data) * percent_permute)
    num_swap_unrel <- round(nrow(unrel_data) * percent_permute)

    rel_shuffle <- slice_sample(normal_data, n=nrow(normal_data))
    unrel_shuffle <- slice_sample(unrel_data, n=nrow(unrel_data))

    # comb <- rbind(rel_shuffle[1:num_swap_rel,], unrel_shuffle[1:num_swap_unrel,])
    # comb <- slice_sample(comb, n = nrow(comb))
    #
    # close_perm <- rbind(comb[1:num_swap_rel, ], rel_shuffle[(num_swap_rel+1):(nrow(rel_shuffle)+1),]) #always includes an NA row to ensure if percent permute is 100% it doesnt add the last column
    # close_perm <- close_perm[complete.cases(close_perm$snp_dist),] #gets rid of NA row
    #
    # distant_perm <- rbind(comb[(num_swap_rel+1):(nrow(comb)), ], unrel_shuffle[(num_swap_unrel+1):(nrow(unrel_shuffle)+1),])
    # distant_perm <- distant_perm[complete.cases(distant_perm$snp_dist),]

    close_perm <- rbind(unrel_shuffle[1:num_swap_rel, ], rel_shuffle[(num_swap_rel+1):(nrow(rel_shuffle)+1),])
    close_perm <- close_perm[complete.cases(close_perm$snp_dist),] #gets rid of NA row

    distant_perm <- rbind(rel_shuffle[1:num_swap_rel, ], unrel_shuffle[(num_swap_rel+1):(nrow(unrel_shuffle)+1),])
    distant_perm <- distant_perm[complete.cases(distant_perm$snp_dist),]


  perm_result <- mixture_snp_cutoff(close_perm$snp_dist, distant_perm$snp_dist, close_perm$time_dist,close_perm$sites)
  perm_ci <- mixture_snp_cutoff_ci(close_perm$snp_dist, distant_perm$snp_dist, close_perm$time_dist,close_perm$sites,
                                       sample_n=sample_n, confidence_level=confidence_level,
                                       start_params = c(perm_result[3], perm_result[2]))
  # Append results to 'result'
  result <- bind_rows(result, tibble(
    method = paste0("Permutation ", i),
    `5%` = perm_ci$confidence_intervals$lambda[1],
    point_est = perm_result$lambda,
    `95%` = perm_ci$confidence_intervals$lambda[2]
  ))

  result_k <- bind_rows(result_k, tibble(
    method = paste0("Permutation ", i),
    `5%` = perm_ci$confidence_intervals$k[1],
    point_est = perm_result$k,
    `95%` = perm_ci$confidence_intervals$k[2]
  ))

  # Append raw results with a method column
  rawperm <- bind_rows(rawperm,
                           perm_ci$raw_results %>% mutate(method = paste0("Permutation ", i)))
  }

    plot <- ggplot(result, aes(x = method, y = point_est)) +
      geom_hline(yintercept = result$point_est[1], color="grey60")+
      geom_point(data=rawperm,aes(x=method, y=lambda),color="grey50",size=0.8, alpha=0.3)+
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
      geom_point(size = 2, color = "red3") +  # Point estimate
      scale_y_continuous(transform = "log10"
                         )+
      labs(title=title,
           x = "Method",
           y = "Rate") +
      theme_minimal()

    plot2 <- ggplot(result_k, aes(x = method, y = point_est)) +
      geom_hline(yintercept = result_k$point_est[1], color="grey60")+
      geom_point(data=rawperm,aes(x=method, y=k),color="grey50",size=0.8, alpha=0.3)+
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
      geom_point(size = 2, color = "red3") +  # Point estimate
      scale_y_continuous(#transform = "log10"
      )+
      labs(title=title,
           x = "Method",
           y = "k") +
      theme_minimal()


  xres_perm <- list(
    result_lambda=result,
    result_k=result_k,
    rawperm=rawperm,
    plot_lambda=plot,
    plot_k=plot2
  )
}
