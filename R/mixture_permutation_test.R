#' Related/Unrelated permutation test for mixture distribution mutation rate estimation
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_n number of bootstrap sampling to conduct
#' @param start_params initial parameters for optim, if NA (as default) will try 3 different start parameters and produce the highest likelyhood result. Specifying the start parameters minimises computing time.
#' @param ci_data optional input for previously calculated CI data (mixture_snp_cutoff_ci) for computational efficiency
#' @param confidence_level confidence level for CIs
#' @param title title for ggplot
#' @param unrelated_time_dist list of time distances from an unrelated data set
#' @param unrelated_sites list of sites considered for each SNP distance in a distant data set
#' @param percent_permute proportion of mixed sample to permute
#'
#' @return data and plot comparing normal results to mixed/distant permuted results
#' @export
#'
#' @examples
mixture_permutation_test <- function(trans_snp_dist, trans_time_dist=NA, trans_sites=NA,
                                      unrelated_snp_dist, unrelated_time_dist=NA, unrelated_sites=NA,
                                  percent_permute=1, sample_n=500, confidence_level=0.95, start_params=NA, ci_data=NA, title=NULL){

  if ((length(trans_snp_dist) < 30) | (length(unrelated_snp_dist) < 30 | sample_n==0)){
    warning("Insufficient data points to fit distributions!")
    return(
      list(
        comparison=c(NA),
        p_value=c(NA),
        normal_point_est=c(NA),
        normal_ci=c(NA),
        perm_ci=c(NA),
        plot=c(NA)
      )
    )
  }
  normal_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
  unrel_data <- tibble(snp_dist=unrelated_snp_dist, time_dist=unrelated_time_dist, sites=unrelated_sites)
  unrel_data <- unrel_data[complete.cases(unrel_data$time_dist), ]
  normal_point_est <-suppressWarnings(mixture_snp_cutoff(normal_data$snp_dist,unrel_data$snp_dist, normal_data$time_dist,normal_data$sites))

  if(anyNA(ci_data)){
    normal_ci <- mixture_snp_cutoff_ci(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites,
                                      sample_n=sample_n, confidence_level=confidence_level,
                                       start_params = c(normal_point_est[3], normal_point_est[2]))
  } else{
    normal_ci <- ci_data
  }

  perm_point_ests <- furrr::future_map_dfr(1:sample_n, ~{

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

    z <- plyr::try_default( #suppresses warnings and errors from mixture_snp_cutoffs
      suppressWarnings(
        mixture_snp_cutoff(
          close_perm$snp_dist, distant_perm$snp_dist, close_perm$time_dist, close_perm$sites, start_params = c(normal_point_est[3], normal_point_est[2])
        ), classes = "warning"),
      data.frame(snp_threshold=NA,k=NA,k=NA,estimated_fp=NA))
    z[1:4]
  },.progress=TRUE, .options = furrr::furrr_options(seed = TRUE))

  lowerres <- perm_point_ests|> #finds quantiles for confidence intervals
    summarise(across(everything(),  ~quantile(.x, 1-confidence_level, na.rm=TRUE)))
  upperres <- perm_point_ests|>
    summarise(across(everything(), ~quantile(.x, confidence_level, na.rm=TRUE)))

  perm_ci <- bind_rows(lowerres, upperres)
  perm_ci <- list(confidence_intervals=perm_ci, raw_results=perm_point_ests)


  comparison <- tibble(
    method = factor(c("Normal", "Mixture Permuted"),levels = c("Normal", "Mixture Permuted")),
    `5%` = c(normal_ci$confidence_intervals$lambda[1], perm_ci$confidence_intervals$lambda[1]),
    point_est = c(normal_point_est$lambda, NA),
    `95%` = c(normal_ci$confidence_intervals$lambda[2], perm_ci$confidence_intervals$lambda[2])
  )

  p_value <- (sum(perm_point_ests$lambda>=normal_point_est$lambda)+1)/
    (length(perm_point_ests$lambda)+1)

  final_result <- list(
    comparison=comparison,
    p_value=p_value,
    normal_point_est=normal_point_est,
    normal_ci=normal_ci,
    perm_ci=perm_ci
  )

  plot <-  ggplot(final_result$comparison, aes(x = method, y = point_est)) +
    geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
    geom_point(data=tibble("point_est"=final_result$normal_ci$raw_results$lambda, "method"="Normal"),color="grey50",size=1, alpha=0.3)+
    geom_point(data=tibble("point_est"=final_result$perm_ci$raw_results$lambda, "method"="Mixture Permuted"),color="grey50",size=1, alpha=0.3)+
    geom_point(y=final_result$normal_point_est$lambda, x="Normal", size = 2, color = "red3") +  # Point estimate
    annotate("label", label=paste0("p=",format(round(final_result$p_value, 4), nsmall = 4)), x=Inf, y=Inf, vjust=1, hjust=1)+
    labs(title=title,
         x = "Method",
         y = "Rate") +
    theme_minimal()


    comparison2 <- tibble(
      method = factor(c("Normal", "Mixture Permuted"),levels = c("Normal", "Mixture Permuted")),
      `5%` = c(normal_ci$confidence_intervals$k[1], perm_ci$confidence_intervals$k[1]),
      point_est = c(normal_point_est$k, NA),
      `95%` = c(normal_ci$confidence_intervals$k[2], perm_ci$confidence_intervals$k[2])
    )

    p_value2 <- (sum(perm_point_ests$k>=normal_point_est$k)+1)/
      (length(perm_point_ests$k)+1)

    plot2 <-  ggplot(comparison2, aes(x = method, y = point_est)) +
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
      geom_point(data=tibble("point_est"=final_result$normal_ci$raw_results$k, "method"="Normal"),color="grey50",size=1, alpha=0.3)+
      geom_point(data=tibble("point_est"=final_result$perm_ci$raw_results$k, "method"="Mixture Permuted"),color="grey50",size=1, alpha=0.3)+
      geom_point(y=final_result$normal_point_est$k, x="Normal", size = 2, color = "red3") +  # Point estimate
      annotate("label", label=paste0("p=",format(round(p_value2, 4), nsmall = 4)), x=Inf, y=Inf, vjust=1, hjust=1)+
      labs(title=title,
           x = "Method",
           y = "k") +
      theme_minimal()

    final_result <- list(
      comparison_lambda=comparison,
      comparison_k=comparison2,
      p_value_lambda=p_value,
      p_value_k=p_value2,
      normal_point_est=normal_point_est,
      normal_ci=normal_ci,
      perm_ci=perm_ci,
      plot_lambda=plot,
      plot_k=plot2
    )


  return(final_result)




}
