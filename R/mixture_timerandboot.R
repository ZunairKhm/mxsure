#' Mixture Timerandboot
#'
#' Creates confidence intervals for outputs from mixture_snp_cutoff using bootstrapping. Utilises the furrr package to allow for parallel computation, to enable use the plan() function from the future package.
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_size size of each bootstrap sample
#' @param sample_n number of bootstrap sampling to conduct
#' @param start_params initial parametrs for optim, if NA (as default) will try 3 different start parameters and produce the highest likelyhood result. Specifying the start parameters minimises computing time.
#' @param truncation_point a SNP distance limit for the data, if set to NA will estimate as if there is no limit
#' @param confidence_level confidence level to produce confidence intervals
#'
#' @importFrom furrr future_map_dfr furrr_options
#'
#' @return Confidence intervals
#'
#' @export
mixture_timerandboot <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA, ci_data=NA, title=NULL,
                                  sample_size=length(trans_snp_dist),truncation_point=NA, sample_n=1000, confidence_level=0.95, start_params=NA){

  if(is.na(truncation_point)){
    truncation_point <- Inf
  }

  #unadjusted result
  normal_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)
  normal_result <- mixture_snp_cutoff(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites, truncation_point=truncation_point)
  if(anyNA(ci_data)){
    normal_ci <- mixture_snp_cutoff_ci(normal_data$snp_dist,unrelated_snp_dist, normal_data$time_dist,normal_data$sites, truncation_point=truncation_point,
                                       sample_size=sample_size, sample_n=sample_n, confidence_level=confidence_level,
                                       start_params = c(normal_result[3], normal_result[2]))
  } else{
    normal_ci <- ci_data
  }

  plot_data <- tibble(
    "method"=factor(c("Normal"), levels=c("Normal", "Time Randomised")),
    "5%"=normal_ci$confidence_intervals$lambda[1],
    "point_est"=normal_result$lambda,
    "95%"=normal_ci$confidence_intervals$lambda[2]
  )

  if (sample_n==0){
    lowerres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)
    upperres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)

    ci <- bind_rows(lowerres, upperres)
    res <- list(confidence_intervals=ci, raw_results=NA)
  }
  else if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){

    mix_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)

    # if (anyNA(start_params)){
    #   test_result <- suppressWarnings(
    #     mixture_snp_cutoff(
    #       mix_data$snp_dist,unrelated_snp_dist, mix_data$time_dist, mix_data$sites,truncation_point=truncation_point, start_params = NA
    #     ), classes = "warning")
    #   start_params <- c(test_result[3], test_result[2])
    # }

    mix_data <- filter(mix_data, snp_dist<truncation_point)
    unrelated_snp_dist <- unrelated_snp_dist[unrelated_snp_dist<truncation_point]

    #bootstrapping both close and distant data sets allowing for parallelisationg
    bootstrapresults <- furrr::future_map_dfr(1:sample_n, ~{
      timerand_data <- tibble(snp_dist=mix_data$snp_dist, time_dist=sample(mix_data$time_dist, size=length(mix_data$time_dist)), sites=mix_data$sites)
      x <- slice_sample(timerand_data, n= sample_size, replace = TRUE)
      y <- sample(unrelated_snp_dist, size=length(unrelated_snp_dist), replace=TRUE)
      z <- plyr::try_default( #suppresses warnings and errors from mixture_snp_cutoffs
        suppressWarnings(
          mixture_snp_cutoff(
            x$snp_dist,y, x$time_dist, x$sites, truncation_point=truncation_point, start_params = start_params
          ), classes = "warning"),
        data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA))
      #return(list(timerand_data, x,y,z ))
                  z[1:4]
    },.progress=TRUE, .options = furrr::furrr_options(seed = TRUE))

    bootstrapresults <- bootstrapresults[complete.cases(bootstrapresults), ] #removes failed results from optim failures

    lowerres <- bootstrapresults|> #finds quantiles for confidence intervals
      summarise(across(everything(),  ~quantile(.x, 1-confidence_level, na.rm=TRUE)))
    upperres <- bootstrapresults|>
      summarise(across(everything(), ~quantile(.x, confidence_level, na.rm=TRUE)))

    ci <- bind_rows(lowerres, upperres)
    p_value <- sum(bootstrapresults$lambda>=normal_result$lambda)/length(bootstrapresults$lambda)


    plot_data <- bind_rows(plot_data, tibble(
      "method"=factor(c("Time Randomised"), levels=c("Normal", "Time Randomised")),
      "5%"=lowerres$lambda,
      "point_est"=NA,
      "95%"=upperres$lambda
    ))

    bootstrapresults$method <- factor(c("Time Randomised"), levels=c("Normal", "Time Randomised"))

    plot_data$method <- factor(plot_data$method, levels = c("Normal", "Time Randomised"))

    plot <- ggplot(plot_data, aes(x = method, y = point_est)) +
      geom_hline(yintercept = plot_data$point_est[1], color="grey60")+
      geom_point(data=bootstrapresults, aes(x=method, y=lambda),color="grey50",size=0.8, alpha=0.3)+
      geom_errorbar(aes(ymin = `5%`, ymax = `95%`), width = 0.2, color = "black") +
      geom_point(size = 2, color = "red3") +  # Point estimate
      scale_x_discrete(limits = c("Normal", "Time Randomised")) +
      annotate("label", label=paste0("p=",format(round(p_value, 4), nsmall = 4)), x=Inf, y=Inf, vjust=1, hjust=1)+
      scale_y_continuous(#transform = "log10"
      )+
      labs(title=title,
           x = NULL,
           y = "Rate") +
      theme_minimal()


    res <- list(confidence_intervals=ci, raw_results=bootstrapresults, plot=plot, p_value=p_value)

  } else {
    warning("Insufficient data points to fit distributions!")

    lowerres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)
    upperres <- data.frame(snp_threshold=NA,lambda=NA,k=NA,estimated_fp=NA)

    ci <- bind_rows(lowerres, upperres)
    res <- list(confidence_intervals=ci, raw_results=NA)
  }


}
