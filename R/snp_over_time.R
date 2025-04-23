library(ggplot2)

#' SNP dist/time plot
#'
#' Creates a plot of snp distance over time
#'
#' @param snp_dist list of SNP distances from a mixed transmission data set
#' @param time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
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
snp_over_time <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites, truncation_point=2000, title="SNPs over Time", jitter=TRUE, p_value=NA, ci_data=NA, time_limits=c(0,NA), under_threshold=TRUE){

  unrelated_snp_dist <- unrelated_snp_dist[unrelated_snp_dist<truncation_point]

  #defining distribution functions for the truncated negative binomial distr, needs to be specific to the truncation point and in the global environment for fitdistplus
  dtruncnbinom <<- function(x, mu, size){
    dnbinom(x = x, size = size, mu = mu) / pnbinom(truncation_point, size = size, mu = mu)
  }
  ptruncnbinom <<- function(q, mu, size){
    pnbinom(q = q, size = size, mu = mu) / pnbinom(truncation_point, size = size, mu = mu)
  }
  qtruncnbinom <<- function(p, mu, size){
    qnbinom(p = p*(pnbinom(truncation_point, size=size, mu=mu)), size = size, mu = mu)
  }


  # distant dataset fitting
  m <- mean(unrelated_snp_dist)
  v <- var(unrelated_snp_dist)
  size <- if (v > m) {
    m^2/(v - m)
  }else{100}

  nb_fit <- fitdistrplus::fitdist(unrelated_snp_dist, dist="truncnbinom", start=list(mu=m, size=size))
  mix_res <- suppressWarnings(mixture_snp_cutoff(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites))

  if(is.na(mean(trans_sites, na.rm=TRUE))){
    trans_sites <- 1
  }

  LH <-  distinct(tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist))
  LH <- LH |>
    mutate(rel_lh = (dpois(snp_dist, (mix_res$lambda*(time_dist/365.25)*mean(trans_sites)+mix_res$intercept), log=TRUE) #/ ppois(truncation_point, (mix_res$lambda*(time_dist/365.25)*mean(trans_sites)+mix_res$intercept))
                     ),
           unrel_lh = (dnbinom(snp_dist, mu = nb_fit$estimate['mu'], size = nb_fit$estimate['size'], log = TRUE) - pnbinom(truncation_point, mu = nb_fit$estimate['mu'], size = nb_fit$estimate['size'], log=TRUE)))|>
    mutate(LHR = (rel_lh-unrel_lh))

  #return(LH)

  rm(ptruncnbinom, envir = .GlobalEnv)
  rm(dtruncnbinom, envir = .GlobalEnv)
  rm(qtruncnbinom, envir = .GlobalEnv)


  data <- data.frame(snp_dist=trans_snp_dist, time_dist=trans_time_dist)
  if(under_threshold){
  data <- filter(data, trans_snp_dist<=mix_res$snp_threshold+1)
  } else {
    data <- filter(data, trans_snp_dist<truncation_point)
  }
  data$LHR <- LH$LHR[match(paste(data$snp_dist, data$time_dist), paste(LH$snp_dist, LH$time_dist))]
  if(jitter==TRUE){
  data$snp_dist <- abs(jitter(data$snp_dist))
  data$time_dist <- abs(jitter(data$time_dist))
  }
    lambda <- mix_res$lambda*mean(trans_sites)
    predictive_intervals <- tibble(time_dist=0:max(c(max(data$time_dist), time_limits[2]), na.rm=TRUE))
    predictive_intervals <- predictive_intervals|>
      mutate(estimate=qpois(0.5, (time_dist/365.25)*lambda))

  if(!anyNA(ci_data)){
    ci <- c(ci_data$confidence_intervals$lambda[1]*mean(trans_sites),
     ci_data$confidence_intervals$lambda[2]*mean(trans_sites))
    }

    predictive_intervals <- predictive_intervals|>
      mutate(low_ci=qpois(0.025, (time_dist/365.25)*lambda),
             high_ci=qpois(0.975, (time_dist/365.25)*lambda))



  if(under_threshold){
  ggplot(data, aes(x=time_dist, y=snp_dist, color=LHR))+
    scale_y_continuous(limits = c(0, mix_res$snp_threshold+1), expand = c(0.01,0.01))+
    scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
      scale_color_gradientn(
        colors = viridis::cividis(256),
        limits=c(-2, 2),
        values = scales::rescale(c(-2, 0 ,2)),
        oob=scales::squish
      ) +
    geom_point()+
    # geom_abline(intercept=0, slope = lambda)
    # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[1], NA), linetype="dotted")+
    # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[2], NA), linetype="dotted")+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=estimate), color="black")+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=low_ci), color="black", linetype="dotted")+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=high_ci), color="black", linetype="dotted")+
      geom_hline(yintercept=mix_res$snp_threshold, linetype="dashed", alpha=0.75, color="black")+
    geom_label(x=Inf, y=Inf, vjust=1, hjust=1,color= "black",
               label = paste0(signif(lambda, digits = 3), " SNPs per year",
                              (ifelse(anyNA(ci_data),"", paste0("\n95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                              (ifelse(anyNA(p_value),"", paste0("\np-value=", format(round(p_value, 4)))))
                              ),
               size=5)+
    labs(title=title,
         y="SNPs",
         x="Time (Days)",
         color="logLHR",
         subtitle=paste0(signif(lambda, digits = 3), " SNPs per year",
                         (ifelse(anyNA(ci_data),"", paste0("; 95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                         (ifelse(anyNA(p_value),"", paste0("; p-value=", format(round(p_value, 4)))))
         ))+
    theme_minimal()
  } else {
    ggplot(data, aes(x=time_dist, y=snp_dist, color=LHR))+
      scale_y_continuous(limits = c(0, NA), expand = c(0.01,0.01), transform = "pseudo_log", breaks = c(signif(exp(seq(0, log(truncation_point), length.out=10)), 1)))+
      scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
      scale_color_gradientn(
        colors = viridis::viridis(256),
        limits=c(-2, 2),
        values = scales::rescale(c(-2, 0, 2)),
        oob=scales::squish
      ) +
      geom_point()+
      # geom_abline(intercept=0, slope = lambda)
      # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[1], NA), linetype="dotted")+
      # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[2], NA), linetype="dotted")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=estimate), color="black")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=low_ci), color="black", linetype="dotted")+
      geom_step(data = predictive_intervals, aes(x=time_dist, y=high_ci), color="black", linetype="dotted")+
      geom_hline(yintercept=mix_res$snp_threshold, linetype="dashed", alpha=0.75, color="black")+
      labs(title=title,
           y="SNPs",
           x="Time (Days)",
           color="logLHR",
           subtitle=paste0(signif(lambda, digits = 3), " SNPs per year",
                           (ifelse(anyNA(ci_data),"", paste0("; 95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                           (ifelse(anyNA(p_value),"", paste0("; p-value=", format(round(p_value, 4)))))
           ))+
      theme_minimal()
  }
    }
