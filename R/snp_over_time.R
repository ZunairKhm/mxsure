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
snp_over_time <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites, truncation_point=2000, title="SNPs over Time", jitter=TRUE, p_value=NA, ci_data=NA, time_limits=c(0,NA), predictive_interval=TRUE){


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

  if(!is.na(mean(trans_sites, na.rm=TRUE))){
    trans_sites <- 1
  }

  LH <-  distinct(tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist))
  LH <- LH |>
    mutate(rel_lh = (dpois(snp_dist, (mix_res$lambda*(time_dist/365.25)*mean(trans_sites)+mix_res$intercept)) / ppois(truncation_point, (mix_res$lambda*(time_dist/365.25)*mean(trans_sites)+mix_res$intercept))),
           unrel_lh = (dnbinom(snp_dist, mu = nb_fit$estimate['mu'], size = nb_fit$estimate['size']) / pnbinom(truncation_point, mu = nb_fit$estimate['mu'], size = nb_fit$estimate['size'])))|>
    mutate(LR = rel_lh/unrel_lh)

  rm(ptruncnbinom, envir = .GlobalEnv)
  rm(dtruncnbinom, envir = .GlobalEnv)
  rm(qtruncnbinom, envir = .GlobalEnv)


  data <- data.frame(snp_dist=trans_snp_dist, time_dist=trans_time_dist)
  data <- filter(data, trans_snp_dist<=mix_res$snp_threshold+1)
  if(jitter==TRUE){
  data$snp_dist <- abs(jitter(data$snp_dist))
  data$time_dist <- abs(jitter(data$time_dist))
  }

    lambda <- (mix_res$lambda*mean(trans_sites))
    predictive_intervals <- tibble(data$time_dist=0:max(c(max(data$time_dist), time_limits[2]), na.rm=TRUE))
    predictive_intervals <- predictive_intervals|>
      mutate(estimate=qpois(0.5, data$time_dist*mix_res$lambda))


  if(!anyNA(ci_data)){
    ci[1] <- (ci$confidence_intervals$lambda[1]*mean(sites))
    ci[2] <- (ci$confidence_intervals$lambda[2]*mean(sites))
    }

    predictive_intervals <- predictive_intervals|>
      mutate(low_ci=qpois(0.025, (data$time_dist/365.25)*mix_res$lambda),
             high_ci=qpois(0.975, (data$time_dist/365.25)*mix_res$lambda))

  if(predictive_intervals){
  ggplot(data, aes(x=time_dist, y=snp_dist))+
    scale_y_continuous(limits = c(0, mix_res$snp_threshold+1), expand = c(0.01,0.01))+
    scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
    geom_point()+
    # geom_abline(intercept=0, slope = lambda)
    # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[1], NA), linetype="dotted")+
    # geom_abline(intercept=0, slope = ifelse(!anyNA(ci), ci[2], NA), linetype="dotted")+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=estimate))+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=low_ci), linetype="dotted")+
    geom_step(data = predictive_intervals, aes(x=time_dist, y=high_ci), linetype="dotted")+
    geom_abline(intercept=snp_threshold, slope=0, linetype="dashed", alpha=0.75)+
    geom_label(x=Inf, y=Inf, vjust=1, hjust=1,
               label = paste0(signif(lambda, digits = 3), " SNPs per year",
                              (ifelse(anyNA(ci),"", paste0("\n95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                              (ifelse(anyNA(p_value),"", paste0("\np-value=", format(round(p_value, 4)))))
                              ),
               size=5)+
    labs(title=title)+
    theme_minimal()
  } else {
    ggplot(data, aes(x=time_dist, y=snp_dist))+
      scale_y_continuous(limits = c(0, mix_res$snp_threshold+1), expand = c(0.01,0.01))+
      scale_x_continuous(limits = c(0, NA), expand = c(0.01,0.01))+
      geom_point()+
      geom_abline(intercept=mix_res$intercept, slope = mix_res$lambda)+
      geom_abline(intercept=snp_threshold, slope=0, linetype="dashed", alpha=0.75)+
      geom_label(x=Inf, y=Inf, vjust=1, hjust=1,
                 label = paste0(signif(lambda, digits = 3), " SNPs per year",
                                (ifelse(anyNA(ci),"", paste0("\n95% CI: ", signif(ci[1], 3), ", ", signif(ci[2], 3)))),
                                (ifelse(anyNA(p_value),"", paste0("\np-value=", format(round(p_value, 4)))))
                 ),
                 size=5)+
      labs(title=title)+
      theme_minimal()
  }
    }
