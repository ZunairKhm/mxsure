library(tidyverse)
library(ggplot2)
library(patchwork)

#### Log Sum Exp ####

log_sum_exp <- function(log_a, log_b) {
  # Ensure log_a is the max
  if (log_a < log_b) {
    tmp <- log_a
    log_a <- log_b
    log_b <- tmp
  }
  # Return the sum in log space
  return(log_a + log(1 + exp(log_b - log_a)))
}


#' Mixture SNP Cutoff
#'
#' Estimates evolutionary rates utilising a mixed transmission dataset and an unrelated data set with a mixed probability distribution approach.Uses these rates to produce a threshold of SNP distances to be considered linked or not.
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param youden whether to produce SNP thresholds using the Youden method
#' @param threshold_range whether to produce a dataset of threshold considering a range of times (from 0.5 to 10 years)
#' @param max_time the maximum time(in days) utilised to calculate SNP thresholds, only applicable when time differences are provided
#'
#' @return SNP threshold, mutation rate, proportion related, estimated false positive and negative rate estimations
#'
#' @export
mixture_snp_cutoff <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                               youden=FALSE,threshold_range=FALSE, max_time= NA, upper.tail=0.95, max_false_positive=0.05, trace=FALSE){

  #### Youden Cutoffs ####
  if(youden==TRUE){
    if ((length(trans_snp_dist) >= 10) && (length(unrelated_snp_dist) >= 10)){

      youden_index <- map_dbl(trans_snp_dist, ~{
        tp <- sum(trans_snp_dist <= .x)
        tn <- sum(unrelated_snp_dist > .x)
        fn <- sum(trans_snp_dist > .x)
        fp <- sum(unrelated_snp_dist <= .x)

        return(tp/(tp+fn) + tn/(tn+fp) - 1)
      })
      youden_snp_threshold <- trans_snp_dist[which.max(youden_index)]

      youden_results <-
        tibble(
          youden_snp_threshold=youden_snp_threshold,
          J=max(youden_index),
          youden_estimated_fp=sum(unrelated_snp_dist<=youden_snp_threshold)/length(unrelated_snp_dist)
        )


    } else {
      warning("Insufficient data points to call youden cutoff-point!")

      youden_results <- tibble(
        youden_snp_threshold=NA,
        J=NA,
        youden_estimated_fp=NA
      )

    }
  }


  #### threshold without time or sites considered ####
  if((anyNA(trans_time_dist) & anyNA(trans_sites))){


    if ((length(trans_snp_dist) >= 20) && (length(unrelated_snp_dist) >= 20)){
      nb_fit <- MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial")

      llk <- function(params, x){
        k <- params[[1]]
        lambda <- params[[2]]

        -sum(map_dbl(x, ~ {log_sum_exp(log(k) + dpois(x = .x,
                                                      lambda =  lambda,
                                                      log = TRUE),
                                       log(1-k) + dnbinom(x = .x,
                                                          size = nb_fit$estimate['size'],
                                                          mu = nb_fit$estimate['mu'],
                                                          log = TRUE))}))
      }

      start_params <- c(0.5, 1) # Initial guesses
      result <- optim(par = start_params, fn = llk, x = trans_snp_dist,
                      method = "L-BFGS-B", lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))

      snp_threshold <- qpois(upper.tail, lambda = result$par[[2]])


    if ((sum(unrelated_snp_dist<=snp_threshold)/length(unrelated_snp_dist)) > max_false_positive){
      warning(paste0("Inferred SNP threshold may have a false positive rate above ",
                     max_false_positive, "!"))
    }

    estimated_fp <- pnbinom(q = snp_threshold,
                            size = nb_fit$estimate['size'],
                            mu = nb_fit$estimate['mu'],
                            lower.tail = TRUE)


    #return results
    results <- tibble(
      snp_threshold=snp_threshold,
      lambda=result$par[[2]],
      k=result$par[[1]],
      estimated_fp=estimated_fp,
      method="base"
    )
    } else {
      warning("Insufficient data points to fit distributions!")
      results <- tibble(
        snp_threshold=NA,
        lambda=NA,
        k=NA,
        estimated_fp=NA,
        estimated_fn=NA,
        method="failure"

      )
    }

    if(youden==TRUE){
      return(list("results" = results,"youden" = youden_results))
    }  else {return(results)
    }


  }

  #### Threshold considering time but not sites
  if(!anyNA(trans_time_dist)&(anyNA(trans_sites))){
    if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){
      nb_fit <- MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial")

      llk <- function(params, x, t){
        k <- params[[1]]
        lambda <- params[[2]]

        -sum(map2_dbl(x, t, ~ {log_sum_exp(log(k) + dpois(x = .x,
                                                          lambda =  lambda*.y,
                                                          log = TRUE),
                                           log(1-k) + dnbinom(x = .x,
                                                              size = nb_fit$estimate['size'],
                                                              mu = nb_fit$estimate['mu'],
                                                              log = TRUE))}))
      }

      start_params <- c(0.5, 1) # Initial guesses
      result <- optim(par = start_params, fn = llk, x = trans_snp_dist, t=trans_time_dist,
                      method = "L-BFGS-B", lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))

      if (is.na(max_time)) {
        max_time <- 2*max(trans_time_dist)
      }

      snp_threshold <- qpois(0.95, lambda=result$par[[2]]*max_time) # qpois(1-max_false_positive, lambda=result$par[[2]]*max_time) #result$par[[2]]*quantile(trans_time_dist, 0.75))
      if(threshold_range==TRUE & !is.na(snp_threshold)){
        threshold_range_df <- data.frame(years=seq(0.5, 10, 0.5), threshold=NA, estimated_fp=NA, prop_pos=NA)
        #snp_threshold$thresholdexpected <- map(snp_threshold$years, ~{result$par[[2]]*365.25*.x*(max(trans_sites)/1000000)})
        threshold_range_df$threshold <- modify(threshold_range_df$years, ~{qpois(0.95, lambda=result$par[[2]]*365.25*.x)})
        threshold_range_df$estimated_fp <-modify(threshold_range_df$threshold, ~{sum(unrelated_snp_dist<=.x)/length(unrelated_snp_dist)})
        threshold_range_df$prop_pos <-modify(threshold_range_df$threshold, ~{sum(trans_snp_dist<=.x)/length(trans_snp_dist)})
      }

    if ((sum(unrelated_snp_dist<=snp_threshold)/length(unrelated_snp_dist)) > max_false_positive){
      warning(paste0("Inferred SNP threshold may have a false positive rate above ",
                     max_false_positive, "!"))
    }

    # results
    results <- tibble(
        snp_threshold=snp_threshold,
        lambda=result$par[[2]],
        k=result$par[[1]],
        estimated_fp=sum(unrelated_snp_dist<=snp_threshold)/length(unrelated_snp_dist),
        method="time"
      )

    } else {
      warning("Insufficient data points to fit distributions!")
      results <-
        tibble(
          snp_threshold=NA,
          lambda=NA,
          k=NA,
          estimated_fp=NA,
          method="failure"

        )
    }
#return results
    if(youden==TRUE & threshold_range==TRUE){
      return(list("results" = results, "youden" = youden_results, "thresholds" = threshold_range_df))
    }
    if (youden==FALSE & threshold_range==TRUE){
      return(list("results" =results, "thresholds" =threshold_range_df))
    }
    if(youden==TRUE & threshold_range==FALSE){
      return(list("results" =results, "youden" =youden_results))
    }
    if(youden==FALSE & threshold_range==FALSE){
      return(results)
    }



  }

  #### Threshold considering time and sites #####
  if(!anyNA(trans_time_dist)&(!anyNA(trans_sites))){
  if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){
    nb_fit <- MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial")

    llk <- function(params, x, t, s){
      k <- params[[1]]
      lambda <- params[[2]]

      -sum(pmap_dbl(list(x, t,s), ~ {log_sum_exp(log(k) + dpois(x = ..1,
                                                                lambda =  lambda*..2*(..3/1000000), #gives rate esimate per day time per million bp
                                                                log = TRUE),
                                                 log(1-k) + dnbinom(x = ..1,
                                                                    size = nb_fit$estimate['size'],
                                                                    mu = nb_fit$estimate['mu'],
                                                                    log = TRUE))}))
    }

    start_params <- c(0.5, 1) # Initial guesses
    result <- optim(par = start_params, fn = llk, x = trans_snp_dist, t=trans_time_dist, s=trans_sites,
                    method = "L-BFGS-B", lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))

    if (is.na(max_time)) {
      max_time <- 2*max(trans_time_dist)
    }

    snp_threshold <- qpois(0.95, lambda=result$par[[2]]*max_time*(mean(trans_sites)/1000000))

    if(threshold_range==TRUE & !is.na(snp_threshold)){
      threshold_range_df <- data.frame(years=seq(0.5, 10, 0.5), threshold=NA, estimated_fp=NA, prop_pos=NA)
    #snp_threshold$thresholdexpected <- map(snp_threshold$years, ~{result$par[[2]]*365.25*.x*(max(trans_sites)/1000000)})
    threshold_range_df$threshold <- modify(threshold_range_df$years, ~{qpois(0.95, lambda=result$par[[2]]*365.25*.x*(mean(trans_sites)/1000000))})
    threshold_range_df$estimated_fp <-modify(threshold_range_df$threshold, ~{sum(unrelated_snp_dist<=.x)/length(unrelated_snp_dist)})
    threshold_range_df$prop_pos <-modify(threshold_range_df$threshold, ~{sum(trans_snp_dist<=.x)/length(trans_snp_dist)})
    }

  if ((sum(unrelated_snp_dist<=snp_threshold)/length(unrelated_snp_dist)) > max_false_positive){
    warning(paste0("Inferred SNP threshold may have a false positive rate above ",
                   max_false_positive, "!"))
  }

  # results
  results <-  tibble(
      snp_threshold=snp_threshold,
      lambda=result$par[[2]],
      k=result$par[[1]],
      estimated_fp=sum(unrelated_snp_dist<=snp_threshold)/length(unrelated_snp_dist),
      method="time+sites"
    )
  } else {
    warning("Insufficient data points to fit distributions!")
    results <-
      tibble(
        snp_threshold=NA,
        lambda=NA,
        k=NA,
        estimated_fp=NA,
        method="failure"

      )
  }
# returning results
  if(youden==TRUE & threshold_range==TRUE){
    return(list("results" = results, "youden" = youden_results, "thresholds" = threshold_range_df))
  }
  if (youden==FALSE & threshold_range==TRUE){
    return(list("results" =results, "thresholds" =threshold_range_df))
  }
  if(youden==TRUE & threshold_range==FALSE){
    return(list("results" =results, "youden" =youden_results))
  }
  if(youden==FALSE & threshold_range==FALSE){
    return(results)
  }
  }


}


#' Bootstrapped confidence intervals
#'
#' Creates confidence intervals for outputs from mixture_snp_cutoff using bootstrapping
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param sample_size size of each bootstrap sample
#' @param sample_n number of bootstrap sampling to conduct
#'
#' @return Confidence intervals
#'
#' @export
mixture_snp_cutoff_ci <- function(trans_snp_dist,unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                                  sample_size=length(trans_snp_dist), sample_n=1000, confidence_level=0.95){

  mix_data <- tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist, sites=trans_sites)

  bootstrapresults <- map_dfr(1:sample_n, ~{
    x <- slice_sample(mix_data, n= sample_size, replace = TRUE)
    y <- mixture_snp_cutoff(x$snp_dist,unrelated_snp_dist, x$time_dist, x$sites)
    y[1:4]
  },.progress=TRUE)

  lowerres <- bootstrapresults|>
    summarise(across(everything(),  ~quantile(.x, 1-confidence_level)))
  upperres <- bootstrapresults|>
    summarise(across(everything(), ~quantile(.x, confidence_level))
    )

  res <- bind_rows(lowerres, upperres)
  return(res)
}

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

#' Simulate Mixed related SNP distance dataset
#'
#' @param lambda mutation rate input (SNPs per day)
#' @param k proportion of mixed dataset that is linked
#' @param nbmu mean of the negative binomial distribution underlying the unrelated SNP distances
#' @param error_param chance that each SNP is missed
#' @param n size of dataset to simulate
#'
#' @return SNP distance dataset with a mixture of related and unrelated pairs with time differences
#' @export
#'
#' @examples
#' #simulates dataset with an E. coli profile
#' x <- simulate_mixsnp_data(lambda=0.00179*4, 0.5, 500, error_param=NA)
simulate_mixsnp_data <- function(lambda, k,nbmu=500, error_param=NA, n=1000){
  mix_snp_dist <- map_dfr(1:n, ~{

    tt <- rexp(1, rate = 0.02) #time distribution
    if (runif(1)<k){
      dd <- rpois(n = 1, lambda = tt*lambda)
      rr <- "Related"
    } else {
      dd <- rnbinom(1, mu = nbmu, size = 1)
      rr <- "Unrelated"
    }
    if (!is.na(error_param)){
      dd <- rbinom(1, dd, 1-error_param) #error param is chance each SNP difference is missed
    }
    return(
      tibble(snp_dist=dd, time_dist=tt, relation=rr)
    )

  })
}

#' Simulate unrelaed SNP distance dataset
#'
#' @param nbmu mean of the underlying negative binomial distribution
#' @param nbsize size parameter of the underlying negative binmial distribution (not the size of the simulated dataset)
#' @param n size/length of the simulated dataset
#'
#' @return unrelated SNP distance dataset
#' @export
#'
#' @examples
#' y <- simulate_unrelsnp_data()
simulate_unrelsnp_data <- function(nbmu=500, nbsize=1, n=1000){
  rnbinom(n, mu=nbmu, size=nbsize)
}



