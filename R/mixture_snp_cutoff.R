library(tidyverse)

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
#' Estimates evolutionary rates from a mixed transmission dataset and an unrelated data set with a mixed probability distribution approach.Uses these rates to produce a threshold of SNP distances to be considered linked or not.
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param youden whether to produce SNP thresholds using the Youden method
#' @param threshold_range whether to produce a dataset of threshold considering a range of times (from 0.5 to 10 years)
#' @param max_time the maximum time(in days) utilised to calculate SNP thresholds, only applicable when time differences are provided
#' @param upper.tail percentile to calculate SNP thresholds
#' @param max_false_positive if the false positive rate from calculated threshold is higher than this value a warning is produced
#' @param trace trace parameter to pass to optim
#' @param start_params initial parametrs for optim, if NA (as default) will try 3 different start parameters and produce the highest likelyhood result
#'
#' @importFrom stats optim pnbinom qpois rnbinom
#' @importFrom dplyr filter
#'
#' @return SNP threshold, mutation rate, proportion related, estimated false positive and negative rate estimations
#'
#' @export
mixture_snp_cutoff <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                               youden=FALSE,threshold_range=FALSE, max_time= NA, upper.tail=0.95, max_false_positive=0.05, trace=FALSE, start_params= NA){

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
      nb_fit <- suppressWarnings(MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial"))

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


      if(anyNA(start_params)){
        result_attempts <- list(
          nlminb1 = nlminb(start=c(0.25,0.0001), objective=llk, x = trans_snp_dist,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace)),
          nlminb2 = nlminb(start=c(0.75,0.005), objective=llk, x = trans_snp_dist,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace)),
          nlminb3 = nlminb(start=c(0.5, 0.001), objective=llk, x = trans_snp_dist,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))
        )

        # Extract the best result
        result <- result_attempts[[which.min(sapply(result_attempts, `[[`, "objective"))]]
        best_result_name <- names(result_attempts)[[which.min(sapply(result_attempts, `[[`, "objective"))]]
      }
      else{
        result <- nlminb(start=start_params, objective=llk, x = trans_snp_dist, t = trans_time_dist, s = trans_sites,
                         lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))
        best_result_name <- "input param nlminb"
      }

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
        method="base",
        parameter_comb=best_result_name
      )
    } else {
      warning("Insufficient data points to fit distributions!")
      results <- tibble(
        snp_threshold=NA,
        lambda=NA,
        k=NA,
        estimated_fp=NA,
        estimated_fn=NA,
        method="failure",
        parameter_comb=NA

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
      nb_fit <- suppressWarnings(MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial"))

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


      if(anyNA(start_params)){
        result_attempts <- list(
          nlminb1 = nlminb(start=c(0.25,0.0001), objective=llk, x = trans_snp_dist, t = trans_time_dist,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace)),
          nlminb2 = nlminb(start=c(0.75,0.005), objective=llk, x = trans_snp_dist, t = trans_time_dist,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace)),
          nlminb3 = nlminb(start=c(0.5, 0.001), objective=llk, x = trans_snp_dist, t = trans_time_dist,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))
        )

        # Extract the best result
        result <- result_attempts[[which.min(sapply(result_attempts, `[[`, "objective"))]]
        best_result_name <- names(result_attempts)[[which.min(sapply(result_attempts, `[[`, "objective"))]]
      }
      else{
        result <- nlminb(start=start_params, objective=llk, x = trans_snp_dist, t = trans_time_dist,
                         lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))
        best_result_name <- "input param nlminb"
      }

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
        method="time",
        parameter_comb=best_result_name
      )

    } else {
      warning("Insufficient data points to fit distributions!")
      results <-
        tibble(
          snp_threshold=NA,
          lambda=NA,
          k=NA,
          estimated_fp=NA,
          method="failure",
          parameter_comb=NA

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
      nb_fit <- suppressWarnings(MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial"))

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

      if(anyNA(start_params)){
        result_attempts <- list(
          nlminb1 = nlminb(start=c(0.25,0.0001), objective=llk, x = trans_snp_dist, t = trans_time_dist, s = trans_sites,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace)),
          nlminb2 = nlminb(start=c(0.75,0.005), objective=llk, x = trans_snp_dist, t = trans_time_dist, s = trans_sites,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace)),
          nlminb3 = nlminb(start=c(0.5, 0.001), objective=llk, x = trans_snp_dist, t = trans_time_dist, s = trans_sites,
                           lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))
        )

        #finds best result
        result <- result_attempts[[which.min(sapply(result_attempts, `[[`, "objective"))]]
        best_result_name <- names(result_attempts)[[which.min(sapply(result_attempts, `[[`, "objective"))]]
      }
      else{
        result <- nlminb(start=start_params, objective=llk, x = trans_snp_dist, t = trans_time_dist, s = trans_sites,
               lower = c(0, 1e-10), upper = c(1, Inf), control = list(trace = trace))
        best_result_name <- "input param nlminb"
      }

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
        method="time+sites",
        parameter_comb=best_result_name
      )
    } else {
      warning("Insufficient data points to fit distributions!")
      results <-
        tibble(
          snp_threshold=NA,
          lambda=NA,
          k=NA,
          estimated_fp=NA,
          method="failure",
          parameter_comb=NA

        )
      threshold_range_df <- NA

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
