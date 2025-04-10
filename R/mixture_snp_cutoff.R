library(tidyverse)
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
#' @param max_time the maximum time(in days) utilised to calculate SNP thresholds, only applicable when time differences are provided, if not provided will utilise maximum of supplied times
#' @param upper.tail percentile to calculate SNP thresholds
#' @param max_false_positive if the false positive rate from calculated threshold is higher than this value a warning is produced
#' @param trace trace parameter to pass to nlminb
#' @param start_params initial parameters for lambda, k, and intercept parameters for optimisation. If NA (as default) will try a range of different start parameters and produce the highest likelyhood result
#' @param truncation_point a SNP distance limit for the data, if set to NA will estimate as if there is no limit. Will be faster with a lower truncation point.
#' @param prior_lambda parameters for a gamma prior distribution for rate estimation, if set to "default" will use default parameters
#' @param prior_k parameters for a beta prior distribution for related proportion estimation, if set to "default" will use default parameters
#' @param lambda_bounds bounds of rate estimation in SNPs/year/site if given time and site data
#' @param k_bounds bounds of related proportion estimation
#'
#' @importFrom stats optim pnbinom qpois rnbinom
#' @importFrom dplyr filter
#' @importFrom fitdistrplus fitdist
#'
#' @return SNP threshold, mutation rate, proportion related, estimated false positive and negative rate estimations
#'
#' @export
mixture_snp_cutoff <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist=NA, trans_sites=NA,
                                  truncation_point=2000,
                               youden=FALSE,threshold_range=FALSE, max_time= NA,
                               prior_lambda=NA, prior_k=NA, lambda_bounds=c(1e-10, 1), k_bounds=c(0,1), intercept_bounds=c(-Inf, Inf),
                               upper.tail=0.95, max_false_positive=0.05, trace=FALSE, start_params= NA){

  #correction to convert to snp/day(/site)
  if(!anyNA(trans_time_dist)){
    lambda_bounds <- lambda_bounds/365.25
  }

  #### Log Sum Exp ####

  log_sum_exp <- function(log_a, log_b) {
    if (anyNA(c(log_a, log_b))){
      return(-Inf)
    }
    # Ensure log_a is the max
     if (log_a < log_b) {
      tmp <- log_a
      log_a <- log_b
      log_b <- tmp
    }
    # Return the sum in log space
    return(log_a + log(1 + exp(log_b - log_a)))
  }

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


  #truncating data
  if(is.na(truncation_point)){
    truncation_point <- Inf
  }
  unrelated_snp_dist_orig <- unrelated_snp_dist
  x <- tibble(trans_snp_dist, trans_time_dist, trans_sites)
  x <- filter(x, trans_snp_dist<truncation_point)
  trans_snp_dist <- x$trans_snp_dist
  trans_time_dist <- x$trans_time_dist
  trans_sites <- x$trans_sites
  unrelated_snp_dist <- unrelated_snp_dist[unrelated_snp_dist<truncation_point]

  #setting default priors, not on by default
  if(length(prior_lambda)==1){
    if(is.element(prior_lambda, "default")){
      prior_lambda <- c(1.06, 0.44)
    }}
  if(length(prior_k)==1){
    if(is.element(prior_k, "default")){
      prior_k <- c(1.60, 1.25)
    }}

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
    warning("No time data inputted, rate will still be estimated but consider whether this is appropriate")

    if ((length(trans_snp_dist) >= 20) && (length(unrelated_snp_dist) >= 20)){
      #distant data fitting
      # distant dataset fitting
      m <- mean(unrelated_snp_dist)
      v <- var(unrelated_snp_dist)
      size <- if (v > m) {
        m^2/(v - m)
      }else{100}

      nb_fit <- fitdistrplus::fitdist(unrelated_snp_dist, dist="truncnbinom", start=list(mu=m, size=size))

      #mixed data fitting
      llk <- function(params, x){
        k <- params[[1]]
        lambda <- params[[2]]
        intercept <- params[[3]]

        -sum(pmap_dbl(list(x), ~ {suppressWarnings(log_sum_exp(log(k) + dpois(x = ..1,
                                                                  lambda =  lambda + intercept, #gives rate esimate per average time of the dataset
                                                                  log = TRUE) -
                                                     ppois(truncation_point,
                                                           lambda =  lambda + intercept,
                                                           log = TRUE ),
                                                   log(1-k) + dnbinom(x = ..1,
                                                                      size = nb_fit$estimate["size"],
                                                                      mu = nb_fit$estimate["mu"],
                                                                      log = TRUE)-
                                                     pnbinom(truncation_point,
                                                             size = nb_fit$estimate["size"],
                                                             mu = nb_fit$estimate["mu"],
                                                             log = TRUE)))
        }))
      }


      if(anyNA(start_params)){
        # Define parameter grid
        start_vals <- expand.grid(k = c(0.25, 0.5, 0.75), lambda = c(0.0001, 0.001, 0.01), intercept = c(0))

        # Run nlminb for each combination
        result_attempts <- pmap(list(start_vals$k, start_vals$lambda, start_vals$intercept),
                                function(k, lambda, intercept) {
                                  nlminb(
                                    start = c(k, lambda, intercept),
                                    objective = llk,
                                    x = trans_snp_dist,
                                    lower = c(k_bounds[1], lambda_bounds[1], intercept_bounds[1]),
                                    upper = c(k_bounds[2], lambda_bounds[2], intercept_bounds[2]),
                                    control = list(trace = trace)
                                  )})

        # Optionally, name the list elements
        names(result_attempts) <- paste0("nlminb", seq_along(result_attempts))

        #finds best result
        result <- result_attempts[[which.min(sapply(result_attempts, `[[`, "objective"))]]
        best_result_name <- names(result_attempts)[[which.min(sapply(result_attempts, `[[`, "objective"))]]
      }
      else{
        result <- nlminb(start=start_params, objective=llk, x = trans_snp_dist,,
                         lower = c(k_bounds[1], lambda_bounds[1], intercept_bounds[1]),
                         upper = c(k_bounds[2], lambda_bounds[2], intercept_bounds[2]),
                         control = list(trace = trace))
        best_result_name <- "input param nlminb"
      }

      snp_threshold <- qpois(upper.tail, lambda = result$par[[2]]*2)


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
        intercept=result$par[[3]],
        estimated_fp=estimated_fp,
        lambda_units="SNPs per average sampling time per genome",
        parameter_comb=best_result_name
      )
    } else {
      warning("Insufficient data points to fit distributions!")
      results <- tibble(
        snp_threshold=NA,
        lambda=NA,
        k=NA,
        intercept=NA,
        estimated_fp=NA,
        estimated_fn=NA,
        lambda_units=NA,
        parameter_comb=NA

      )
    }

    if(youden==TRUE){
      return(list("results" = results,"youden" = youden_results))
    }  else {return(results)
    }


  }

  #### Threshold considering time but not sites ####
  if(!anyNA(trans_time_dist)&(anyNA(trans_sites))){
    if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){

      # distant dataset fitting
      m <- mean(unrelated_snp_dist)
      v <- var(unrelated_snp_dist)
      size <- if (v > m) {
        m^2/(v - m)
      }else{100}

      nb_fit <- fitdistrplus::fitdist(unrelated_snp_dist, dist="truncnbinom", start=list(mu=m, size=size))

      #mixed data fitting
      llk <- function(params, x, t){
        k <- params[[1]]
        lambda <- params[[2]]
        intercept <- params[[3]]

        -sum(pmap_dbl(list(x, t), ~ {suppressWarnings(log_sum_exp(log(k) + dpois(x = ..1,
                                                                  lambda =  lambda*(..2) + intercept, #gives rate esimate per day
                                                                  log = TRUE)
                                                 -
                                                     ppois(truncation_point,
                                                           lambda =  lambda*..2 + intercept,
                                                           log = TRUE )
                                                 ,
                                                   log(1-k) + dnbinom(x = ..1,
                                                                      size = nb_fit$estimate["size"],
                                                                      mu = nb_fit$estimate["mu"],
                                                                      log = TRUE)
                                                 -
                                                     pnbinom(truncation_point,
                                                             size = nb_fit$estimate["size"],
                                                             mu = nb_fit$estimate["mu"],
                                                             log = TRUE)
                                                 ))+
                                         ifelse(!anyNA(prior_k),  dbeta(k, prior_k[1], prior_k[2], log = TRUE), 0)+
                                         ifelse(!anyNA(prior_lambda), dgamma(lambda, prior_lambda[1], prior_lambda[2], log = TRUE),0)
          }))
      }

      if(anyNA(start_params)){
        # Define parameter grid
        start_vals <- expand.grid(k = c(0.25, 0.5, 0.75), lambda = c(0.01, 0.1, 1), intercept = c(0))

        # Run nlminb for each combination
        result_attempts <- pmap(list(start_vals$k, start_vals$lambda, start_vals$intercept),
                                function(k, lambda, intercept) {
                                  nlminb(
                                    start = c(k, lambda, intercept),
                                    objective = llk,
                                    x = trans_snp_dist,
                                    t = trans_time_dist,
                                    lower = c(k_bounds[1], lambda_bounds[1], intercept_bounds[1]),
                                    upper = c(k_bounds[2], lambda_bounds[2], intercept_bounds[2]),
                                    control = list(trace = trace)
                                  )})


        # Extract the best result
        result <- result_attempts[[which.min(sapply(result_attempts, `[[`, "objective"))]]

      }
      else{
        result <- nlminb(start=c(start_params, 0), objective=llk, x = trans_snp_dist, t = trans_time_dist,
                         lower = c(k_bounds[1], lambda_bounds[1], intercept_bounds[1]),
                         upper = c(k_bounds[2], lambda_bounds[2], intercept_bounds[2]),
                         control = list(trace = trace))

      }

      if (is.na(max_time)) {
        max_time <- max(trans_time_dist)
      }

      snp_threshold <- qpois(upper.tail, lambda=result$par[[2]]*(max_time))
      if(threshold_range==TRUE & !is.na(snp_threshold)){
        threshold_range_df <- data.frame(years=seq(0.5, 10, 0.5), threshold=NA, estimated_fp=NA, prop_pos=NA)
        threshold_range_df$threshold <- modify(threshold_range_df$years, ~{qpois(upper.tail, lambda=result$par[[2]]*365.25*.x)})
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
        lambda=result$par[[2]]*365.25,
        k=result$par[[1]],
        intercept=result$par[[3]],
        estimated_fp=sum(unrelated_snp_dist<=snp_threshold)/length(unrelated_snp_dist),
        lambda_units="SNPs per year per genome"

      )

    } else {
      warning("Insufficient data points to fit distributions!")
      results <-
        tibble(
          snp_threshold=NA,
          lambda=NA,
          k=NA,
          intercept=NA,
          estimated_fp=NA,
          lambda_units=NA

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

      # distant dataset fitting
      m <- mean(unrelated_snp_dist)
      v <- var(unrelated_snp_dist)
      size <- if (v > m) {
        m^2/(v - m)
      }else{100}

      nb_fit <- fitdistrplus::fitdist(unrelated_snp_dist, dist="truncnbinom", start=list(mu=m, size=size))

      #mixed dataset fitting
      llk <- function(params, x, t, s){
        k <- params[[1]]
        lambda <- params[[2]]
        intercept <- params[[3]]

        -sum(pmap_dbl(list(x, t, s), ~ {suppressWarnings(log_sum_exp(log(k) + dpois(x = ..1,
                                                                  lambda =  lambda*(..2)*..3 + intercept, #gives rate esimate per day per bp
                                                                  log = TRUE) -
                                                     ppois(truncation_point,
                                                           lambda =  lambda*..2*(..3) + intercept,
                                                           log = TRUE ),
                                                   log(1-k) + dnbinom(x = ..1,
                                                                      size = nb_fit$estimate["size"],
                                                                      mu = nb_fit$estimate["mu"],
                                                                      log = TRUE)-
                                                     pnbinom(truncation_point,
                                                             size = nb_fit$estimate["size"],
                                                             mu = nb_fit$estimate["mu"],
                                                             log = TRUE)))
                                      # +
                                      # ifelse(!anyNA(prior_k),  dbeta(k, prior_k[1], prior_k[2], log = TRUE), 0)+
                                      # ifelse(!anyNA(prior_lambda), dgamma(lambda, prior_lambda[1], prior_lambda[2], log = TRUE),0)
          }))
      }

      if(anyNA(start_params)){
        # Define parameter grid
        start_vals <- expand.grid(k = c(0.5, 0.75, 0.9), lambda = c(1e-11, 1e-10, 1e-9, 1e-8), intercept=c(0))

        # Run nlminb for each combination
        result_attempts <- pmap(list(start_vals$k, start_vals$lambda, start_vals$intercept),
                                function(k, lambda, intercept) {
                                  nlminb(
                                    start = c(k, lambda, intercept),
                                    objective = llk,
                                    x = trans_snp_dist,
                                    t = trans_time_dist,
                                    s = trans_sites,
                                    lower = c(k_bounds[1], lambda_bounds[1], intercept_bounds[1]),
                                    upper = c(k_bounds[2], lambda_bounds[2], intercept_bounds[2]),
                                    control = list(trace = trace)
                                  )})



        #finds best result
        result <- result_attempts[[which.min(sapply(result_attempts, `[[`, "objective"))]]

      }
      else{
        result <- nlminb(start=c(start_params), objective=llk, x = trans_snp_dist, t = trans_time_dist, s = trans_sites,
                         lower = c(k_bounds[1], lambda_bounds[1], intercept_bounds[1]),
                         upper = c(k_bounds[2], lambda_bounds[2], intercept_bounds[2]),
                         control = list(trace = trace))

      }

      if (is.na(max_time)) {
        max_time <- max(trans_time_dist)
      }

      snp_threshold <- qpois(upper.tail, lambda=result$par[[2]]*(max_time)*(mean(trans_sites)))

      if(threshold_range==TRUE & !is.na(snp_threshold)){
        threshold_range_df <- data.frame(years=seq(0.5, 10, 0.5), threshold=NA, estimated_fp=NA, prop_pos=NA)
        threshold_range_df$threshold <- modify(threshold_range_df$years, ~{qpois(upper.tail, lambda=result$par[[2]]*.x*365.25*(mean(trans_sites)))})
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
        lambda=result$par[[2]]*365.25,
        k=result$par[[1]],
        intercept=result$par[[3]],
        estimated_fp=sum(unrelated_snp_dist<=snp_threshold)/length(unrelated_snp_dist),
        lambda_units="SNPs per year per site",
        convergence=result$convergence,
        message=result$message,
        iterations=result$iterations
      )
    } else {
      warning("Insufficient data points to fit distributions!")
      results <-
        tibble(
          snp_threshold=NA,
          lambda=NA,
          k=NA,
          intercept=NA,
          estimated_fp=NA,
          lambda_units=NA

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

  rm(ptruncnbinom, envir = .GlobalEnv)
  rm(dtruncnbinom, envir = .GlobalEnv)
  rm(qtruncnbinom, envir = .GlobalEnv)

}
