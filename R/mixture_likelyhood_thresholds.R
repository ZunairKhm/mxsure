#' Mixture Likelyhood Thresholds
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param year_range range of years to calculate thresholds over
#'
#' @return Threshold corresponding to the highest SNP distance where it is more likely to be related than unrelated
#' @export
#'
#' @examples
mixture_likelyhood_thresholds <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites, truncation_point=2000){
  if ((length(trans_snp_dist) >= 30) && (length(unrelated_snp_dist) >= 30)){

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


    result <-  distinct(tibble(snp_dist=trans_snp_dist, time_dist=trans_time_dist))
    result <- result |>
      mutate(rel_lh = (dpois(snp_dist, (mix_res$lambda*(time_dist/365.25)*mean(trans_sites)+mix_res$intercept)) / ppois(truncation_point, (mix_res$lambda*(time_dist/365.25)*mean(trans_sites)+mix_res$intercept))),
             unrel_lh = (dnbinom(snp_dist, mu = nb_fit$estimate['mu'], size = nb_fit$estimate['size']) / pnbinom(truncation_point, mu = nb_fit$estimate['mu'], size = nb_fit$estimate['size'])))|>
    mutate(LR = rel_lh/unrel_lh)

    rm(ptruncnbinom, envir = .GlobalEnv)
    rm(dtruncnbinom, envir = .GlobalEnv)
    rm(qtruncnbinom, envir = .GlobalEnv)

    return(result)

  }else {
    warning("Insufficient data points to fit distributions!")
    return(tibble(NA))
  }

}
