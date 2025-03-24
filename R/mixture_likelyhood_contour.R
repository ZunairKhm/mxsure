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

negllk <- function(k, lambda, x, t, s, nbfitmu, nbfitsize,truncation_point, prior_k, prior_lambda){

  -log(-sum(pmap_dbl(list(x, t,s), ~ {log_sum_exp(log(k) + dpois(x = ..1,
                                                                 lambda =  lambda*..2*(..3/1000000), #gives rate esimate per day time per million bp
                                                                 log = TRUE) -
                                                    ppois(truncation_point,
                                                          lambda =  lambda*..2*(..3/1000000),
                                                          log = TRUE ),
                                                  log(1-k) + dnbinom(x = ..1,
                                                                     size = nbfitsize,
                                                                     mu = nbfitmu,
                                                                     log = TRUE)-
                                                    pnbinom(truncation_point,
                                                            size = nbfitsize,
                                                            mu = nbfitmu,
                                                            log = TRUE))+
      ifelse(!anyNA(prior_k),  dbeta(k, prior_k[1], prior_k[2], log = TRUE), 0)+
      ifelse(!anyNA(prior_lambda), dgamma(lambda, prior_lambda[1], prior_lambda[2], log = TRUE),0)
  })))
}

#' Contour plot for likelyhood estimation of mixture distribution mutation rate estimation
#'
#' @param trans_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param trans_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param trans_sites list of sites considered for each SNP distance in mixed data set
#' @param resolution resolution of plot (main computation limitation)
#' @param lambda_limits bounds of lambda axis (plots on a log10 scale)
#' @param k_limits bounds of k axis
#' @param bins number of bins to pass to ggplot contour function
#' @param title title to pass to ggplot title
#'
#' @return ggplot contour map of likelyhood space
#' @export
#'
#' @examples
#'
mixture_likelyhood_contour <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites, threshold=NA, truncation_point=NA,  prior_lambda=NA, prior_k=NA,
                                       resolution=100, lambda_limits=c(0.0001, 1), k_limits=c(0.2,0.99), estimate=NULL, low_ci=NULL, high_ci=NULL,
                                       bins=NULL, title="Mixture Likelihood Contour Plot"){

  #defining distribution functions for the truncated negative binomial distr
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

  #setting priors
  if(length(prior_lambda)==1){
    if(is.element(prior_lambda, "default")){
      prior_lambda <- c(1.06, 0.44)
    }}
  if(length(prior_k)==1){
    if(is.element(prior_k, "default")){
      prior_k <- c(1.60, 1.25)
    }}

  # distant dataset fitting
  m <- mean(unrelated_snp_dist)
  v <- var(unrelated_snp_dist)
  size <- if (v > m) {
    m^2/(v - m)
  }else{100}

  nb_fit <- fitdistrplus::fitdist(unrelated_snp_dist, dist="truncnbinom", start=list(mu=m, size=size))

  lambda <- exp(seq(log(lambda_limits[1]),log(lambda_limits[2]), length.out = resolution))
  k=seq(k_limits[1],k_limits[2],length.out=resolution)

  # Generate grid of k and lambda combinations
  data <- expand.grid(k = k, lambda = lambda)

  # Compute likelihood for each k, lambda combination
  data$likelyhood <- future.apply::future_apply(data, 1, function(row) {
    k_val <- row["k"]
    lambda_val <- row["lambda"]

    # Compute likelihood using the negllk function
    negllk(k_val, lambda_val,
           x = trans_snp_dist,
           t = trans_time_dist,
           s = trans_sites,
           nbfitmu = nb_fit$estimate["mu"],
           nbfitsize = nb_fit$estimate["size"],
           truncation_point=truncation_point,
           prior_k=prior_k,
           prior_lambda=prior_lambda)
  })

  # Plot the results using ggplot
  ggplot(data, aes(y = k, x = lambda, z = likelyhood)) +
    scale_x_continuous(limits=lambda_limits, expand=c(0,0), transform = "log",
                       breaks = signif(10^(seq(log10(lambda_limits[1]),log10(lambda_limits[2]), length.out = 5)),1),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    scale_y_continuous(limits=k_limits, expand=c(0,0))+
    geom_contour_filled(bins=bins) +
    #geom_contour(colour="black", bins=bins)+
    geomtextpath::geom_textcontour(bins=bins, size = 2.5, padding = unit(0.05, "in"))+
    annotate("point",x=estimate[1], y=estimate[2],shape=4, color="red3")+
    annotate("rect", xmin = low_ci[1], ymin = low_ci[2],
             xmax = high_ci[1], ymax = high_ci[2],
             fill = NA, color = "red3", linetype = "dotted", size = 0.5)+
    labs(title = title,
         x="Rate",
         y="Proportion Related")+
    theme(legend.position="none")
}
