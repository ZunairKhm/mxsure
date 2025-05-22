#' Contour plot for likelyhood estimation of mixture distribution mutation rate estimation
#'
#' @param mixed_snp_dist list of SNP distances from a mixed transmission data set
#' @param unrelated_snp_dist list of SNP distances from an unrelated data set
#' @param mixed_time_dist list of time differences between samples from each SNP distance in the mixed data set (in days)
#' @param mixed_sites list of sites considered for each SNP distance in mixed data set
#' @param resolution resolution of plot, how many points along k and lambda to calculate (main computation limitation)
#' @param lambda_limits bounds of lambda axis (plots on a log10 scale)
#' @param k_limits bounds of k axis
#' @param bins number of bins to pass to ggplot contour function
#' @param title title to pass to ggplot title
#' @param truncation_point a SNP distance limit for the data, if set to NA will estimate as if there is no limit
#' @param prior_lambda parameters for a gamma prior distribution for rate estimation, if set to "default" will use some default parameters
#' @param prior_k parameters for a beta prior distribution for related proportion estimation, if set to "default" will use some default parameters
#' @param lambda_bounds bounds of rate estimation in SNPs/year/site if given time and site data
#' @param k_bounds bounds of related proportion estimation
#' @param intercept_bounds bounds of intercept estimation
#' @param low_ci 2 element vector for the low CIs for lambda and k respectively. If supplied, alongside high_ci, will produce a box around the estimate
#' @param high_ci 2 element vector for the high CIs for lambda and k respectively. If supplied, alongside low_ci, will produce a box around the estimate
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous geom_contour_filled unit annotate labs theme
#' @importFrom tibble tibble
#'
#' @return ggplot contour map of likelyhood space
#' @export
mxsure_likelyhood_contour <- function(mixed_snp_dist, unrelated_snp_dist, mixed_time_dist, mixed_sites, truncation_point=NA,
                                       prior_lambda=NA, prior_k=NA, lambda_bounds=c(0, 1), k_bounds=c(0,1), intercept_bounds=c(-Inf, Inf),
                                       resolution=100, lambda_limits=c(1e-8, 1e-4), k_limits=c(0.2,0.99), low_ci=NULL, high_ci=NULL,
                                       bins=NULL, title="Likelihood Contour Plot"){

  likelyhood<- ".x" <- NULL

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

  negllk <- function(k, lambda, x, t, s, intercept, nbfitmu, nbfitsize,truncation_point, prior_k, prior_lambda){

    -log(-sum(pmap_dbl(list(x, t,s), ~ {log_sum_exp(log(k) + dpois(x = ..1,
                                                                   lambda =  lambda*..2*(..3/1000000)+intercept, #gives rate esimate per day time per million bp
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


  lambda_limits <- lambda_limits*1e6/365.25

  result <- suppressWarnings(
    mxsure_estimate(
      mixed_snp_dist,unrelated_snp_dist, mixed_time_dist, mixed_sites, truncation_point=truncation_point, start_params = NA,
      lambda_bounds = lambda_bounds, k_bounds=k_bounds,  intercept_bounds=intercept_bounds)
    , classes = "warning")


  #truncating data
  if(is.na(truncation_point)){
    truncation_point <- Inf
  }
  unrelated_snp_dist_orig <- unrelated_snp_dist
  x <- tibble(mixed_snp_dist, mixed_time_dist, mixed_sites)
  x <- filter(x, mixed_snp_dist<truncation_point)
  mixed_snp_dist <- x$mixed_snp_dist
  mixed_time_dist <- x$mixed_time_dist
  mixed_sites <- x$mixed_sites
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
           x = mixed_snp_dist,
           t = mixed_time_dist,
           s = mixed_sites,
           intercept = result$intercept,
           nbfitmu = result$nb_mu,
           nbfitsize = result$nb_size,
           truncation_point=truncation_point,
           prior_k=prior_k,
           prior_lambda=prior_lambda)
  })

  # Plot the results using ggplot
  ggplot(data, aes(y = k, x = lambda*365.25/1e6, z = likelyhood)) +
    scale_x_continuous(limits=lambda_limits*365.25/1e6, expand=c(0,0), transform = "log",
                       breaks = signif(10^(seq(log10(lambda_limits[1]*365.25/1e6),log10(lambda_limits[2]*365.25/1e6), length.out = 5)),1),
                       labels = scales::math_format(10^.x))+
    scale_y_continuous(limits=k_limits, expand=c(0,0))+
    geom_contour_filled(bins=bins) +
    #geom_contour(colour="black", bins=bins)+
    geomtextpath::geom_textcontour(bins=bins, size = 2.5, padding = unit(0.05, "in"))+
    annotate("point",x=result$lambda, y=result$k,shape=4, color="red3", size=3)+
    annotate("rect", xmin = low_ci[1], ymin = low_ci[2],
             xmax = high_ci[1], ymax = high_ci[2],
             fill = NA, color = "red3", linetype = "dotted", size = 0.5)+
    labs(title = title,
         x="Rate (SNPs/year/site)",
         y="Proportion Related")+
    theme(legend.position="none")
}

