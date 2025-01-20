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

negllk <- function(k, lambda, x, t, s, nbfitmu, nbfitsize){

  -log(-sum(pmap_dbl(list(x, t,s), ~ {log_sum_exp(log(k) + dpois(x = ..1,
                                                           lambda =  lambda*..2*(..3/1000000), #gives rate esimate per day time per million bp
                                                           log = TRUE),
                                            log(1-k) + dnbinom(x = ..1,
                                                               size = nbfitsize,
                                                               mu = nbfitmu,
                                                               log = TRUE))})))
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
mixture_likelyhood_contour <- function(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites,
                                       resolution=100, lambda_limits=NA, k_limits=NA, bins=NULL, title="Mixture Likelihood Contour Plot"){
  if(anyNA(lambda_limits)){
    temp_res <- mixture_snp_cutoff(trans_snp_dist, unrelated_snp_dist, trans_time_dist, trans_sites)
    lambda_limits <- c(temp_res$lambda*0.1,temp_res$lambda*10)
    k_limits <- c(max(c (temp_res$k-0.1, 0)),min(c(temp_res$k+0.1,1)))
  }
  nb_fit <- suppressWarnings( MASS::fitdistr(x=unrelated_snp_dist, densfun = "negative binomial"
  ))
  lambda <- exp(seq(log(lambda_limits[1]),log(lambda_limits[2]), length.out = resolution))
  k=seq(k_limits[1],k_limits[2],length.out=resolution)

  # Generate grid of k and lambda combinations
  data <- expand.grid(k = k, lambda = lambda)

  # Compute likelihood for each k, lambda combination
  data$likelyhood <- future.apply::future_apply(data, 1, function(row) {
    k_val <- row["k"]
    lambda_val <- row["lambda"]

    # Compute likelihood using the negllk function (assuming it is defined elsewhere)
    negllk(k_val, lambda_val,
           x = trans_snp_dist,
           t = trans_time_dist,
           s = trans_sites,
           nbfitmu = nb_fit$estimate["mu"],
           nbfitsize = nb_fit$estimate["size"])
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
    labs(title = title,
         x="Rate",
         y="Proportion Related")+
    theme(legend.position="none")
}
