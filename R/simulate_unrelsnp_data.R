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



