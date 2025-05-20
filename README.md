
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MxSure

<!-- badges: start -->

<!-- badges: end -->

This package implements a mixture distribution approach to estimating
mutation rates and SNP thresholds for paired metagenomic data. It also
provides some plotting functions and simulation functions.

## Installation

You can install the development version of mixturesnpcutoff from
[GitHub](https://github.com/) with:

``` r
install.packages("pak")
pak::pak("ZunairKhm/mixturesnpcutoff")
```

## Example

This shows the basic use case using data produced from the simulation
functions.

``` r
library(mxsure)
## basic example code
x <- simulate_mixsnp_data(1, 0.8)
y <- simulate_mixsnp_data(1, 0, n=1000)

mxsure_estimate(x$snp_dist, y$snp_dist, x$time_dist)
#> # A tibble: 1 × 8
#>   snp_threshold lambda     k intercept estimated_fp lambda_units   nb_size nb_mu
#>           <dbl>  <dbl> <dbl>     <dbl>        <dbl> <chr>            <dbl> <dbl>
#> 1             3   1.02 0.809  -0.00150        0.009 SNPs per year…    4.35  25.2
```

mxsure_estimate will produce a mutation rate in SNPs/year or
SNPs/year/site if sites considered is provided. It will also produce a
cutoff that considers a max time for SNP distances that can be
considered related/transmitted or unrelated.
