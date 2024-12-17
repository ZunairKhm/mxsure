
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixturesnpcutoff

<!-- badges: start -->
<!-- badges: end -->

This package implements a mixture distribution approach to estimating
mutation rates and SNP thresholds for paired metagenomic data. It also
provides some plotting functions and a simulation functions.

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
library(mixturesnpcutoff)
#> Loading required package: tidyverse
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
## basic example code
x <- simulate_mixsnp_data(0.00179*4, 0.8, 500)
y <- simulate_unrelsnp_data()

z<- mixture_snp_cutoff(x$snp_dist, y, x$time_dist)
#> Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
#> Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

mixture_snp_cutoff will produce a mutation rate in SNPs/day or
SNPs/day/million basepairs if sites considered is provided. It will also
produce a cutoff that considers a max time (default 2 years) for SNP
distances that can be considered related/transmitted or unrelated.
