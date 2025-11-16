
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MxSure

<!-- badges: start -->

[![R-CMD-check](https://github.com/zunairkhm/mxsure/workflows/R-CMD-check/badge.svg)](https://github.com/zunairkhm/mxsure/actions)

<!-- [![DOI](https://zenodo.org/badge/XXXX.svg)](https://zenodo.org/badge/latestdoi/XXXX) -->

<!-- badges: end -->

This package implements a mixture distribution approach to estimating
substitution rates and SNP thresholds from pairwise SNP comparison data,
in a longitudinal/within-host or transmission setting. It also provides
some plotting and simulation functions.

## Installation

You can install the development version of mxsure from
[GitHub](https://github.com/) with:

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("ZunairKhm/mxsure", quiet=T)
library(mxsure)
```

## Loading in data

MxSure can input data from inStrain compare output, TRACS distance
output, or from a distance matrix. Either provide an R object or
filepath to the respective function and the output will be an R object
ready for inputting into MxSure. These examples only show loading in one
dataset for each method however a mixed and distant dataset is required
for MxSure.

### TRACS distance output

``` r
#loading example data provided with package
tracs_example_path <- system.file("extdata", "example_tracs_output.csv", package = "mxsure")
sample_dates <- readr::read_csv(system.file("extdata", "example_sampling_dates.csv", package = "mxsure"))
#> Rows: 30 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr  (1): sample_id
#> date (1): date
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
head(sample_dates)
#> # A tibble: 6 × 2
#>   sample_id date      
#>   <chr>     <date>    
#> 1 sample1   2025-01-01
#> 2 sample2   2025-01-07
#> 3 sample3   2025-01-13
#> 4 sample4   2025-01-19
#> 5 sample5   2025-01-25
#> 6 sample6   2025-01-31

#using function to ready data for input into mxsure
mxsure_input_tracs <-  mxsure_input_tracs(file_path=tracs_example_path, dates=sample_dates, dates_sample_col = "sample_id", dates_date_col = "date")
head(mxsure_input_tracs)
#>   snp_dist time_dist   sites
#> 1     1653       174  804627
#> 2     2001         6 1660022
#> 3     1557         6 1453158
#> 4     1805         6 2079670
#> 5     1598         6  858564
#> 6     2125         6 1039812
```

### inStrain compare output

``` r
instrain_example_path <- system.file("extdata", "example_instrain_output.csv", package = "mxsure")
sample_dates <- readr::read_csv(system.file("extdata", "example_sampling_dates.csv", package = "mxsure"))
#> Rows: 30 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr  (1): sample_id
#> date (1): date
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
mxsure_input_instrain <-  mxsure_input_instrain(file_path=instrain_example_path, dates=sample_dates, dates_sample_col = "sample_id", dates_date_col = "date")
head(mxsure_input_instrain)
#>   snp_dist time_dist   sites
#> 1        2       174 4074096
#> 2     8851         6 4152305
#> 3     6832         6 4104415
#> 4     6029         6 4059303
#> 5     7111         6 4113956
#> 6     6171         6 4138456
```

### distance matrix

note: this function does not allow for ‘sites considered’ or equivalent;
this will need to be inputed manually

``` r
distmat_example_path <- system.file("extdata", "example_distance_matrix.csv", package = "mxsure")
sample_dates <- readr::read_csv(system.file("extdata", "example_sampling_dates.csv", package = "mxsure"))
#> Rows: 30 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr  (1): sample_id
#> date (1): date
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
mxsure_input_distmat <-  mxsure_input_distmatrix(file_path=distmat_example_path, dates=sample_dates, dates_sample_col = "sample_id", dates_date_col = "date")
head(mxsure_input_distmat)
#>    sampleA sampleB snp_dist time_dist
#>     <char>  <char>    <int>     <num>
#> 1: sample2 sample1      467         6
#> 2: sample3 sample1      183        12
#> 3: sample4 sample1       18        18
#> 4: sample5 sample1      199        24
#> 5: sample6 sample1      430        30
#> 6: sample7 sample1      310        36
```

## Basic Operation

MxSure in an R package designed to estimate substitution rates and infer
SNP thresholds from longitudinally collected pairwise SNP comparison
data. This data is likely to contain a mixture of related pairs (from
longitudinal samples of the same microbe showing evolution over time)
and unrelated pairs (co-captured samples from the same strain/species).

To demonstrate its use we can first simulate some example data using the
`simulate_mixsnp_data()` command. Here we simulate 100 data points (each
representing a pairwise SNP comparison between two sequences) with a
true substitution rate of 5 SNPs/year with each data point having an 80%
chance to be a related data point. By default sampling time differences
for each of these data points are between 0 and 1 year.

``` r
library(mxsure)
set.seed(123)
mixed_data <- simulate_mixsnp_data(lambda=5, k=0.8, n=100)
head(mixed_data, 5)
#> # A tibble: 5 × 3
#>   snp_dist time_dist relation 
#>      <dbl>     <dbl> <chr>    
#> 1      181     105.  Unrelated
#> 2        0      16.6 Related  
#> 3        8     349.  Related  
#> 4       34      89.9 Unrelated
#> 5        3     234.  Related
```

MxSure’s methodology is to fit a mixture distribution, which describes
the related and unrelated SNP distances in separate models, to this
dataset. However, to do this we need to fit the unrelated SNP distance
model (a right truncated negative binomial distribution) on another
dataset of pairwise SNP comparisons that has no related pairs. This
could be comparisons between different individuals and ideally
individuals that are geographically separated to exclude transmission
events. Here we again simulate the data using the same true substitution
rate but with related proportion, k, set to 0.

``` r
distant_data <- simulate_mixsnp_data(lambda=5, k=0, n=1000)
head(distant_data, 5)
#> # A tibble: 5 × 3
#>   snp_dist time_dist relation 
#>      <dbl>     <dbl> <chr>    
#> 1       36      5.90 Unrelated
#> 2       96    239.   Unrelated
#> 3      129    331.   Unrelated
#> 4      163    151.   Unrelated
#> 5      137    103.   Unrelated
```

We can then use `mxsure_estimate()` to estimate the substitution rate,
estimate the related proportion, and infer a SNP threshold.

``` r
result <- mxsure_estimate(
  mixed_snp_dist = mixed_data$snp_dist,
  unrelated_snp_dist = distant_data$snp_dist,
  mixed_time_dist = mixed_data$time_dist
)
result
#> # A tibble: 1 × 8
#>   snp_threshold lambda     k intercept estimated_fp lambda_units   nb_size nb_mu
#>           <dbl>  <dbl> <dbl>     <dbl>        <dbl> <chr>            <dbl> <dbl>
#> 1            10   5.71 0.770    -0.146        0.002 SNPs per year…    4.07  125.
```

Here lambda is the estimated substitution rate and k is the estimated
related proportion. Intercept is a time independent parameter that
describes the expected SNPs at time 0, which can be used to model
sequencing error. Estimated fp displayed the proportion of the distant
dataset that has a SNP distance below the inferred SNP threshold.
nb_size and nb_mu are the parameters of the unrelated SNP distance
model.

As in this example, if sites considered are not supplied for each
comparison in the mixed dataset the estimated rate will have units
SNPs/year/genome. We can display what this by simulating some amount of
sites.

``` r
mixed_data_sites <- mixed_data
mixed_data_sites$sites <- runif(100, 0.5e6, 1.5e6)

result_sites <- mxsure_estimate(
  mixed_snp_dist = mixed_data_sites$snp_dist,
  unrelated_snp_dist = distant_data$snp_dist,
  mixed_time_dist = mixed_data_sites$time_dist,
  mixed_sites = mixed_data_sites$sites
)
result_sites
#> # A tibble: 1 × 12
#>   snp_threshold     lambda     k intercept estimated_fp lambda_units convergence
#>           <dbl>      <dbl> <dbl>     <dbl>        <dbl> <chr>              <int>
#> 1            10 0.00000565 0.770   -0.0324        0.002 SNPs per ye…           0
#> # ℹ 5 more variables: message <chr>, iterations <int>, nb_size <dbl>,
#> #   nb_mu <dbl>, lambda_per_genome <dbl>
```

## Youden and Threshold Range

For threshold assessment it may be useful to compare the MxSure
threshold with what would be inferred from the Youden method.

``` r
result <- mxsure_estimate(
  mixed_snp_dist = mixed_data$snp_dist,
  unrelated_snp_dist = distant_data$snp_dist,
  mixed_time_dist = mixed_data$time_dist,
  youden = TRUE
)
result$youden
#> # A tibble: 1 × 3
#>   youden_snp_threshold     J youden_estimated_fp
#>                  <dbl> <dbl>               <dbl>
#> 1                   11 0.768               0.002
```

It may also be useful to examine how a different threshold time affects
the MxSure threshold.

``` r
result <- mxsure_estimate(
  mixed_snp_dist = mixed_data$snp_dist,
  unrelated_snp_dist = distant_data$snp_dist,
  mixed_time_dist = mixed_data$time_dist,
  threshold_range = TRUE
)
result$threshold_range
#>    years threshold estimated_fp prop_pos
#> 1    0.5         6        0.000     0.69
#> 2    1.0        10        0.002     0.76
#> 3    1.5        13        0.002     0.77
#> 4    2.0        17        0.003     0.77
#> 5    2.5        21        0.010     0.77
#> 6    3.0        24        0.013     0.77
#> 7    3.5        27        0.019     0.77
#> 8    4.0        31        0.026     0.77
#> 9    4.5        34        0.036     0.78
#> 10   5.0        37        0.043     0.78
#> 11   5.5        41        0.058     0.79
#> 12   6.0        44        0.066     0.79
#> 13   6.5        47        0.077     0.79
#> 14   7.0        50        0.090     0.79
#> 15   7.5        54        0.103     0.79
#> 16   8.0        57        0.121     0.80
#> 17   8.5        60        0.139     0.81
#> 18   9.0        63        0.159     0.81
#> 19   9.5        66        0.178     0.82
#> 20  10.0        70        0.200     0.82
```

## Confidence Intervals

Confidence intervals for substitution rate estimation are produced via
bootstrapping. This can be computed faster with multi-threading using
the future package. 20 bootstraps is very likely too few, you would
likely need at least 200 but ideally 500 or 1000.

``` r
future::plan("future::multisession", workers=2)
ci <- mxsure_ci( mixed_snp_dist = mixed_data$snp_dist,
  unrelated_snp_dist = distant_data$snp_dist,
  mixed_time_dist = mixed_data$time_dist,
  bootstraps=20
  )
future::plan("future::sequential")
ci$confidence_intervals
#> # A tibble: 2 × 5
#>   snp_threshold lambda     k intercept estimated_fp
#>           <dbl>  <dbl> <dbl>     <dbl>        <dbl>
#> 1             9   5.21 0.760    -0.159       0     
#> 2            10   6.07 0.820    -0.130       0.0041
```

## Time Randomisation Test

This is a test that validates there is temporal signal in the underlying
data. It achieves this by randomly permuting the time column in the
mixed dataset and produces estimates and 95% CIs on this time randomised
dataset. This is repeated on multiple time randomisation datasets, 10 is
recommended, and if more than 5% of these bootstrap samples are above or
equal to the original estimate this test is failed. Estimates should not
be considered to be valid if this test fails. The output of this
function is a graph showing this result and an outcome table.

``` r
future::plan("future::multisession", workers=2)
timerand <- mxsure_timerandtest( mixed_snp_dist = mixed_data$snp_dist,
  unrelated_snp_dist = distant_data$snp_dist,
  mixed_time_dist = mixed_data$time_dist,
  bootstraps=20,
  permutations=2,
  original_result = result$results,
  ci_data = ci
  )
#> [1] "Processing Permutation: 1"
#> [1] "Processing Permutation: 2"
future::plan("future::sequential")
timerand$plot
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

``` r
timerand$outcome
#> # A tibble: 1 × 8
#>   n_permutations any_overlapping_est n_overlapping_est perc_overlapping_est
#>            <dbl> <lgl>                           <int>                <dbl>
#> 1              2 FALSE                               0                    0
#> # ℹ 4 more variables: any_overlapping_lowci <lgl>, n_overlapping_lowci <int>,
#> #   perc_overlapping_lowci <dbl>, failure_perc <dbl>
```

## Likelihood ratios, and SNP over Time plot

For each data point in the mixed dataset, after fitting, we can assess
the likelihood that any is related or unrelated. This is achieved using
the `mxsure_likelihood()` function which reports likelihoods, log
likelihoods, and likelihood ratios.

``` r
likelihoods <- mxsure_likelihood(mixed_snp_dist = mixed_data$snp_dist,
  unrelated_snp_dist = distant_data$snp_dist,
  mixed_time_dist = mixed_data$time_dist)
head(likelihoods, 5)
#> # A tibble: 5 × 9
#>   snp_dist time_dist sites rel_loglh unrel_loglh  logLHR    rel_lh    unrel_lh
#>      <dbl>     <dbl> <dbl>     <dbl>       <dbl>   <dbl>     <dbl>       <dbl>
#> 1      181     105.      1  -692.          -5.76 -686.   3.00e-301 0.00315    
#> 2        0      16.6     1    -0.114      -14.1    13.9  8.92e-  1 0.000000788
#> 3        8     349.      1    -2.55        -9.12    6.57 7.78e-  2 0.000109   
#> 4       34      89.9     1   -82.0         -6.02  -76.0  2.46e- 36 0.00242    
#> 5        3     234.      1    -1.54       -11.1     9.58 2.15e-  1 0.0000149  
#> # ℹ 1 more variable: LHR <dbl>
```

Based on these likelihoods we can plot the original data with SNP
distance versus time distance for each comparisons, using the
`snp_over_time()` function. We colour the points as their likelihoods.
We can also plot the predictive intervals and inferred SNP threshold.
`snp_over_time()` calls `mxsure_likelihood()` under the hood.

``` r
snp_over_time(mixed_snp_dist = mixed_data$snp_dist,
  unrelated_snp_dist = distant_data$snp_dist,
  mixed_time_dist = mixed_data$time_dist
  )
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_step()`).
#> Removed 9 rows containing missing values or values outside the scale range
#> (`geom_step()`).
#> Removed 9 rows containing missing values or values outside the scale range
#> (`geom_step()`).
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />
