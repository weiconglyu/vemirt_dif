
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

The goal of `vemirt` is to detect differential item functioning (DIF) in
two-parameter logistic (2PL) models via (Gaussian variational)
expectation-maximization estimation.

## Installation

You can install the development version of `vemirt` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("weiconglyu/vemirt")
```

## Example

This is a basic example which shows you how to apply Gaussian
variational expectation-maximization estimation with importance sampling
to the simulated data from the package:

``` r
library(vemirt)
result <- with(DIF_simdata, DIF_GVEMM(Y, D, X, 'IWGVEMM', Lambda0 = c(0.3, 0.5)))
#> Running GVEMM for initial values...
#> Fitting the model using different lambdas...
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
round(result[[1]]$beta, 3)
#>       [,1]  [,2]  [,3]  [,4] [,5]  [,6]   [,7]   [,8]   [,9]  [,10]  [,11]
#> [1,] 0.000 0.000 0.000 0.000 0.00 0.000  0.000  0.000  0.000  0.000  0.000
#> [2,] 0.370 0.000 0.345 0.423 0.00 0.000 -0.305 -0.613 -0.429  0.000 -0.563
#> [3,] 0.892 0.362 0.927 1.211 0.96 0.654 -0.906 -1.105 -0.932 -0.681 -1.061
#>       [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
#> [1,]  0.000     0     0     0     0     0     0 0.000     0
#> [2,] -0.512     0     0     0     0     0     0 0.192     0
#> [3,] -0.909     0     0     0     0     0     0 0.000     0
round(result[[2]]$beta, 3)
#>      [,1] [,2] [,3]  [,4]  [,5]  [,6]  [,7]   [,8]   [,9]  [,10] [,11]  [,12]
#> [1,] 0.00    0 0.00 0.000 0.000 0.000  0.00  0.000  0.000  0.000  0.00  0.000
#> [2,] 0.00    0 0.00 0.000 0.000 0.000  0.00 -0.612  0.000  0.000  0.00  0.000
#> [3,] 0.69    0 0.74 0.938 0.956 0.609 -0.89 -1.128 -0.737 -0.707 -1.05 -0.685
#>      [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
#> [1,]     0     0     0     0     0     0     0     0
#> [2,]     0     0     0     0     0     0     0     0
#> [3,]     0     0     0     0     0     0     0     0
```
