GLASSOO
================

[![Build
Status](https://travis-ci.org/MGallow/GLASSOO.svg?branch=master)](https://travis-ci.org/MGallow/GLASSOO)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/GLASSOO)](https://cran.r-project.org/package=GLASSOO)

## Overview

<center>

![](https://github.com/MGallow/GLASSOO/raw/master/vignettes/images/gif.gif)

</center>

`GLASSOO` is an R package that estimates a lasso-penalized precision
matrix via block-wise coordinate descent – also known as the graphical
lasso (glasso) algorithm. This package is similar to
[CVglasso](https://mgallow.github.io/CVglasso/) – but rather than being
a wrapper around the
[glasso](https://cran.r-project.org/web/packages/glasso/index.html)
package, the code is completely re-written in C++. A (possibly
incomplete) list of functions contained in the package can be found
below:

  - `GLASSO()` computes the estimated precision matrix

  - `plot.GLASSO()` produces a heat map or line graph for cross
    validation errors

See [vignette](https://mgallow.github.io/GLASSOO/) or
[manual](https://github.com/MGallow/GLASSOO/blob/master/GLASSOO.pdf).

## Installation

``` r
# The easiest way to install is from GitHub:
# install.packages("devtools")
devtools::install_github("MGallow/GLASSOO")
```

If there are any issues/bugs, please let me know:
[github](https://github.com/MGallow/GLASSOO/issues). You can also
contact me via my [website](https://mgallow.github.io/). Pull requests
are welcome\!

## Usage

``` r
library(GLASSOO)

# generate data from a sparse matrix
# first compute covariance matrix
S = matrix(0.7, nrow = 5, ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    S[i, j] = S[i, j]^abs(i - j)
  }
}

# print oracle precision matrix (shrinkage might be useful)
(Omega = round(qr.solve(S), 3))
```

    ##        [,1]   [,2]   [,3]   [,4]   [,5]
    ## [1,]  1.961 -1.373  0.000  0.000  0.000
    ## [2,] -1.373  2.922 -1.373  0.000  0.000
    ## [3,]  0.000 -1.373  2.922 -1.373  0.000
    ## [4,]  0.000  0.000 -1.373  2.922 -1.373
    ## [5,]  0.000  0.000  0.000 -1.373  1.961

``` r
# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
set.seed(123)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt

# calculate sample covariance
Sample = (nrow(X) - 1)/nrow(X)*cov(X)

# print sample precision matrix (perhaps a bad estimate)
round(qr.solve(cov(X)), 5)
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.30646 -1.53483  0.21884 -0.08521  0.24066
    ## [2,] -1.53483  3.24286 -1.66346 -0.14134  0.18760
    ## [3,]  0.21884 -1.66346  3.16698 -1.23906 -0.10906
    ## [4,] -0.08521 -0.14134 -1.23906  2.74022 -1.35853
    ## [5,]  0.24066  0.18760 -0.10906 -1.35853  2.03323

``` r
# GLASSO (lam = 0.5)
GLASSO(S = Sample, lam = 0.5)
```

    ## 
    ## Call: GLASSO(S = Sample, lam = 0.5)
    ## 
    ## Iterations:
    ## [1] 3
    ## 
    ## Tuning parameter:
    ##       log10(lam)  lam
    ## [1,]      -0.301  0.5
    ## 
    ## Log-likelihood: -10.44936
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.34080 -0.00973  0.00000  0.00000  0.00000
    ## [2,] -0.00973  1.19263 -0.09615  0.00000  0.00000
    ## [3,]  0.00000 -0.09615  1.21895 -0.11424  0.00000
    ## [4,]  0.00000  0.00000 -0.11424  1.06968 -0.13534
    ## [5,]  0.00000  0.00000  0.00000 -0.13534  1.12473

``` r
# GLASSO cross validation
(GLASSO = GLASSO(X))
```

    ## 
    ## Call: GLASSO(X = X)
    ## 
    ## Iterations:
    ## [1] 3
    ## 
    ## Tuning parameter:
    ##       log10(lam)    lam
    ## [1,]      -1.544  0.029
    ## 
    ## Log-likelihood: -110.16675
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.13226 -1.24669  0.00000  0.00000  0.18709
    ## [2,] -1.24669  2.75122 -1.29912 -0.07341  0.00000
    ## [3,]  0.00000 -1.29912  2.81744 -1.15682 -0.00113
    ## [4,]  0.00000 -0.07341 -1.15682  2.46461 -1.17086
    ## [5,]  0.18709  0.00000 -0.00113 -1.17086  1.86326

``` r
# produce line graph for CV errors for GLASSO
plot(GLASSO)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# produce CV heat map for GLASSO
plot(GLASSO, type = "heatmap")
```

![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->
