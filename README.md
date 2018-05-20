GLASSOO
================

See [manual](https://github.com/MGallow/GLASSO/blob/master/GLASSOO.pdf).

Overview
--------

<br>

<p align="center">
<img src="lik.gif">
</p>
<br>

`GLASSO` is an R package that estimates a penalized precision matrix via block-wise coordinate descent -- also known as the graphical lasso (glasso) algorithm. A (possibly incomplete) list of functions contained in the package can be found below:

-   `GLASSO()` computes the estimated precision matrix

-   `plot.GLASSO()` produces a heat map or line graph for cross validation errors

Installation
------------

``` r
# The easiest way to install is from GitHub:
# install.packages("devtools")
devtools::install_github("MGallow/GLASSOO")
```

If there are any issues/bugs, please let me know: [github](https://github.com/MGallow/GLASSOO/issues). You can also contact me via my [website](http://users.stat.umn.edu/~gall0441/). Pull requests are welcome!

Usage
-----

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
(Omega = qr.solve(S) %>% round(3))
```

    ##        [,1]   [,2]   [,3]   [,4]   [,5]
    ## [1,]  1.961 -1.373  0.000  0.000  0.000
    ## [2,] -1.373  2.922 -1.373  0.000  0.000
    ## [3,]  0.000 -1.373  2.922 -1.373  0.000
    ## [4,]  0.000  0.000 -1.373  2.922 -1.373
    ## [5,]  0.000  0.000  0.000 -1.373  1.961

``` r
# generate 1000 x 5 matrix with rows drawn from iid N_p(0, S)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt

# calculate sample covariance
Sample = (nrow(X) - 1)/nrow(X)*cov(X)

# print sample precision matrix (perhaps a bad estimate)
(qr.solve(cov(X)) %>% round(5))
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.98365 -1.54113 -0.03896 -0.04392  0.11939
    ## [2,] -1.54113  3.13590 -1.14656  0.07428 -0.27851
    ## [3,] -0.03896 -1.14656  2.44800 -1.08953  0.00029
    ## [4,] -0.04392  0.07428 -1.08953  2.40163 -1.10604
    ## [5,]  0.11939 -0.27851  0.00029 -1.10604  1.87394

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
    ## Log-likelihood: -11.94762
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  0.99244 -0.23296  0.00000  0.00000  0.00000
    ## [2,] -0.23296  1.16354 -0.14944  0.00000  0.00000
    ## [3,]  0.00000 -0.14944  1.04543 -0.15200  0.00000
    ## [4,]  0.00000  0.00000 -0.15200  1.05441 -0.13787
    ## [5,]  0.00000  0.00000  0.00000 -0.13787  1.08982

``` r
# GLASSO cross validation
GLASSO(X, lam = 10^seq(-5, 5, 0.5))
```

    ## 
    ## Call: GLASSO(X = X, lam = 10^seq(-5, 5, 0.5))
    ## 
    ## Iterations:
    ## [1] 7
    ## 
    ## Tuning parameter:
    ##       log10(lam)    lam
    ## [1,]        -1.5  0.032
    ## 
    ## Log-likelihood: -137.19636
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.82323 -1.31666 -0.06110  0.00000  0.00000
    ## [2,] -1.31666  2.76907 -0.99074  0.00000 -0.15513
    ## [3,] -0.06110 -0.99074  2.22827 -0.94831 -0.03176
    ## [4,]  0.00000  0.00000 -0.94831  2.19159 -0.97400
    ## [5,]  0.00000 -0.15513 -0.03176 -0.97400  1.75654

``` r
# produce CV heat map for GLASSO
GLASSO = GLASSO(X, lam = 10^seq(-5, 5, 0.01))
GLASSO %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
# produce line graph for CV errors for GLASSO
GLASSO %>% plot(type = "line")
```

![](README_files/figure-markdown_github/unnamed-chunk-2-2.png)
