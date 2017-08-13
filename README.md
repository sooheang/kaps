<!-- README.md is generated from README.Rmd. Please edit that file -->
kaps
====

*K*-adaptive partitioning for survival data

[![Build Status](https://travis-ci.org/sooheang/kaps.svg?branch=master)](https://travis-ci.org/sooheang/kaps) [![Build status](https://ci.appveyor.com/api/projects/status/x7frit6xwy688d5e/branch/master?svg=true)](https://ci.appveyor.com/project/sooheang/kaps/branch/master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kaps)](http://cran.r-project.org/package=kaps) [![Coverage Status](https://img.shields.io/codecov/c/github/sooheang/kaps/master.svg)](https://codecov.io/github/sooheang/kaps?branch=master) [![CRAN\_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/kaps)](http://cran.rstudio.com/web/packages/kaps/index.html)

kaps is a multi-way partitioning algorithm, which divides the data into *K* heterogeneous subgroups based on the information from a prognostic factor. Such a multi-way partition is found by maximizing the minimum of the subgroup pairwise test statistics. An optimal number of subgroups is determined by a permutation test.

You can install:

-   the latest released version from CRAN with

    ``` r
    install.packages("kaps")
    ```

-   the latest development version from github with

    ``` r
    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("sooheang/kaps")
    ```

If you encounter a clear bug, please file a minimal reproducible example on [github](https://github.com/sooheang/kaps/issues).

Tutorial
--------

To illustrate the function kaps() with various options, we use an artificial data, toy, which consists of 150 artificial observations of the survival time (time), its censoring status (status) and 6 covariates: the patient's race (Race), age at the initial time (Age), pathological grade (Grade), early onset depression (EOD), the number of metastasis LNs (meta) and the number of examined LNs (exam). The data can be called up from the package `kaps`:

``` r
library(kaps)
#> Loading required package: survival
data('toy', package = 'kaps')
head(toy)
#>   meta status time
#> 1    1      0    0
#> 2    4      1   26
#> 3    0      1   22
#> 4    9      1   15
#> 5    0      1   70
#> 6    1      0   96
```

Here we utilize just 3 variables: meta, status, and time. The number of metastasis LNs, meta, is used as an ordered prognostic factor for finding heterogeneous subgroups.

``` r
toy <- toy[, c("meta", "status", "time")]
```

### Selecting a set of cut-off points for given *K*

Suppose we specify the number of subgroups in advance. For instance, *K* = 3. To select an optimal set of two cut-off points when *K* = 3, the function kaps is called via the following statements

``` r
fit1 <- kaps(survival::Surv(time, status) ~ meta, data = toy, K = 3)
#> Now, finding optimal number of subgroups (K) by KAPS-permutation. 
#> Now, selecting a set of cut-off points...
fit1 
#> Call:
#> kaps(formula = survival::Surv(time, status) ~ meta, data = toy, 
#>     K = 3)
#> 
#>  K-Adaptive Partitioning for Survival Data
#> 
#> Samples= 150                 Optimal K=3 
#> 
#> 
#> Selecting a set of cut-off points:
#>       Xk df Pr(>|Xk|)  X1 df Pr(>|X1|) adj.Pr(|X1|) cut-off points  
#> K=3 36.8  2         0 7.2  1    0.0073     0.014001          0, 10 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> P-values of pairwise comparisons
#>             0<=meta<=0 0<meta<=10
#> 0<meta<=10       1e-04          -
#> 10<meta<=38     <.0000     0.0073
```

References
----------

-   Kang, H. J.<sup>+</sup>, Eo, S-H<sup>+</sup>, Kim, S. C., Park, K. M., Lee, Y. J., Lee, S. K., Yu, E., Cho, H., and Hong, S-M (2014). Increased Number of Metastatic Lymph Nodes in Adenocarcinoma of the Ampulla of Vater as a Prognostic Factor: A Proposal of New Nodal Classification, *Surgery*, **155**(1), 74-84. (+:co-first author) [Link](http://dx.doi.org/10.1016/j.surg.2013.08.004)

-   Eo, S-H, Kang, H. J., Hong, S-M, and Cho, H. (2014). K-Adaptive Partitioning for Survival Data, with an Application to Cancer Staging, *arXiv*, **1306.4615**. [Link](https://arxiv.org/abs/1306.4615)

ChangeLog
=========

-   v1.1.3: Update survival function
-   v1.1.0: Intitial release on GitHub
