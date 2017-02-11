<!-- README.md is generated from README.Rmd. Please edit that file -->
kaps
====

*K*-adaptive partitioning for survival data

[![Build Status](https://travis-ci.org/sooheang/kaps.svg?branch=master)](https://travis-ci.org/sooheang/kaps) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/sooheang/sooheang?branch=master&svg=true)](https://ci.appveyor.com/project/sooheang/kaps) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kaps)](http://cran.r-project.org/package=kaps) [![Coverage Status](https://img.shields.io/codecov/c/github/sooheang/kaps/master.svg)](https://codecov.io/github/sooheang/kaps?branch=master)

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

CHANELOG
========

-   v1.1.0: Intitial release

References
==========
