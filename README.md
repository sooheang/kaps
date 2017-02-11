<!-- README.md is generated from README.Rmd. Please edit that file -->
kaps
====

*K*-adaptive partitioning for survival data

[![Build Status](https://travis-ci.org/sooheang/kaps.svg?branch=master)](https://travis-ci.org/sooheang/kaps) [![Build status](https://ci.appveyor.com/api/projects/status/x7frit6xwy688d5e/branch/master?svg=true)](https://ci.appveyor.com/project/sooheang/kaps/branch/master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kaps)](http://cran.r-project.org/package=kaps) [![Coverage Status](https://img.shields.io/codecov/c/github/sooheang/kaps/master.svg)](https://codecov.io/github/sooheang/kaps?branch=master)

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

Tutorial (Under Construction)
=============================

References
==========

-   Kang, H. J.<sup>+</sup>, Eo, S-H<sup>+</sup>, Kim, S. C., Park, K. M., Lee, Y. J., Lee, S. K., Yu, E., Cho, H., and Hong, S-M (2014). Increased Number of Metastatic Lymph Nodes in Adenocarcinoma of the Ampulla of Vater as a Prognostic Factor: A Proposal of New Nodal Classification, *Surgery*, **155**(1), 74-84. (+:co-first author) [Link](http://dx.doi.org/10.1016/j.surg.2013.08.004)

-   Eo, S-H, Kang, H. J., Hong, S-M, and Cho, H. (2014). K-Adaptive Partitioning for Survival Data, with an Application to Cancer Staging, *arXiv*, **1306.4615**. [Link](https://arxiv.org/abs/1306.4615)

ChangeLog
=========

-   v1.1.0: Intitial release from GitHub
