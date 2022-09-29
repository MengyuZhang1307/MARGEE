# MARGEE (Mixed-model Association robust tests for GEne–Environment interactions)

## Description

MARGEE is an R package for gene-environment interaction (GEI) tests and joint tests for GWAS. Perform generalized linear mixed robust single-variant main effect tests, gene-environment interaction tests.

## Installing

MARGEE imports R packages 
[Rcpp](https://cran.r-project.org/web/packages/CompQuadForm/index.html), 
[Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), 
[parallel](https://cran.r-project.org/web/views/HighPerformanceComputing.html), 
[MASS](https://cran.r-project.org/web/packages/MASS/index.html), 
[SeqArray](http://bioconductor.org/packages/release/bioc/html/SeqArray.html), 
[SeqVarTools](https://bioconductor.org/packages/release/bioc/html/SeqVarTools.html), 
[foreach](https://cran.r-project.org/web/packages/foreach/index.html), 
[GMMAT](https://cran.r-project.org/web/packages/GMMAT/index.html), 
[CompQuadForm](https://cran.r-project.org/web/packages/CompQuadForm/index.html). 
[doMC](https://cran.r-project.org/web/packages/doMC/index.html) is required to run parallel computing.

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the [instructions](https://www.intel.com/content/www/us/en/developer/articles/technical/using-onemkl-with-r.html) 
on building R with Intel MKL.

To install MARGEE from GitHub, please use

```
devtools::install_github("MengyuZhang1307/MARGEE", ref = "main")
```

## Version

The current version is 0.1.0 (Sep 28, 2022).

## License

## Contact
Please refer to the R help document of MARGEE for specific questions about each function. 
For comments, suggestions, bug reports and questions, please contact Mengyu Zhang (megnyu.zhang@uth.tmc.edu). 
For bug reports, please include an example to reproduce the problem without having to access your confidential data.

