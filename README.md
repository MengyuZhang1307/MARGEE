# RoM (Robust Mixed-model Association tests for Gene–Environment Interaction)

## Description

RoM is an R package for Mixed-model association test. Perform generalized linear mixed robust single-variant gene-environment interaction tests and joint tests for common variants.

## Installing

RoM imports R packages 
[Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), 
[parallel](https://cran.r-project.org/web/views/HighPerformanceComputing.html), 
[MASS](https://cran.r-project.org/web/packages/MASS/index.html), 
[SeqArray](http://bioconductor.org/packages/release/bioc/html/SeqArray.html), 
[SeqVarTools](https://bioconductor.org/packages/release/bioc/html/SeqVarTools.html), 
[foreach](https://cran.r-project.org/web/packages/foreach/index.html), 
[GMMAT](https://cran.r-project.org/web/packages/GMMAT/index.html),.
[doMC](https://cran.r-project.org/web/packages/doMC/index.html) is required to run parallel computing.

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the [instructions](https://www.intel.com/content/www/us/en/developer/articles/technical/using-onemkl-with-r.html) 
on building R with Intel MKL.

To install RoM from GitHub, please use

```
devtools::install_github("MengyuZhang1307/RoM", ref = "main")
```

## Version

The current version is 0.1.0 (Oct 10, 2022).

## License

This software is licensed under GPL-3.

## Contact

Please refer to the R help document of RoM for specific questions about each function. 
For comments, suggestions, bug reports and questions, please contact Mengyu Zhang (mengyu.zhang@uth.tmc.edu). 
For bug reports, please include an example to reproduce the problem without having to access your confidential data.

