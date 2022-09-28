# MARGEE (Mixed-model Association robust tests for GEne–Environment interactions)

## Description

## Installing

MARGEE imports R packages [Rcpp], Matrix, parallel, MASS, SeqArray, SeqVarTools, foreach, GMMAT, CompQuadForm. 
[doMC] is required to run parallel computing.

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel 
Library (or other fast BLAS/LAPACK libraries). See the [instructions](https://www.intel.com/content/www/us/en/developer/articles/technical/using-onemkl-with-r.html) 
on building R with Intel MKL.

To install FiMAP from GitHub, please use

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

