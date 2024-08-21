# GEDI: Gene expression decomposition and integration

## Content

- [Overview](#Overview)
- [Requirements](#Requirements)
- [Installation](#Installation)
- [Usage](#Usage)
- [License](./LICENSE.md)
- [Preprint](#Preprint)
- [Reproducible analysis for manuscript](https://github.com/csglab/GEDI_manuscript)
- [Contact](#Contact)

# Overview

A generative model that unifies integration, cluster-free differential expression, pathway and regulatory network analysis, data normalization and imputation. 

# Requirements

The following R packages are required: 

  * [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
  * [RcppEigen](https://cran.r-project.org/web/packages/RcppEigen/index.html)
  * [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
  * [metR](https://cran.r-project.org/web/packages/metR/index.html)
  * [rsvd](https://cran.r-project.org/web/packages/rsvd/index.html)
  * [scales](https://cran.r-project.org/web/packages/scales/index.html)

Note, RcppEigen also needs [gfortran](https://fortran-lang.org/learn/os_setup/install_gfortran/) to be installed on the system.

Other R dependencies (used for the notebooks):

  * [uwot](https://cran.r-project.org/web/packages/uwot/index.html)
  * [scran](https://bioconductor.org/packages/release/bioc/html/scran.html)
  * [scater](https://bioconductor.org/packages/release/bioc/html/scater.html)
  
The code was tested using R 4.0.0 running on CentOS Linux 7.

`GEDI` should be compatible with Windows, Mac, and Linux operating systems.

# Installation

From an `R` session, type:

```{r}
devtools::install_github("csglab/GEDI")
```
  
# Usage

Check the following notebooks for examples on how to run `GEDI`: 

* [Quick intro](https://csglab.github.io/GEDI/vignettes/GEDI_quick.html)
* [Analysis of sample to sample variability](https://csglab.github.io/GEDI/vignettes/GEDI_sample_PBMC.html)
* [TF analysis in PBMC](https://csglab.github.io/GEDI/vignettes/GEDI_tf_analysis_PBMC.html)
* [Ratio mode: splicing](https://csglab.github.io/GEDI/vignettes/GEDI_splicing.html)

# Article

[Madrigal, A., Lu, T., Soto, L. M., & Najafabadi, H. S. (2024). A unified model for interpretable latent embedding of multi-sample, multi-condition single-cell data. Nature Communications, 15(1), 6573.](https://www.nature.com/articles/s41467-024-50963-0)

# Contact

We use GitHub [issues](https://github.com/csglab/GEDI/issues) for tracking requests and bugs. Please submit a new issue if you have any comment or you would like to report a software bug. 
