# GEDI: Gene expression decomposition and integration

## Content

- [Overview](#Overview)
- [System Requirements](#System-requirements)
- [Installation](#Installation)
- [Usage](#instructions-for-use)
- [License](./LICENSE.md)
- [Issues](https://github.com/csglab/GEDI/issues)
- [Pre-print](https://www.biorxiv.org/content/10.1101/2023.08.15.553327v1)
- [Reproducible analysis for manuscript](https://github.com/csglab/GEDI_manuscript)

# Overview

A generative model that unifies integration, cluster-free differential expression, pathway and regulatory network analysis,  data normalization and imputation. 

# System Requirements

## **Dependencies:** 

R dependencies for GEDI:

  * Rcpp
  * RcppEigen
  * ggplot2
  * methods
  * metR
  * rsvd
  * scales	

Note, RcppEigen also needs gfortran to be installed on the system.

Other R dependencies (used for the notebooks):

  * uwot
  * scran
  * scater

# Installation

From an R session `R` session, type:

```{r}

devtools::install_github("csglab/GEDI")

```
  
# **Usage:**

Load GEDI

```{r}

library(GEDI)

```

#### **Arguments:**

* Samples:  Batch variable to use. It should be a character vector.
* expression matrix: Could be one of the following:
	+ Y: The log-transformed (and possibly normalized) gene expression matrix.
	+ M: The raw read count matrix.  It can also be a list of two matrices, in which case they are considered as paired observations whose log-ratio must be modelled.
* K: The number of latent variables.
* mode: Two values are allowed: 
	+ Bl2: L2 norm of the entire B matrix is fixed. Interpretation is that we’re projecting the data on a lower-dimensional hyperplane with dimension K. 
	+ Bsphere: L2 norms of B columns are fixed. Interpretation is that we’re projecting the data on a hyperellipsoid of dimension K. 	
* itelim: Number of iterations for optimization.

**Optional arguments**:

* C: The gene-level biological prior. If NULL, it means that there is no prior for Z
* H: Sample-level prior for sources of variation. If NULL, there will be no prior for Qi 
* oi_shrinkage: Shrinkage multiplier for oi (offset vector per sample i)

For this example, we are going to use GEDI with the raw counts, Bsphere mode and a matrix of cell type markers as prior biological information.

```{r}

K<- 100
itelim <- 150
mode <- "Bsphere"

model <- new("GEDI") # Initialize GEDI object
model$setup( Samples = sce$Sample, M = raw_counts, C=c_mat, mode = mode, K = K ) # set up parameters
model$initialize.LVs(randomSeed = 1) # initialize LVs
model$optimize(itelim) # run model

```
