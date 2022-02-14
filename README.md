# GEDI: Gene expression decomposition and integration

Data and code used to run GEDI on a pbmc data from [https://doi.org/10.1038/s41587-020-0465-8] This data is peripheral mono nuclear cells profiled with single cell RNA seq across different sequencing technologies. We have multiple batches ( each sequencing technology ) for the same set of cell types.

1. process_data.Rmd: Create a Single Cell Experiment Object from a raw count matrix and metadata and obtain most variable genes. 
2. run_gedi.Rmd : Running GEDI with raw counts and a C matrix. Performs integration, interpretation and imputation of gene expression. 

#### **Dependencies:** 

R dependencies for GEDI:

  * Rcpp
  * RcppEigen
  * quadprog
  * ClusterR
  * Matrix
  * rsvd
  * ggplot2
  * scales

Note, RcppEigen also needs gfortran to be installed on the system.

Other R dependencies (used for the notebooks):

  * batchelor
  * uwot
  * scran
  * scater
  * ggrepel
  * ggrastr
  * tictoc
  * pheatmap
  * viridis
  * RColorBrewer

#### **Data:** 

* raw_counts_initial.rds [ raw counts of pbmcs ]
* meta_initial.rds [ metadata of pbmcs ]
* C_matrices/CellMarker.rds [http://bio-bigdata.hrbmu.edu.cn/CellMarker/]
* C_matrices/C_hallmark_msigdb.rds [ http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C2 ]
* C_matrices/C_HCL.rds [ https://doi.org/10.1038/s41586-020-2157-4 ]

Note: The original pbmc data contained 2 experiments, pbmc1 and pbmc2. The data for this tutorial comes from pbmc1.

## **Usage:**  

Load GEDI functions

```{r}
source("gedi/scIntegration.v98.svdC.R") # load GEDI functions
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
* H: The gene-level prior for unwanted sources of variation. If NULL, there will be no prior for Qi 

For this example, we are going to use GEDI with the raw counts, Bsphere mode and a matrix of cell type markers ( from cellmarker database) as prior biological information.

```{r}

K<- 100
itelim <- 150
mode <- "Bsphere"
oi_shrinkage<- 0.001

model <- new("GEDI") # Initialize GEDI object
model$setup( Samples = sce$Sample, M = raw_counts, C=c_mat, mode = mode, K = K, oi_shrinkage=oi_shrinkage ) # set up parameters
model$initialize.LVs(randomSeed = 1) # initialize LVs
tic("Optimization")
model$optimize(itelim) # run model
toc()

```

For a full tutorial, please see the run_gedi.html example. 

