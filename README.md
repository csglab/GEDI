# GEDI

Data and code used to run GEDI on a pbmc data from [https://doi.org/10.1038/s41587-020-0465-8] This data is peripheral mono nuclear cells profiled with single cell RNA seq across different sequencing technologies. We have multiple batches ( each sequencing technology ) for the same set of cell types.

1. process_data.Rmd: Create a Single Cell Experiment Object from a raw count matrix and metadata and obtain most variable genes. 
2. run_gedi.Rmd : Running GEDI with raw counts and a C matrix. Performs integration, interpretation and imputation of gene expression. 

Dependencies for GEDI:

R dependencies:

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

Data:

* raw_counts_initial.rds [ raw counts of pbmcs ]
* meta_initial.rds [ metadata of pbmcs ]
* C_matrices/CellMarker.rds [http://bio-bigdata.hrbmu.edu.cn/CellMarker/]
* C_matrices/C_hallmark_msigdb.rds [ http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C2 ]
* C_matrices/C_HCL.rds [ https://doi.org/10.1038/s41586-020-2157-4 ]


Note: The original pbmc data contained 2 experiments, pbmc1 and pbmc2. The data for this tutorial comes from pbmc1.
