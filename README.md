# highSpaClone

`highSpaClone` is a computational framework for spatial CNA inference and tumor subclone classification from high-resolution spatial transcriptomics (SRT) data. The method integrates two key features: CNAs, which define genetically distinct tumor subclones, and spatial proximity, which reflects the tendency of related clones to aggregate within the tumor microenvironment. By jointly modeling these features, highSpaClone assigns subclonal labels to individual cells and reconstructs the spatial CNA landscape of whole tissue sections. The algorithm is implemented within an efficient optimization framework that scales to datasets comprising hundreds of thousands of cells, enabling accurate and interpretable analysis across diverse SRT platforms.
<img width="4500" height="3000" alt="Figure 1" src="https://github.com/user-attachments/assets/2fd11815-8f40-4692-9465-90296ebf4abe" />

## Dependencies
- R 4.3+
- Package dependencies: Rcpp, RcppArmadillo, Matrix, dplyr, magrittr, ComplexHeatmap, RANN, circlize, cluster, ggplot2, stats, utils, methods, parallel

## Installation
The R package can be installed from github:
```R
# Install devtools, if necessary
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github('CZang409/highSpaClone')

# load package
library(highSpaClone)
```

## Tutorial
Please see the [highSpaClone tutorial website](https://czang409.github.io/highSpaClone/).

