# highSpaClone

`highSpaClone` is a computational framework for spatial CNV inference and tumor subclone classification from high-resolution spatial transcriptomics (SRT) data. The method integrates two key features: CNVs, which define genetically distinct tumor subclones, and spatial proximity, which reflects the tendency of related clones to aggregate within the tumor microenvironment. By jointly modeling these features, highSpaClone assigns subclonal labels to individual cells and reconstructs the spatial CNV landscape of whole tissue sections. The algorithm is implemented within an efficient optimization framework that scales to datasets comprising hundreds of thousands of cells, enabling accurate and interpretable analysis across diverse SRT platforms.

<img width="4500" height="3000" alt="IRIS Figure 1 (10)" src="https://github.com/user-attachments/assets/416efb15-fe2c-47aa-8331-272c06922a59" />

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

