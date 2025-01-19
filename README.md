# MDBrainT
DNA methylation deconvolution for brain TME

CIBERSORT
The goal of CIBERSORT is to run the CIBERSORT flow.

Installation
You can install the released version of CIBERSORT from github with:

``` r
install.packages("devtools")
devtools::install_github("mackaay/MDBrainT")
```

Example
``` r
library(MDBrainT)
data(brain_sig_matrix)
data(example_brain)
data(celltype_anno)
results <- MDBrainT(sig_matrix = brain_sig_matrix, mixture_file = example_brain, cell_annotation = celltype_anno)
```
