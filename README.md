# MDBrainT
DNA methylation deconvolution for brain TME



Installation
You can install the released version of MDBrainT from github with:

``` r
install.packages("devtools")
devtools::install_github("mackaay/MDBrainT")
install.packages("nnls")
```

Example
``` r
library(nnls)
library(MDBrainT)
data(brain_sig_matrix)
data(example_brain)
data(celltype_anno)
results <- MDBrainT(sig_matrix = brain_sig_matrix, mixture_file = example_brain, cell_annotation = celltype_anno)
```
