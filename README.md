# SeuratH5AD

**SeuratH5AD** is a lightweight, zero-Python dependency R package for lossless conversion between Seurat V5 and AnnData (.h5ad) formats.

## Installation

```r
# devtools needs to be installed
devtools::install_github("P1679826660/SeuratH5AD",subdir = "SeuratH5AD")

library(SeuratH5AD)
# read
seu <- h5ad_to_seurat("test.h5ad")
# exp
seurat_to_h5ad(seu, "final_product.h5ad")
