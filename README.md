# SeuratH5AD

**SeuratH5AD** is a lightweight, zero-Python dependency R package for lossless conversion between Seurat V5 and AnnData (.h5ad) formats.

## Core Strengths (Why use this?)

1. * Pure R : Based on 'hdf5r', no need to install Python/Anacoda, no need to configure reticulate, no need to create conda env.
2. * Smart Transpose : Built-in DNA feature scoring algorithm to automatically identify and correct the common dimension confusion of "Genes x Cells" vs "Cells x Genes" in H5AD.
3. * High robustness :
* Automatic processing of CSR/CSC/Dense matrix format interconversion.
* Support 'raw' slot reading and writing.
* Keep 'obsm' (PCA/UMAP) and 'layers' intact.
4. * Efficient compression : GZIP Level 4 is used by default when exporting, which significantly reduces the file size without losing accuracy.

## Installation

```r
# devtools needs to be installed
devtools::install_github("P1679826660/SeuratH5AD",subdir = "SeuratH5AD")

library(SeuratH5AD)
# read
seu <- h5ad_to_seurat("test.h5ad")
# exp
seurat_to_h5ad(seu, "final_product.h5ad")
