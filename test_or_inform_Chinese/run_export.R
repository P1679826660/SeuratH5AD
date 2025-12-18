# run_export.R

library(hdf5r)
library(Seurat)
library(Matrix)


source("utils_write_matrix.R")
source("utils_write_metadata.R")
source("utils_write_components.R")
source("export_h5ad.R")


output_file <- "export_test.h5ad"
seurat_to_h5ad(seurat_obj, output_file)
