#' Export Seurat V5 to H5AD
#' @param seurat_obj Seurat Object
#' @param file_path Output path
#' @param compression_level 0-9
#' @param write_raw_slot boolean
#' @export
seurat_to_h5ad <- function(seurat_obj, file_path, compression_level = 4, write_raw_slot = TRUE) {
  message('>>> [Export] Start: ', file_path)
  if (file.exists(file_path)) file.remove(file_path)
  file_h5 <- hdf5r::H5File$new(file_path, mode = 'w')
  on.exit(file_h5$close_all(), add = TRUE)
  assay_name <- Seurat::DefaultAssay(seurat_obj)
  layers_avail <- names(seurat_obj[[assay_name]]@layers)
  main_layer <- if ('data' %in% layers_avail) 'data' else 'counts'
  main_mat <- Seurat::GetAssayData(seurat_obj, layer = main_layer, assay = assay_name)
  x_grp <- file_h5$create_group('X')
  write_matrix_h5(x_grp, main_mat, compression_level)
  write_dataframe_h5(file_h5, 'obs', seurat_obj@meta.data)
  feat_meta <- seurat_obj[[assay_name]][[]]
  if(nrow(feat_meta)==0) feat_meta <- data.frame(row.names=rownames(seurat_obj))
  write_dataframe_h5(file_h5, 'var', feat_meta)
  if (write_raw_slot && 'counts' %in% layers_avail) {
    raw_grp <- file_h5$create_group('raw'); raw_x_grp <- raw_grp$create_group('X')
    counts_mat <- Seurat::GetAssayData(seurat_obj, layer = 'counts', assay = assay_name)
    write_matrix_h5(raw_x_grp, counts_mat, compression_level)
    write_dataframe_h5(raw_grp, 'var', feat_meta)
  }
  write_obsm_h5(file_h5, seurat_obj)
  write_layers_h5(file_h5, seurat_obj)
  message('>>> Export Success!')
}
