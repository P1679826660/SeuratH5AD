#' Export Seurat V5 to H5AD
#' @param seurat_obj Seurat Object
#' @param file_path Output path
#' @param compression_level Integer 0-9
#' @param write_raw_slot boolean
#' @export
seurat_to_h5ad <- function(seurat_obj, file_path, compression_level = 4, write_raw_slot = TRUE) {
  message('>>> [Export] Start: ', file_path)

  if (!inherits(seurat_obj, "Seurat")) stop("Input is not a Seurat object")

  if (file.exists(file_path)) {
    tryCatch(file.remove(file_path), error = function(e) stop("Cannot remove existing file: ", file_path))
  }

  file_h5 <- hdf5r::H5File$new(file_path, mode = 'w')
  on.exit(file_h5$close_all(), add = TRUE)

  assay_name <- Seurat::DefaultAssay(seurat_obj)

  # Determine main layer for X
  # For Seurat V5, check layers explicitly
  if (inherits(seurat_obj[[assay_name]], "Assay5")) {
     layers_avail <- SeuratObject::Layers(seurat_obj, assay = assay_name)
     main_layer <- if ('data' %in% layers_avail) 'data' else 'counts'
  } else {
     # Seurat V3/V4
     main_layer <- 'data'
  }

  main_mat <- tryCatch({
    Seurat::GetAssayData(seurat_obj, layer = main_layer, assay = assay_name)
  }, error = function(e) {
    Seurat::GetAssayData(seurat_obj, slot = 'data', assay = assay_name)
  })

  x_grp <- file_h5$create_group('X')
  write_matrix_h5(x_grp, main_mat, compression_level)

  write_dataframe_h5(file_h5, 'obs', seurat_obj@meta.data)

  feat_meta <- tryCatch(seurat_obj[[assay_name]][[]], error = function(e) data.frame())
  if (nrow(feat_meta) == 0) feat_meta <- data.frame(row.names = rownames(seurat_obj))
  write_dataframe_h5(file_h5, 'var', feat_meta)

  # Write Raw counts if requested and available
  has_counts <- tryCatch({
    !is.null(Seurat::GetAssayData(seurat_obj, slot = 'counts', assay = assay_name))
  }, error = function(e) FALSE)

  if (write_raw_slot && has_counts) {
    counts_mat <- Seurat::GetAssayData(seurat_obj, slot = 'counts', assay = assay_name)
    raw_grp <- file_h5$create_group('raw')
    raw_x_grp <- raw_grp$create_group('X')
    write_matrix_h5(raw_x_grp, counts_mat, compression_level)
    write_dataframe_h5(raw_grp, 'var', feat_meta)
  }

  write_obsm_h5(file_h5, seurat_obj)
  write_layers_h5(file_h5, seurat_obj)

  message('>>> Export Success!')
}
