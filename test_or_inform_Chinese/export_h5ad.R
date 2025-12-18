# ==============================================================================
# Module Name: Seurat to H5AD Exporter (Main V2 Fixed)
# Function Description: Main controller program, fixes namespace export errors
# ==============================================================================

#' Export Seurat V5 object to .h5ad file
#' @param seurat_obj Seurat Object
#' @param file_path Output path (e.g., "output.h5ad")
#' @export
seurat_to_h5ad <- function(seurat_obj, file_path) {
  
  message(">>> [Export] Starting export of Seurat object to: ", file_path)
  
  # 0. Prepare environment
  if (file.exists(file_path)) file.remove(file_path)
  file_h5 <- hdf5r::H5File$new(file_path, mode = "w")
  on.exit(file_h5$close_all(), add = TRUE)
  
  assay_name <- Seurat::DefaultAssay(seurat_obj)
  message("    - Using Default Assay: ", assay_name)
  
  # 1. Write X (Main Matrix)
  message(">>> [Phase 1] Writing main matrix X ...")
  
  # [Fix Point 1] Get available layers
  layers_avail <- names(seurat_obj[[assay_name]]@layers)
  
  # Decision: Prioritize 'data' if available; otherwise use 'counts'
  main_layer <- "counts" # Default fallback
  if ("data" %in% layers_avail) {
    main_layer <- "data"
  }
  
  message(sprintf("    - Selected Layer '%s' as X", main_layer))
  
  # [Fix Point 2] Use GetAssayData instead of LayerData
  main_mat <- Seurat::GetAssayData(seurat_obj, layer = main_layer, assay = assay_name)
  
  # Create X group
  x_grp <- file_h5$create_group("X")
  # Call matrix writing module (handles transposition and CSC format automatically)
  write_matrix_h5(x_grp, main_mat)
  
  # 2. Write obs (Cell Metadata)
  message(">>> [Phase 2] Writing obs (Cell Metadata)...")
  write_dataframe_h5(file_h5, "obs", seurat_obj@meta.data)
  
  # 3. Write var (Feature Metadata)
  message(">>> [Phase 3] Writing var (Feature Metadata)...")
  # In Seurat V5, feature metadata is located in the assay's meta.data
  feat_meta <- seurat_obj[[assay_name]][[]]
  
  # Ensure row names exist (Gene names)
  if (nrow(feat_meta) == 0) {
    # If no metadata, create a DF with indices only
    feat_meta <- data.frame(row.names = rownames(seurat_obj))
  }
  write_dataframe_h5(file_h5, "var", feat_meta)
  
  # 4. Write obsm (Embeddings)
  message(">>> [Phase 4] Writing obsm (Embeddings)...")
  write_obsm_h5(file_h5, seurat_obj)
  
  # 5. Write layers (Additional Layers)
  message(">>> [Phase 5] Writing layers (Additional Layers)...")
  write_layers_h5(file_h5, seurat_obj)
  
  message(">>> Export successful!")
  return(invisible(TRUE))
}
