# ==============================================================================
# 模块名称：Seurat to H5AD Exporter (Main V2 Fixed)
# 功能描述：主控程序，修复 namespace 导出报错
# ==============================================================================

#' 将 Seurat V5 对象导出为 .h5ad 文件
#' @param seurat_obj Seurat 对象
#' @param file_path 输出路径 (例如 "output.h5ad")
#' @export
seurat_to_h5ad <- function(seurat_obj, file_path) {
  
  message(">>> [Export] 开始导出 Seurat 对象到: ", file_path)
  
  # 0. 准备环境
  if (file.exists(file_path)) file.remove(file_path)
  file_h5 <- hdf5r::H5File$new(file_path, mode = "w")
  on.exit(file_h5$close_all(), add = TRUE)
  
  assay_name <- Seurat::DefaultAssay(seurat_obj)
  message("    - 使用 Default Assay: ", assay_name)
  
  # 1. 写入 X (主矩阵)
  message(">>> [Phase 1] 写入主矩阵 X ...")
  
  # [修复点 1] 获取可用 layers
  layers_avail <- names(seurat_obj[[assay_name]]@layers)
  
  # 决策：如果有 'data'，优先用 'data'；否则用 'counts'
  main_layer <- "counts" # 默认 fallback
  if ("data" %in% layers_avail) {
    main_layer <- "data"
  }
  
  message(sprintf("    - 选定 Layer '%s' 作为 X", main_layer))
  
  # [修复点 2] 使用 GetAssayData 替代 LayerData
  main_mat <- Seurat::GetAssayData(seurat_obj, layer = main_layer, assay = assay_name)
  
  # 创建 X 组
  x_grp <- file_h5$create_group("X")
  # 调用矩阵写入模块 (自动处理转置和 CSC 格式)
  write_matrix_h5(x_grp, main_mat)
  
  # 2. 写入 obs (细胞元数据)
  message(">>> [Phase 2] 写入 obs (Cell Metadata)...")
  write_dataframe_h5(file_h5, "obs", seurat_obj@meta.data)
  
  # 3. 写入 var (基因元数据)
  message(">>> [Phase 3] 写入 var (Feature Metadata)...")
  # Seurat V5 中 feature metadata 在 assay 的 meta.data 里
  feat_meta <- seurat_obj[[assay_name]][[]]
  
  # 确保行名存在 (基因名)
  if (nrow(feat_meta) == 0) {
    # 如果没元数据，至少创建一个只有索引的 DF
    feat_meta <- data.frame(row.names = rownames(seurat_obj))
  }
  write_dataframe_h5(file_h5, "var", feat_meta)
  
  # 4. 写入 obsm (降维)
  message(">>> [Phase 4] 写入 obsm (Embeddings)...")
  write_obsm_h5(file_h5, seurat_obj)
  
  # 5. 写入 layers (其他矩阵)
  message(">>> [Phase 5] 写入 layers (Additional Layers)...")
  write_layers_h5(file_h5, seurat_obj)
  
  message(">>> 导出成功！")
  return(invisible(TRUE))
}