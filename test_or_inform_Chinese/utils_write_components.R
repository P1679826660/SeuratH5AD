# ==============================================================================
# 模块名称：H5AD Components Writer (V2 Fixed)
# 功能描述：写入 obsm (Embeddings) 和 layers，修复了函数导出问题
# ==============================================================================

#' 写入降维信息 (Reductions -> obsm)
#' @param h5_file HDF5 文件对象
#' @param seurat_obj Seurat 对象
write_obsm_h5 <- function(h5_file, seurat_obj) {
  reductions <- names(seurat_obj@reductions)
  if (length(reductions) == 0) return()
  
  # 如果 obsm 组已存在先删除，防止重名报错
  if (h5_file$exists("obsm")) h5_file$link_delete("obsm")
  obsm_grp <- h5_file$create_group("obsm")
  
  for (red in reductions) {
    # 获取 Embeddings: Cells x Dims
    emb <- Seurat::Embeddings(seurat_obj[[red]])
    
    # 构造键名: pca -> X_pca, umap -> X_umap
    key_name <- paste0("X_", red)
    
    # 写入稠密矩阵 (需要转置以适配 Python 读序)
    obsm_grp$create_dataset(key_name, t(emb))
  }
}

#' 写入 Layers (Assay Data -> layers)
#' @param h5_file HDF5 文件对象
#' @param seurat_obj Seurat 对象
write_layers_h5 <- function(h5_file, seurat_obj) {
  # 假设只处理 Default Assay
  assay_name <- Seurat::DefaultAssay(seurat_obj)
  assay_obj <- seurat_obj[[assay_name]]
  
  # [修复点 1] 获取 layers 名称的鲁棒写法
  # 尝试从对象内部直接获取 names(assay@layers)
  layer_names <- names(assay_obj@layers)
  
  if (length(layer_names) == 0) return()
  
  # 创建 layers 组
  if (!h5_file$exists("layers")) layers_grp <- h5_file$create_group("layers")
  else layers_grp <- h5_file[["layers"]]
  
  for (lname in layer_names) {
    # 跳过 data (因为它通常去 X)
    if (lname == "data") next
    
    # [修复点 2] 使用 GetAssayData 替代 LayerData
    mat <- Seurat::GetAssayData(seurat_obj, layer = lname, assay = assay_name)
    
    # 命名映射
    save_name <- lname
    if (lname == "scale.data") save_name <- "scale_data" # 避免点号
    
    # 写入
    # scale.data 通常是稠密的
    if (inherits(mat, "matrix")) {
      # 稠密写入
      if(layers_grp$exists(save_name)) layers_grp$link_delete(save_name)
      layers_grp$create_dataset(save_name, t(mat))
    } else {
      # 稀疏写入 (利用 matrix 模块的函数)
      if(layers_grp$exists(save_name)) layers_grp$link_delete(save_name)
      sub_grp <- layers_grp$create_group(save_name)
      write_matrix_h5(sub_grp, mat)
    }
  }
}