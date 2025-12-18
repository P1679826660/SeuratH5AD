# ==============================================================================
# 模块名称：H5AD Embeddings Reader
# 功能描述：读取 obsm 组中的降维坐标 (PCA, UMAP等) 并转换为 Seurat DimReduc 对象
# ==============================================================================

#' 读取 obsm 组并转换为 Seurat 降维对象列表
#' @param h5_file H5File 对象
#' @param seurat_object 已经创建好的 Seurat 对象 (用于校验细胞名)
#' @param transpose_mode 布尔值，指示是否发生了转置 (如果转置了，obsm 也需要适配)
#' @return 更新后的 Seurat 对象
read_h5ad_obsm <- function(h5_file, seurat_object, transpose_mode = FALSE) {
  if (!h5_file$exists("obsm")) return(seurat_object)
  
  grp <- h5_file[["obsm"]]
  dataset_names <- grp$ls(recursive = FALSE)$name
  
  # 获取细胞名用于校验
  cell_names <- colnames(seurat_object)
  n_cells <- length(cell_names)
  
  message(sprintf(">>> [Phase 6] 发现 obsm (降维信息): %s", paste(dataset_names, collapse = ", ")))
  
  for (name in dataset_names) {
    # 读取坐标矩阵
    # obsm 通常是 (N_cells x N_dims) 的稠密矩阵
    coords <- grp[[name]]$read()
    
    # 维度检查与修正
    # 如果发生了 transpose_mode (原文件 obs是基因)，那么理论上 obsm 应该跟着 cells (原 var) 走
    # 但 H5AD 标准中 obsm 总是与 obs 对齐。
    # 这是一个极其罕见的边界情况：如果 h5ad 整体反了，obsm 可能也需要特殊处理。
    # 但通常 obsm 的行数 = 细胞数。
    
    # 确保矩阵是 matrix 类型
    if (!is.matrix(coords)) coords <- as.matrix(coords)
    
    # 检查行数是否匹配细胞数
    if (nrow(coords) != n_cells) {
      # 尝试转置
      if (ncol(coords) == n_cells) {
        coords <- t(coords)
      } else {
        warning(sprintf("    跳过 %s: 维度 (%d x %d) 与细胞数 (%d) 不匹配。", name, nrow(coords), ncol(coords), n_cells))
        next
      }
    }
    
    # 清洗名称：X_umap -> umap, X_pca -> pca
    key_clean <- gsub("^X_", "", name)
    key_clean <- tolower(key_clean)
    
    # Seurat Key 要求 (例如 "umap_")
    seurat_key <- paste0(key_clean, "_")
    
    # 给坐标赋予行名 (必须与细胞名一致)
    rownames(coords) <- cell_names
    # 给坐标赋予列名 (PC_1, UMAP_1...)
    colnames(coords) <- paste0(seurat_key, 1:ncol(coords))
    
    # 创建 DimReduc 对象
    # 注意：我们这里没有 feature.loadings (varm)，通常 h5ad 不一定有，或者在 varm 里
    # 为了兼容性，这里只创建 Embeddings
    dim_reduc <- Seurat::CreateDimReducObject(
      embeddings = coords,
      key = seurat_key,
      assay = Seurat::DefaultAssay(seurat_object)
    )
    
    # 注入到 Seurat 对象
    seurat_object[[key_clean]] <- dim_reduc
  }
  
  return(seurat_object)
}