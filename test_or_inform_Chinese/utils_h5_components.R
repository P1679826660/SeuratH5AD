# ==============================================================================
# 模块名称：H5AD Components Loader
# 功能描述：读取 Layers, Embeddings (obsm), Loadings (varm), Graphs (obsp)
# ==============================================================================

#' 读取并注入降维信息 (Reductions)
#' @param seu Seurat 对象
#' @param h5_file H5File 句柄
#' @param transpose_mode 也就是 convert_h5ad 中的 "转置模式" 标记
add_reductions <- function(seu, h5_file, transpose_mode) {
  
  # === 1. 确定谁存了 Cell Embeddings (坐标) ===
  # 标准模式: obsm 存坐标 (Cells x Dims)
  # 转置模式: varm 存坐标 (因为 var 变成了 Cells)
  emb_group_name <- if(transpose_mode) "varm" else "obsm"
  loadings_group_name <- if(transpose_mode) "obsm" else "varm"
  
  cell_names <- colnames(seu)
  
  if (h5_file$exists(emb_group_name)) {
    grp <- h5_file[[emb_group_name]]
    keys <- grp$ls(recursive=FALSE)$name
    
    for (k in keys) {
      # 过滤掉非矩阵数据
      if (!inherits(grp[[k]], "H5D") && !inherits(grp[[k]], "H5Group")) next 
      
      # 读取坐标
      # 注意：obsm 读取通常是 (N_cells, N_dims)。
      # 如果是 Dense Dataset，读入 R 后可能是 (N_dims, N_cells) 因为 R 是列优先读取 HDF5
      raw_emb <- tryCatch(grp[[k]]$read(), error=function(e) NULL)
      if (is.null(raw_emb)) next
      
      # 确保它是矩阵
      if (!is.matrix(raw_emb)) raw_emb <- as.matrix(raw_emb)
      
      # 维度对齐检查
      # 我们需要行数 = 细胞数
      if (ncol(raw_emb) == length(cell_names)) {
        raw_emb <- t(raw_emb)
      }
      
      if (nrow(raw_emb) != length(cell_names)) {
        # 维度依然对不上，跳过
        next
      }
      
      # 命名规范化 X_pca -> pca
      clean_name <- tolower(gsub("^X_", "", k))
      seurat_key <- paste0(clean_name, "_")
      
      rownames(raw_emb) <- cell_names
      colnames(raw_emb) <- paste0(seurat_key, 1:ncol(raw_emb))
      
      # 创建 DR 对象
      dr_obj <- Seurat::CreateDimReducObject(
        embeddings = raw_emb,
        key = seurat_key,
        assay = Seurat::DefaultAssay(seu)
      )
      
      # === 2. 尝试读取对应的 Feature Loadings ===
      # 如果有对应的 varm (基因权重)
      if (h5_file$exists(loadings_group_name)) {
        # 猜测名称：PCA 对应 PCs, usually keys match or are similar
        # Scanpy convention: obsm['X_pca'] <-> varm['PCs']
        loadings_key <- NULL
        if (clean_name == "pca" && h5_file[[loadings_group_name]]$exists("PCs")) loadings_key <- "PCs"
        else if (h5_file[[loadings_group_name]]$exists(k)) loadings_key <- k
        
        if (!is.null(loadings_key)) {
          raw_load <- h5_file[[loadings_group_name]][[loadings_key]]$read()
          if (!is.matrix(raw_load)) raw_load <- as.matrix(raw_load)
          
          # Loadings 需要是 (Features x Dims)
          # 检查维度
          n_features <- nrow(seu)
          if (ncol(raw_load) == n_features) raw_load <- t(raw_load)
          
          if (nrow(raw_load) == n_features) {
            rownames(raw_load) <- rownames(seu)
            colnames(raw_load) <- paste0(seurat_key, 1:ncol(raw_load))
            
            # 注入 Loadings
            Seurat::Loadings(dr_obj) <- raw_load
          }
        }
      }
      
      seu[[clean_name]] <- dr_obj
    }
  }
  return(seu)
}

#' 读取并注入额外 Layers (data, scale.data)
add_layers <- function(seu, h5_file, transpose_mode) {
  if (!h5_file$exists("layers")) return(seu)
  
  grp <- h5_file[["layers"]]
  layer_names <- grp$ls(recursive=FALSE)$name
  
  default_assay <- Seurat::DefaultAssay(seu)
  
  for (lname in layer_names) {
    # 决定是否需要转置
    # 如果标准模式 (obs=cells), layers 通常是 Cells x Genes。Seurat Layer 需要 Genes x Cells。
    # 所以标准模式下，read_matrix_node 需要 transpose_output=FALSE (因为内部 logic 会自动转 CSR)
    # 还是直接利用 read_matrix_node 的智能逻辑？
    # 简单点：我们希望得到 Genes x Cells。
    
    # 这里的逻辑必须跟 convert_h5ad.R 里的 Main Matrix 读取保持一致。
    # 如果 Main Matrix 转置了，Layer 也要转置。
    
    # 这里的 transpose_output 参数意为：是否需要在读取结果基础上再 T 一次
    # 如果文件本身是 Standard (Cells x Genes)，读入 CSR 自动变 (Genes x Cells)，我们不需要再 T。
    # 如果文件是 Transposed (Genes x Cells)，读入 CSR 自动变 (Cells x Genes)，我们需要 T 变成 Genes x Cells。
    
    do_transpose <- transpose_mode 
    
    mat <- read_matrix_node(grp[[lname]], transpose_output = do_transpose)
    
    # 维度二次校验
    if (ncol(mat) != ncol(seu) || nrow(mat) != nrow(seu)) {
      # 如果维度不对，尝试转置
      mat <- Matrix::t(mat)
    }
    
    if (ncol(mat) == ncol(seu) && nrow(mat) == nrow(seu)) {
      # 命名映射
      # Scanpy: 'logcounts' -> Seurat: 'data'
      # Scanpy: 'scale' -> Seurat: 'scale.data'
      
      target_slot <- lname
      if (lname %in% c("logcounts", "norm_data", "normalize")) target_slot <- "data"
      if (grepl("scale", lname)) target_slot <- "scale.data"
      
      # 注入
      if (target_slot == "scale.data") {
        # scale.data 必须是 dense matrix
        seu <- Seurat::SetAssayData(seu, slot = "scale.data", new.data = as.matrix(mat), assay = default_assay)
      } else if (target_slot == "data") {
        seu <- Seurat::SetAssayData(seu, slot = "data", new.data = mat, assay = default_assay)
      } else {
        # 其他 layer，如 'spliced'
        Seurat::LayerData(seu, assay = default_assay, layer = lname) <- mat
      }
      message(paste("    已添加 Layer:", lname, "->", target_slot))
    }
  }
  return(seu)
}