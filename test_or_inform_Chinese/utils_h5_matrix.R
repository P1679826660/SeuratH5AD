# ==============================================================================
# 模块名称：H5AD Matrix Loader (V2 Generic)
# 功能描述：通用的 HDF5 矩阵读取器，支持 X, raw/X, layers/*
# ==============================================================================

#' 通用矩阵读取函数
#' @param group_node H5Group 或 H5Dataset (例如 file[["X"]] 或 file[["layers"]][["norm_data"]])
#' @param transpose_output 逻辑值。如果为 TRUE，将结果转置 (用于适配 Seurat Genes x Cells)
#' @return Matrix 对象 (CsparseMatrix 或 Dense)
read_matrix_node <- function(node, transpose_output = FALSE) {
  
  mat <- NULL
  
  # 1. 识别存储类型
  if (inherits(node, "H5Group")) {
    # === 稀疏矩阵 (CSR/CSC) ===
    attrs <- hdf5r::h5attr_names(node)
    encoding_attr <- if("encoding-type" %in% attrs) hdf5r::h5attr(node, "encoding-type") else "csr_matrix"
    shape_attr <- if("shape" %in% attrs) hdf5r::h5attr(node, "shape") else NULL
    
    data <- node[["data"]]$read()
    indices <- node[["indices"]]$read()
    indptr <- node[["indptr"]]$read()
    
    # 假设它是 CSR (Python默认)
    # Python CSR: (Cells, Genes). indices=cols, indptr=rows
    # R sparseMatrix 默认构造是 CSC.
    # 如果我们将 CSR 数据填入 CSC 构造函数：
    #   i (indices) -> 变成了行索引
    #   p (indptr) -> 变成了列指针
    #   结果：矩阵被物理转置了。即 (Cells x Genes) -> (Genes x Cells)
    
    if (encoding_attr == "csr_matrix") {
      # 此时 mat 已经是 Genes x Cells (如果原数据是 Cells x Genes)
      mat <- Matrix::sparseMatrix(
        i = indices + 1,
        p = indptr,
        x = data,
        dims = c(shape_attr[2], shape_attr[1]), # 反转维度
        index1 = TRUE
      )
      
      # 如果这里本身就已经转置了(变成了 Genes x Cells)，而用户还要求 transpose_output
      # 这里的逻辑比较绕。
      # 现状：CSR读入Seurat是"天然转置"的。
      # 如果 transpose_output = TRUE (意味着我们需要 Genes x Cells)，那么这里已经不需要动了。
      # 如果 transpose_output = FALSE (意味着我们需要保持 Cells x Genes)，则需要手动 t() 回去。
      
      if (!transpose_output) {
        mat <- Matrix::t(mat)
      }
      
    } else if (encoding_attr == "csc_matrix") {
      # CSC 直接读入是原样 (Cells x Genes)
      mat <- Matrix::sparseMatrix(
        i = indices + 1,
        p = indptr,
        x = data,
        dims = c(shape_attr[1], shape_attr[2]),
        index1 = TRUE
      )
      if (transpose_output) {
        mat <- Matrix::t(mat)
      }
    }
    
  } else {
    # === 稠密矩阵 (Dataset) ===
    # HDF5 Dataset 读取到 R 会保留维度顺序，但 R 是列优先。
    # 通常读出来是 (Genes, Cells) 如果原存是 (Cells, Genes) 且无特殊属性。
    # 为稳妥起见，先读成标准矩阵
    raw_data <- node[,]
    
    # 转换为稀疏以统一处理 (除非特别要求 dense，这里先统一转)
    mat <- as(raw_data, "CsparseMatrix")
    
    # 如果用户要求转置
    if (transpose_output) {
      mat <- Matrix::t(mat)
    }
  }
  
  return(mat)
}

# 辅助：获取节点属性
h5_node_attributes <- function(node) {
  if(!exists("h5attr_names", where = asNamespace("hdf5r"))) return(list())
  
  n <- hdf5r::h5attr_names(node)
  res <- list()
  for (i in n) res[[i]] <- hdf5r::h5attr(node, i)
  return(res)
}