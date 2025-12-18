# ==============================================================================
# 模块名称：H5AD Matrix Writer
# 功能描述：将 Seurat 稀疏矩阵转换为 AnnData 兼容的 HDF5 格式 (自动转置)
# ==============================================================================

#' 写入稀疏矩阵到 HDF5 Group
#' @param h5_group HDF5 组对象 (例如 file[["X"]] 或 file[["layers"]]$create_group("norm"))
#' @param mat Seurat 的稀疏矩阵 (Genes x Cells)
#' @description 该函数会自动转置矩阵为 (Cells x Genes) 并以 CSC 格式写入
write_matrix_h5 <- function(h5_group, mat) {
  
  # 1. 物理转置：Seurat (Genes x Cells) -> AnnData (Cells x Genes)
  # 注意：Matrix::t() 对于 dgCMatrix 会返回一个 dgCMatrix (表示转置后的矩阵)
  # 此时，行变成 Cell，列变成 Gene。
  # 存储格式依然是 CSC (Compressed Sparse Column)，即按 Gene 压缩。
  # 这对 AnnData 是合法的 (encoding-type = "csc_matrix")。
  mat_t <- Matrix::t(mat)
  
  # 2. 提取底层数据 (R 的 dgCMatrix slot 已经是 0-based，这是完美的)
  # @x: 数值
  # @i: 行索引 (Row indices)
  # @p: 列指针 (Column pointers)
  data <- mat_t@x
  indices <- mat_t@i
  indptr <- mat_t@p
  
  # 3. 写入数据集
  # 使用 create_dataset 并启用 gzip 压缩，减少体积
  h5_group$create_dataset("data", data, gzip_level = 4)
  h5_group$create_dataset("indices", indices, gzip_level = 4)
  h5_group$create_dataset("indptr", indptr, gzip_level = 4)
  
  # 4. 写入关键属性 (AnnData 识别标准)
  # shape: (n_cells, n_genes)
  h5_group$create_attr("encoding-type", "csc_matrix")
  h5_group$create_attr("encoding-version", "0.1.0")
  # 注意 HDF5 shape 顺序。h5attr 写入向量时通常是一维的。
  # Python 读取时 shape 应该是 (Cells, Genes)
  h5_group$create_attr("shape", dim(mat_t))
}

#' 写入稠密矩阵 (用于 Embeddings 等)
#' @param h5_dataset HDF5 数据集路径
write_dense_h5 <- function(h5_file, path, mat) {
  # mat: Cells x Dims
  # HDF5R 写入矩阵时，默认会保留 R 的列优先。
  # Python 读取时，维度会反转。
  # 如果我们在 R 里是 (Cells x Dims)，直接写，Python 读出来是 (Dims x Cells)。
  # 所以写入稠密矩阵前，通常需要转置一下，或者利用 hdf5r 的特性。
  # 经验法则：为了让 Python 读到 (N, M)，R 应该写 t(mat)。
  
  h5_file$create_dataset(path, t(mat), gzip_level = 4)
}