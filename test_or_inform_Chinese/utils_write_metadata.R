# ==============================================================================
# 模块名称：H5AD Metadata Writer
# 功能描述：将 DataFrame 写入 HDF5 Group，自动处理行名索引
# ==============================================================================

#' 写入 DataFrame 到 HDF5 Group
#' @param h5_file HDF5 文件对象
#' @param group_name 组名 ("obs" 或 "var")
#' @param df R 的 data.frame
#' @param index_name 索引列的名称 (通常是 "_index")
write_dataframe_h5 <- function(h5_file, group_name, df, index_name = "_index") {
  
  if (is.null(df)) return()
  
  # 创建组
  if (h5_file$exists(group_name)) h5_file$link_delete(group_name)
  grp <- h5_file$create_group(group_name)
  
  # 1. 写入行名 (Index)
  # AnnData 要求行名作为一列存储，并在属性中引用它
  row_names <- rownames(df)
  if (is.null(row_names)) row_names <- as.character(1:nrow(df))
  
  # 写入索引列
  grp$create_dataset(index_name, row_names)
  
  # 标记这是一个数据框，并且索引是 index_name
  grp$create_attr("_index", index_name)
  grp$create_attr("encoding-type", "dataframe")
  grp$create_attr("encoding-version", "0.2.0")
  
  # 写入列顺序 (可选，但在 Python Pandas 中有用)
  col_order <- c(index_name, colnames(df))
  grp$create_attr("column-order", col_order)
  
  # 2. 写入每一列
  for (col in colnames(df)) {
    val <- df[[col]]
    
    # 类型转换：Factor -> Character (为了最大兼容性)
    if (is.factor(val)) {
      val <- as.character(val)
    }
    
    # 处理 NA：HDF5 不支持字符串 NA。转为空字符串或 "NA"
    if (is.character(val)) {
      val[is.na(val)] <- ""
    }
    
    # 写入数据集
    # 注意：字符串向量在 hdf5r 中自动处理为变长字符串
    grp$create_dataset(col, val)
  }
}