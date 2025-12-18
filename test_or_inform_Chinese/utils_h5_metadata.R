# ==============================================================================
# 模块名称：H5AD Metadata Reader (V2 Fixed)
# 功能描述：读取 obs 或 var 组并转换为 DataFrame，增强了对空列情况的处理
# ==============================================================================

#' 读取 obs 或 var 组并转换为 DataFrame
read_h5ad_dataframe <- function(h5_file, group_name = "obs") {
  if (!h5_file$exists(group_name)) return(NULL)
  
  grp <- h5_file[[group_name]]
  
  # 1. 寻找索引 (Index)
  # AnnData 通常在 group 属性中存储 '_index' 指向索引列名
  attrs <- h5_node_attributes(grp)
  index_col_name <- if("_index" %in% names(attrs)) attrs[["_index"]] else "_index"
  
  # 2. 遍历读取列
  df_list <- list()
  dataset_names <- grp$ls(recursive = FALSE)$name
  
  # 过滤掉特殊的 __categories (如果有)
  dataset_names <- dataset_names[!grepl("^__", dataset_names)]
  
  row_names_vals <- NULL
  
  for (name in dataset_names) {
    obj <- grp[[name]]
    
    # 只处理 Dataset
    if (inherits(obj, "H5D")) {
      val <- obj$read()
      
      if (name == index_col_name) {
        row_names_vals <- val
      } else {
        # 确保读取的是向量，如果是多维数组(除了1D)，可能需要处理
        if (is.null(dim(val)) || length(dim(val)) == 1) {
          df_list[[name]] <- val
        } else {
          # 如果遇到多维数组作为metadata列(比较少见)，暂且跳过或取第一列
          # 为了工程稳定性，这里选择警告并跳过，防止破坏 data.frame 结构
          warning(paste("跳过复杂结构的元数据列:", name))
        }
      }
    }
  }
  
  # 3. 构建 DataFrame (核心修复部分)
  
  # 情况 A: 只有索引，没有其他数据列 (df_list 为空)
  if (length(df_list) == 0) {
    if (!is.null(row_names_vals)) {
      # 创建一个只有行名，没有列的 DataFrame
      # 必须显式指定 row.names 来确立行数
      df <- data.frame(row.names = row_names_vals)
    } else {
      # 既没数据也没索引，返回 NULL 或空对象
      return(NULL)
    }
  } else {
    # 情况 B: 有数据列
    # 先构建数据框，暂不设行名，避免长度不匹配直接报错
    df <- data.frame(df_list, stringsAsFactors = FALSE)
    
    # 安全地设置行名
    if (!is.null(row_names_vals)) {
      # 1. 检查去重
      if (anyDuplicated(row_names_vals)) {
        warning(paste0(group_name, " 索引包含重复项，正在通过 make.unique 去重..."))
        row_names_vals <- make.unique(as.character(row_names_vals))
      }
      
      # 2. 检查长度一致性 (这是你报错的核心原因)
      if (length(row_names_vals) == nrow(df)) {
        rownames(df) <- row_names_vals
      } else {
        # 严重警告：索引长度与数据长度不符
        warning(sprintf(
          "警告：在 %s 中，索引长度 (%d) 与数据列长度 (%d) 不匹配。将忽略原有索引，使用数字索引。",
          group_name, length(row_names_vals), nrow(df)
        ))
        # 此时保留数字索引，不强行赋值，防止 crash
      }
    }
  }
  
  return(df)
}

# 确保辅助函数存在 (如果你的 convert_h5ad.R 没有独立 source 这个，建议保留在这里或确保 matrix 模块已加载)
# 这里为了安全起见，如果不报错重复定义，可以不写。
# 但为了确保独立性，假设它在 utils_h5_matrix.R 中定义了，这里就不重复定义了。