read_h5ad_dataframe <- function(h5_file, group_name = 'obs') {
  if (!h5_file$exists(group_name)) return(NULL)

  grp <- h5_file[[group_name]]
  attrs <- h5_node_attributes(grp)
  index_col_name <- if ('_index' %in% names(attrs)) attrs[['_index']] else '_index'

  df_list <- list()
  
  # 获取该组下所有数据集的名称
  dataset_names <- grp$ls(recursive = FALSE)
  dataset_names <- dataset_names[dataset_names$link.type == 'H5L_TYPE_HARD', 'name']
  # 排除以双下划线开头的（通常是内部元数据）
  dataset_names <- dataset_names[!grepl('^__', dataset_names)]

  row_names_vals <- NULL

  for (name in dataset_names) {
    obj <- grp[[name]]
    
    # 确保是数据集
    if (inherits(obj, 'H5D')) {
      val <- tryCatch(obj$read(), error = function(e) NULL)
      if (is.null(val)) next

      # --- [FIX START] Categorical/Factor 自动解码逻辑 ---
      # 检查该列是否有 'categories' 属性
      obj_attrs <- h5_node_attributes(obj)
      if ('categories' %in% names(obj_attrs)) {
        cat_ref <- obj_attrs[['categories']]
        
        # AnnData 通常将 categories 存储为对象引用 (Object Reference)
        # 我们需要解析这个引用来获取真正的文本标签
        categories <- NULL
        tryCatch({
          if (inherits(cat_ref, "H5Ref")) {
             # 通过引用找到对应的 HDF5 对象（通常在 /obs/__categories/xxx 或类似位置）
             cat_dset <- h5_file[[cat_ref]]
             if (inherits(cat_dset, 'H5D')) {
               categories <- cat_dset$read()
             }
          }
        }, error = function(e) {
          # 仅在调试时输出，生产环境保持静默或记录日志
          message("Warning: Failed to resolve categories for ", name)
        })

        if (!is.null(categories)) {
          # 执行映射：Python Index (0-based) -> R Index (1-based) -> Labels
          # 1. 处理数值类型，确保是整数
          val_indices <- as.integer(val)
          
          # 2. Python 中 -1 通常代表 NA/NaN，或者单纯的 NaN
          # 将 Python 的 0-based 索引转换为 R 的 1-based
          val_r_indices <- val_indices + 1
          
          # 3. 映射值。任何超出范围的索引（如原本是 -1 变成 0）都会自动变成 NA
          # 只有有效的索引能取到值
          val_mapped <- categories[val_r_indices]
          
          # 4. 如果映射成功，覆盖原始数值
          val <- val_mapped
          
          # 可选：如果希望在 R 里也是 factor 类型，取消下面注释
          # val <- factor(val, levels = categories)
        }
      }
      # --- [FIX END] ---

      if (name == index_col_name) {
        row_names_vals <- val
      } else if (is.null(dim(val)) || length(dim(val)) == 1) {
        df_list[[name]] <- val
      }
    }
  }

  # 组装 DataFrame (逻辑保持不变)
  if (length(df_list) == 0 && is.null(row_names_vals)) return(NULL)

  df <- tryCatch({
    if (length(df_list) > 0) {
      max_len <- if (!is.null(row_names_vals)) length(row_names_vals) else max(sapply(df_list, length))

      safe_list <- lapply(df_list, function(v) {
        if (length(v) == max_len) return(v)
        if (length(v) < max_len) return(c(v, rep(NA, max_len - length(v))))
        return(v[1:max_len])
      })
      data.frame(safe_list, stringsAsFactors = FALSE)
    } else {
      data.frame(row.names = seq_along(row_names_vals))
    }
  }, error = function(e) {
    message("Warning: Failed to construct dataframe for ", group_name)
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)

  if (!is.null(row_names_vals)) {
    row_names_vals <- as.character(row_names_vals)
    if (anyDuplicated(row_names_vals)) {
      row_names_vals <- make.unique(row_names_vals)
    }
    if (length(row_names_vals) == nrow(df)) {
      rownames(df) <- row_names_vals
    }
  }
  return(df)
}
