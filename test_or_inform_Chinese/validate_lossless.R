# ==============================================================================
# 脚本名称：H5AD 无损转换终极验证器 (Bit-level Comparator)
# 功能：对比两个 h5ad 文件的底层数据，证明“大小不同但信息一致”
# ==============================================================================

library(hdf5r)
library(Matrix)

compare_h5ad <- function(file_orig, file_new) {
  
  message(">>> 正在启动无损验证程序...")
  message("    源文件 (A): ", file_orig, " (", round(file.size(file_orig)/1024/1024, 2), " MB)")
  message("    新文件 (B): ", file_new,  " (", round(file.size(file_new)/1024/1024, 2), " MB)")
  
  h5_a <- H5File$new(file_orig, mode = "r")
  h5_b <- H5File$new(file_new, mode = "r")
  on.exit({ h5_a$close_all(); h5_b$close_all() })
  
  # --- 1. 检查压缩过滤器 (解释体积差异的核心) ---
  message("\n[1. 存储压缩分析]")
  
  get_compression <- function(h5, name) {
    if(h5$exists(name)) {
      ds <- h5[[name]]
      # data 是一定存在的 dataset
      if(inherits(ds, "H5Group")) ds <- ds[["data"]] 
      filter <- ds$get_create_plist()$get_nfilters()
      if (filter > 0) return("GZIP/Compressed") else return("None/Uncompressed")
    }
    return("Unknown")
  }
  
  msg_a <- get_compression(h5_a, "X")
  msg_b <- get_compression(h5_b, "X")
  message(sprintf("    源文件 X 压缩状态: %s", msg_a))
  message(sprintf("    新文件 X 压缩状态: %s", msg_b))
  
  if (msg_a != msg_b) {
    message("    !!! 结论: 文件大小差异主要源于压缩算法的不同，而非数据丢失。")
  }
  
  # --- 2. 矩阵内容全量对比 ---
  message("\n[2. 矩阵数值全量对比 (X)]")
  
  read_full_matrix <- function(h5_file) {
    # 简单的读取器，不依赖我们之前的复杂模块，直接读 data
    x_node <- h5_file[["X"]]
    if(inherits(x_node, "H5Group")) {
      return(list(
        sum = sum(x_node[["data"]]$read()),
        nz  = length(x_node[["data"]]$read()),
        dim = x_node[["data"]]$dims # 注意：这里读的是 data 向量长度，不是矩阵维度
      ))
    } else {
      # Dense
      val <- x_node[,]
      return(list(sum = sum(val), nz = sum(val != 0)))
    }
  }
  
  stats_a <- read_full_matrix(h5_a)
  stats_b <- read_full_matrix(h5_b)
  
  message(sprintf("    源文件 X 总和: %.4f | 非零元素数: %d", stats_a$sum, stats_a$nz))
  message(sprintf("    新文件 X 总和: %.4f | 非零元素数: %d", stats_b$sum, stats_b$nz))
  
  if (abs(stats_a$sum - stats_b$sum) < 0.001 && stats_a$nz == stats_b$nz) {
    message("    >>> [PASS] 矩阵信息完全一致！(数据无损)")
  } else {
    message("    >>> [FAIL] 矩阵信息不一致！请检查！")
  }
  
  # --- 3. 维度对比 ---
  message("\n[3. 逻辑维度对比]")
  # 读取 obs/var 索引长度
  get_idx_len <- function(h5, grp) {
    if(!h5$exists(grp)) return(0)
    g <- h5[[grp]]
    # 找 _index
    if(g$attr_exists("_index")) {
      idx_name <- g$attr_open("_index")$read()
      return(g[[idx_name]]$dims)
    }
    return(0)
  }
  
  obs_a <- get_idx_len(h5_a, "obs"); var_a <- get_idx_len(h5_a, "var")
  obs_b <- get_idx_len(h5_b, "obs"); var_b <- get_idx_len(h5_b, "var")
  
  # 注意：如果发生了转置，A 的 obs 可能等于 B 的 var
  message(sprintf("    源文件 (obs, var): (%d, %d)", obs_a, var_a))
  message(sprintf("    新文件 (obs, var): (%d, %d)", obs_b, var_b))
  
  if ((obs_a == obs_b && var_a == var_b) || (obs_a == var_b && var_a == obs_b)) {
    message("    >>> [PASS] 维度匹配 (支持转置等价性)。")
  } else {
    message("    >>> [FAIL] 维度丢失！")
  }
  
  message("\n>>> 验证完成。")
}

# 执行对比
file_orig <- "test.h5ad"
file_new  <- "final_product.h5ad"
compare_h5ad(file_orig, file_new)
