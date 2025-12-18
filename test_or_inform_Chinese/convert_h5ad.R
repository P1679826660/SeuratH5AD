# ==============================================================================
# 模块名称：H5AD Converter Main (V6 Full Stack)
# 功能描述：最终版调度器。处理 Raw/X 关系，组装 Layers, Reductions, Meta.features
# ==============================================================================

#' 主函数：将 H5AD 转换为 Seurat V5 Object (Full)
#' @export
h5ad_to_seurat <- function(file_path) {
  message(">>> [Phase 1] 初始化 HDF5: ", file_path)
  file_h5 <- hdf5r::H5File$new(file_path, mode = "r")
  on.exit(file_h5$close_all(), add = TRUE)
  
  # === Phase 2: 元数据读取与方向侦测 ===
  message(">>> [Phase 2] 侦测数据方向 (DNA Check)...")
  obs <- read_h5ad_dataframe(file_h5, "obs")
  var <- read_h5ad_dataframe(file_h5, "var")
  
  # 简化的 DNA 探测逻辑 (基于之前的 V5)
  obs_names <- if(!is.null(obs)) rownames(obs) else character(0)
  var_names <- if(!is.null(var)) rownames(var) else character(0)
  
  # 如果 obs 名字像基因 (非 DNA)，且 var 像 DNA -> Transpose Mode
  obs_dna <- calculate_dna_score(obs_names)
  var_dna <- calculate_dna_score(var_names)
  
  transpose_mode <- FALSE
  if (var_dna > 0.8 && obs_dna < 0.5) {
    message("!!! 判定为转置模式 (obs=Genes, var=Cells)")
    transpose_mode <- TRUE
  }
  
  # 确定最终元数据
  final_cells_meta <- if(transpose_mode) var else obs
  final_features_meta <- if(transpose_mode) obs else var
  
  # === Phase 3: 确定 Counts 来源 (Raw vs X) ===
  # 逻辑：如果有 file[['raw/X']]，那通常是 counts。此时 file[['X']] 通常是 data。
  # 如果没有 raw，file[['X']] 就是 counts。
  
  counts_node <- NULL
  data_node <- NULL # 用于存储 normalized data
  
  has_raw <- file_h5$exists("raw") && file_h5[["raw"]]$exists("X")
  
  if (has_raw) {
    message(">>> 发现 'raw' 节点，使用 raw/X 作为 Counts，根目录 X 作为 Data...")
    counts_node <- file_h5[["raw/X"]]
    # 注意：如果使用了 raw，features 元数据通常在 raw/var 中
    raw_var <- read_h5ad_dataframe(file_h5, "raw/var")
    if (!is.null(raw_var) && !transpose_mode) final_features_meta <- raw_var
    
    # 根目录 X 存为 data (如果它存在)
    if (file_h5$exists("X")) data_node <- file_h5[["X"]]
    
  } else {
    message(">>> 未发现 'raw'，使用根目录 X 作为 Counts...")
    if (file_h5$exists("X")) counts_node <- file_h5[["X"]]
  }
  
  if (is.null(counts_node)) stop("无法找到表达矩阵 (X 或 raw/X)")
  
  # === Phase 4: 读取主矩阵并构建对象 ===
  message(">>> [Phase 4] 读取 Counts 矩阵...")
  # 这里的 transpose_output 逻辑：
  # 如果 transpose_mode=TRUE (obs是基因)，意味着文件里是 (Genes x Cells)。
  # read_matrix_node 默认读 CSR 会转置。如果读 CSR -> 变 (Cells x Genes)。
  # 我们需要 (Genes x Cells)。所以需要 transpose_output = TRUE。
  # 这里的逻辑非常绕，最好的办法是：读出来，检查维度，不对就转。
  
  counts_mat <- read_matrix_node(counts_node, transpose_output = FALSE)
  
  # 强制维度对齐
  target_cells <- nrow(final_cells_meta)
  target_feats <- nrow(final_features_meta)
  
  if (ncol(counts_mat) != target_cells) {
    counts_mat <- Matrix::t(counts_mat)
  }
  
  # 赋予名字
  if(!is.null(final_cells_meta)) colnames(counts_mat) <- rownames(final_cells_meta)
  if(!is.null(final_features_meta)) rownames(counts_mat) <- rownames(final_features_meta)
  
  message(">>> 构建基础 Seurat 对象...")
  seu <- Seurat::CreateSeuratObject(counts = counts_mat, meta.data = final_cells_meta, project = "H5AD")
  
  # === Phase 5: 注入 Data (如果有) ===
  if (!is.null(data_node)) {
    message(">>> [Phase 5] 注入 Normalized Data...")
    data_mat <- read_matrix_node(data_node, transpose_output = FALSE)
    if (ncol(data_mat) != target_cells) data_mat <- Matrix::t(data_mat)
    
    # 赋予名字以确保安全
    colnames(data_mat) <- colnames(seu)
    rownames(data_mat) <- rownames(seu)
    
    seu <- Seurat::SetAssayData(seu, slot = "data", new.data = data_mat)
  }
  
  # === Phase 6: 注入 Feature Metadata ===
  if (!is.null(final_features_meta)) {
    message(">>> [Phase 6] 注入 Feature Metadata...")
    default_assay <- Seurat::DefaultAssay(seu)
    for (i in colnames(final_features_meta)) {
      seu[[default_assay]][[i]] <- final_features_meta[[i]]
    }
    
    # 检测高变基因
    if ("highly_variable" %in% colnames(final_features_meta)) {
      hv_genes <- rownames(final_features_meta)[which(final_features_meta$highly_variable == TRUE)]
      if (length(hv_genes) > 0) {
        Seurat::VariableFeatures(seu) <- hv_genes
        message(paste("    已标记高变基因数量:", length(hv_genes)))
      }
    }
  }
  
  # === Phase 7: 注入其他 Layers, Reductions ===
  # 必须先 source utils_h5_components.R
  if (exists("add_reductions")) {
    message(">>> [Phase 7] 注入 Reductions & Layers...")
    seu <- add_layers(seu, file_h5, transpose_mode)
    seu <- add_reductions(seu, file_h5, transpose_mode)
  }
  
  message(">>> 全部转换完成！")
  return(seu)
}