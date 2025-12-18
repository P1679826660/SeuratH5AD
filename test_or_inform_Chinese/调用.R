# ==============================================================================
# SeuratH5AD 终极修复版构建脚本 (Windows Safe Mode)
# ==============================================================================

# 1. 清理环境
pkg_name <- "SeuratH5AD"
if (dir.exists(pkg_name)) unlink(pkg_name, recursive = TRUE, force = TRUE)
dir.create(pkg_name)
dir.create(file.path(pkg_name, "R"))
message(">>> [1/6] 目录已重置: ", pkg_name)

# 2. 写入 DESCRIPTION
desc <- c(
  paste0("Package: ", pkg_name),
  "Title: High-Performance H5AD and Seurat V5 Converter",
  "Version: 0.1.0",
  "Authors@R: person(\"Your\", \"Name\", , \"your@email.com\", role = c(\"aut\", \"cre\"))",
  "Description: A lightweight R package for converting .h5ad to Seurat V5.",
  "License: MIT",
  "Encoding: UTF-8",
  "LazyData: true",
  "Roxygen: list(markdown = TRUE)",
  "RoxygenNote: 7.2.3",
  "Imports:",
  "    hdf5r (>= 1.3.8),",
  "    Seurat (>= 5.0.0),",
  "    SeuratObject,",
  "    Matrix,",
  "    methods"
)
writeLines(desc, file.path(pkg_name, "DESCRIPTION"))

# 3. 写入 utils_core.R (无中文，纯净版)
code_core <- c(
  "#' @importFrom Matrix sparseMatrix t",
  "#' @importFrom hdf5r h5attr_names h5attr",
  "NULL",
  "",
  "read_matrix_node <- function(node, transpose_output = FALSE) {",
  "  mat <- NULL",
  "  if (inherits(node, 'H5Group')) {",
  "    if(!exists('h5attr_names', where = asNamespace('hdf5r'))) return(NULL)",
  "    attrs <- hdf5r::h5attr_names(node)",
  "    encoding_attr <- if('encoding-type' %in% attrs) hdf5r::h5attr(node, 'encoding-type') else 'csr_matrix'",
  "    shape_attr <- if('shape' %in% attrs) hdf5r::h5attr(node, 'shape') else NULL",
  "    data <- node[['data']]$read()",
  "    indices <- node[['indices']]$read()",
  "    indptr <- node[['indptr']]$read()",
  "    if (encoding_attr == 'csr_matrix') {",
  "      mat <- Matrix::sparseMatrix(i = indices + 1, p = indptr, x = data, dims = c(shape_attr[2], shape_attr[1]), index1 = TRUE)",
  "      if (!transpose_output) mat <- Matrix::t(mat)",
  "    } else if (encoding_attr == 'csc_matrix') {",
  "      mat <- Matrix::sparseMatrix(i = indices + 1, p = indptr, x = data, dims = c(shape_attr[1], shape_attr[2]), index1 = TRUE)",
  "      if (transpose_output) mat <- Matrix::t(mat)",
  "    }",
  "  } else {",
  "    raw_data <- node[,]",
  "    mat <- methods::as(raw_data, 'CsparseMatrix')",
  "    if (transpose_output) mat <- Matrix::t(mat)",
  "  }",
  "  return(mat)",
  "}",
  "",
  "calculate_dna_score <- function(identifiers) {",
  "  if (length(identifiers) == 0) return(0)",
  "  n_sample <- min(length(identifiers), 100)",
  "  sample_ids <- identifiers[1:n_sample]",
  "  clean_ids <- toupper(gsub('[^a-zA-Z]', '', sample_ids))",
  "  total_chars <- sum(nchar(clean_ids))",
  "  if (total_chars == 0) return(0)",
  "  non_acgt_chars <- gsub('[ACGT]', '', clean_ids)",
  "  return(1 - (sum(nchar(non_acgt_chars)) / total_chars))",
  "}",
  "",
  "h5_node_attributes <- function(node) {",
  "  n <- hdf5r::h5attr_names(node)",
  "  res <- list(); for (i in n) res[[i]] <- hdf5r::h5attr(node, i); return(res)",
  "}"
)
writeLines(code_core, file.path(pkg_name, "R/utils_core.R"))
message(">>> [2/6] 核心工具写入完成")

# 4. 写入 utils_read.R
code_read <- c(
  "read_h5ad_dataframe <- function(h5_file, group_name = 'obs') {",
  "  if (!h5_file$exists(group_name)) return(NULL)",
  "  grp <- h5_file[[group_name]]; attrs <- h5_node_attributes(grp)",
  "  index_col_name <- if('_index' %in% names(attrs)) attrs[['_index']] else '_index'",
  "  df_list <- list(); dataset_names <- grp$ls(recursive = FALSE)$name",
  "  dataset_names <- dataset_names[!grepl('^__', dataset_names)]",
  "  row_names_vals <- NULL",
  "  for (name in dataset_names) {",
  "    obj <- grp[[name]]",
  "    if (inherits(obj, 'H5D')) {",
  "      val <- obj$read()",
  "      if (name == index_col_name) row_names_vals <- val",
  "      else if (is.null(dim(val)) || length(dim(val)) == 1) df_list[[name]] <- val",
  "    }",
  "  }",
  "  if (length(df_list) == 0) {",
  "    if (!is.null(row_names_vals)) df <- data.frame(row.names = row_names_vals) else return(NULL)",
  "  } else {",
  "    df <- data.frame(df_list, stringsAsFactors = FALSE)",
  "    if (!is.null(row_names_vals)) {",
  "       if (anyDuplicated(row_names_vals)) row_names_vals <- make.unique(as.character(row_names_vals))",
  "       if (length(row_names_vals) == nrow(df)) rownames(df) <- row_names_vals",
  "    }",
  "  }",
  "  return(df)",
  "}",
  "",
  "add_reductions <- function(seu, h5_file, transpose_mode) {",
  "  emb_group_name <- if(transpose_mode) 'varm' else 'obsm'",
  "  cell_names <- colnames(seu)",
  "  if (h5_file$exists(emb_group_name)) {",
  "    grp <- h5_file[[emb_group_name]]",
  "    for (k in grp$ls(recursive=FALSE)$name) {",
  "      if (!inherits(grp[[k]], 'H5D') && !inherits(grp[[k]], 'H5Group')) next",
  "      raw_emb <- tryCatch(grp[[k]]$read(), error=function(e) NULL)",
  "      if (is.null(raw_emb)) next",
  "      if (!is.matrix(raw_emb)) raw_emb <- as.matrix(raw_emb)",
  "      if (ncol(raw_emb) == length(cell_names)) raw_emb <- t(raw_emb)",
  "      if (nrow(raw_emb) != length(cell_names)) next",
  "      clean_name <- tolower(gsub('^X_', '', k))",
  "      seurat_key <- paste0(clean_name, '_')",
  "      rownames(raw_emb) <- cell_names",
  "      colnames(raw_emb) <- paste0(seurat_key, 1:ncol(raw_emb))",
  "      dr_obj <- Seurat::CreateDimReducObject(embeddings = raw_emb, key = seurat_key, assay = Seurat::DefaultAssay(seu))",
  "      seu[[clean_name]] <- dr_obj",
  "    }",
  "  }",
  "  return(seu)",
  "}",
  "",
  "add_layers <- function(seu, h5_file, transpose_mode) {",
  "  if (!h5_file$exists('layers')) return(seu)",
  "  grp <- h5_file[['layers']]",
  "  for (lname in grp$ls(recursive=FALSE)$name) {",
  "    mat <- read_matrix_node(grp[[lname]], transpose_output = transpose_mode)",
  "    if (is.null(mat)) next",
  "    if (ncol(mat) != ncol(seu)) mat <- Matrix::t(mat)",
  "    if (ncol(mat) == ncol(seu) && nrow(mat) == nrow(seu)) {",
  "      target_slot <- if(lname %in% c('logcounts','norm_data')) 'data' else if(grepl('scale', lname)) 'scale.data' else lname",
  "      da <- Seurat::DefaultAssay(seu)",
  "      if (target_slot == 'scale.data') {",
  "        Seurat::SetAssayData(seu, slot = 'scale.data', new.data = as.matrix(mat), assay = da)",
  "      } else if (target_slot == 'data') {",
  "        Seurat::SetAssayData(seu, slot = 'data', new.data = mat, assay = da)",
  "      } else {",
  "        Seurat::LayerData(seu, assay = da, layer = lname) <- mat",
  "      }",
  "    }",
  "  }",
  "  return(seu)",
  "}"
)
writeLines(code_read, file.path(pkg_name, "R/utils_read.R"))
message(">>> [3/6] 读取模块写入完成")

# 5. 写入 utils_write.R
code_write <- c(
  "write_matrix_h5 <- function(h5_group, mat, compression_level = 4) {",
  "  mat_t <- Matrix::t(mat)",
  "  h5_group$create_dataset('data', mat_t@x, gzip_level = compression_level)",
  "  h5_group$create_dataset('indices', mat_t@i, gzip_level = compression_level)",
  "  h5_group$create_dataset('indptr', mat_t@p, gzip_level = compression_level)",
  "  h5_group$create_attr('encoding-type', 'csc_matrix')",
  "  h5_group$create_attr('encoding-version', '0.1.0')",
  "  h5_group$create_attr('shape', dim(mat_t))",
  "}",
  "",
  "write_dataframe_h5 <- function(h5_file, group_name, df, index_name = '_index') {",
  "  if (is.null(df)) return()",
  "  if (h5_file$exists(group_name)) h5_file$link_delete(group_name)",
  "  grp <- h5_file$create_group(group_name)",
  "  row_names <- rownames(df)",
  "  if (is.null(row_names)) row_names <- as.character(1:nrow(df))",
  "  grp$create_dataset(index_name, row_names)",
  "  grp$create_attr('_index', index_name)",
  "  grp$create_attr('encoding-type', 'dataframe')",
  "  grp$create_attr('encoding-version', '0.2.0')",
  "  grp$create_attr('column-order', c(index_name, colnames(df)))",
  "  for (col in colnames(df)) {",
  "    val <- df[[col]]",
  "    if (is.factor(val)) val <- as.character(val)",
  "    if (is.character(val)) val[is.na(val)] <- ''",
  "    grp$create_dataset(col, val)",
  "  }",
  "}",
  "",
  "write_obsm_h5 <- function(h5_file, seurat_obj) {",
  "  reductions <- names(seurat_obj@reductions)",
  "  if (length(reductions) == 0) return()",
  "  if (h5_file$exists('obsm')) h5_file$link_delete('obsm')",
  "  obsm_grp <- h5_file$create_group('obsm')",
  "  for (red in reductions) {",
  "    emb <- Seurat::Embeddings(seurat_obj[[red]])",
  "    key_name <- paste0('X_', red)",
  "    obsm_grp$create_dataset(key_name, t(emb))",
  "  }",
  "}",
  "",
  "write_layers_h5 <- function(h5_file, seurat_obj) {",
  "  assay_name <- Seurat::DefaultAssay(seurat_obj)",
  "  layer_names <- names(seurat_obj[[assay_name]]@layers)",
  "  if (length(layer_names) == 0) return()",
  "  if (!h5_file$exists('layers')) layers_grp <- h5_file$create_group('layers') else layers_grp <- h5_file[['layers']]",
  "  for (lname in layer_names) {",
  "    if (lname == 'data') next",
  "    mat <- Seurat::GetAssayData(seurat_obj, layer = lname, assay = assay_name)",
  "    save_name <- if (lname == 'scale.data') 'scale_data' else lname",
  "    if(layers_grp$exists(save_name)) layers_grp$link_delete(save_name)",
  "    if (inherits(mat, 'matrix')) {",
  "       layers_grp$create_dataset(save_name, t(mat))",
  "    } else {",
  "       sub_grp <- layers_grp$create_group(save_name); write_matrix_h5(sub_grp, mat)",
  "    }",
  "  }",
  "}"
)
writeLines(code_write, file.path(pkg_name, "R/utils_write.R"))
message(">>> [4/6] 写入模块写入完成")

# 6. 写入 convert_h5ad.R (Import function)
code_convert <- c(
  "#' Convert H5AD to Seurat V5",
  "#' @param file_path Path to .h5ad file",
  "#' @return Seurat object",
  "#' @export",
  "h5ad_to_seurat <- function(file_path) {",
  "  message('>>> [Import] Loading: ', file_path)",
  "  file_h5 <- hdf5r::H5File$new(file_path, mode = 'r')",
  "  on.exit(file_h5$close_all(), add = TRUE)",
  "  obs <- read_h5ad_dataframe(file_h5, 'obs')",
  "  var <- read_h5ad_dataframe(file_h5, 'var')",
  "  obs_dna <- calculate_dna_score(if(!is.null(obs)) rownames(obs) else character(0))",
  "  var_dna <- calculate_dna_score(if(!is.null(var)) rownames(var) else character(0))",
  "  transpose_mode <- FALSE",
  "  if (var_dna > 0.8 && obs_dna < 0.5) {",
  "    message('!!! Transpose detected (obs=Genes, var=Cells)')",
  "    transpose_mode <- TRUE",
  "  }",
  "  final_cells_meta <- if(transpose_mode) var else obs",
  "  final_features_meta <- if(transpose_mode) obs else var",
  "  counts_node <- NULL; data_node <- NULL",
  "  has_raw <- file_h5$exists('raw') && file_h5[['raw']]$exists('X')",
  "  if (has_raw) {",
  "    counts_node <- file_h5[['raw/X']]",
  "    raw_var <- read_h5ad_dataframe(file_h5, 'raw/var')",
  "    if (!is.null(raw_var) && !transpose_mode) final_features_meta <- raw_var",
  "    if (file_h5$exists('X')) data_node <- file_h5[['X']]",
  "  } else { if (file_h5$exists('X')) counts_node <- file_h5[['X']] }",
  "  if (is.null(counts_node)) stop('No expression matrix found (X or raw/X)')",
  "  counts_mat <- read_matrix_node(counts_node, transpose_output = FALSE)",
  "  target_cells <- nrow(final_cells_meta)",
  "  if (ncol(counts_mat) != target_cells) counts_mat <- Matrix::t(counts_mat)",
  "  if(!is.null(final_cells_meta)) colnames(counts_mat) <- rownames(final_cells_meta)",
  "  if(!is.null(final_features_meta)) rownames(counts_mat) <- rownames(final_features_meta)",
  "  seu <- Seurat::CreateSeuratObject(counts = counts_mat, meta.data = final_cells_meta, project = 'H5AD')",
  "  if (!is.null(data_node)) {",
  "    data_mat <- read_matrix_node(data_node, transpose_output = FALSE)",
  "    if (ncol(data_mat) != target_cells) data_mat <- Matrix::t(data_mat)",
  "    colnames(data_mat) <- colnames(seu); rownames(data_mat) <- rownames(seu)",
  "    seu <- Seurat::SetAssayData(seu, slot = 'data', new.data = data_mat)",
  "  }",
  "  if (!is.null(final_features_meta)) {",
  "    da <- Seurat::DefaultAssay(seu)",
  "    for (i in colnames(final_features_meta)) seu[[da]][[i]] <- final_features_meta[[i]]",
  "    if ('highly_variable' %in% colnames(final_features_meta)) {",
  "      hv <- rownames(final_features_meta)[which(final_features_meta$highly_variable == TRUE)]",
  "      if(length(hv)>0) Seurat::VariableFeatures(seu) <- hv",
  "    }",
  "  }",
  "  seu <- add_layers(seu, file_h5, transpose_mode)",
  "  seu <- add_reductions(seu, file_h5, transpose_mode)",
  "  message('>>> Conversion Complete!')",
  "  return(seu)",
  "}"
)
writeLines(code_convert, file.path(pkg_name, "R/convert_h5ad.R"))
message(">>> [5/6] 转换入口写入完成")

# 7. 写入 export_h5ad.R (Export function)
code_export <- c(
  "#' Export Seurat V5 to H5AD",
  "#' @param seurat_obj Seurat Object",
  "#' @param file_path Output path",
  "#' @param compression_level 0-9",
  "#' @param write_raw_slot boolean",
  "#' @export",
  "seurat_to_h5ad <- function(seurat_obj, file_path, compression_level = 4, write_raw_slot = TRUE) {",
  "  message('>>> [Export] Start: ', file_path)",
  "  if (file.exists(file_path)) file.remove(file_path)",
  "  file_h5 <- hdf5r::H5File$new(file_path, mode = 'w')",
  "  on.exit(file_h5$close_all(), add = TRUE)",
  "  assay_name <- Seurat::DefaultAssay(seurat_obj)",
  "  layers_avail <- names(seurat_obj[[assay_name]]@layers)",
  "  main_layer <- if ('data' %in% layers_avail) 'data' else 'counts'",
  "  main_mat <- Seurat::GetAssayData(seurat_obj, layer = main_layer, assay = assay_name)",
  "  x_grp <- file_h5$create_group('X')",
  "  write_matrix_h5(x_grp, main_mat, compression_level)",
  "  write_dataframe_h5(file_h5, 'obs', seurat_obj@meta.data)",
  "  feat_meta <- seurat_obj[[assay_name]][[]]",
  "  if(nrow(feat_meta)==0) feat_meta <- data.frame(row.names=rownames(seurat_obj))",
  "  write_dataframe_h5(file_h5, 'var', feat_meta)",
  "  if (write_raw_slot && 'counts' %in% layers_avail) {",
  "    raw_grp <- file_h5$create_group('raw'); raw_x_grp <- raw_grp$create_group('X')",
  "    counts_mat <- Seurat::GetAssayData(seurat_obj, layer = 'counts', assay = assay_name)",
  "    write_matrix_h5(raw_x_grp, counts_mat, compression_level)",
  "    write_dataframe_h5(raw_grp, 'var', feat_meta)",
  "  }",
  "  write_obsm_h5(file_h5, seurat_obj)",
  "  write_layers_h5(file_h5, seurat_obj)",
  "  message('>>> Export Success!')",
  "}"
)
writeLines(code_export, file.path(pkg_name, "R/export_h5ad.R"))
message(">>> [6/6] 导出入口写入完成")

# 8. 安装
message(">>> 正在安装包...")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
current_wd <- getwd()
setwd(pkg_name)
tryCatch({
  devtools::document()
  devtools::install(upgrade = "never", force = TRUE, quiet = FALSE)
}, finally = {
  setwd(current_wd)
})

message("\n========================================================")
message(">>> 成功！请重启 R (Ctrl+Shift+F10) 后运行 library(SeuratH5AD)")
message("========================================================")




#######################
library(SeuratH5AD)

# 测试读取
seu <- h5ad_to_seurat("test.h5ad")

# 测试导出
seurat_to_h5ad(seu, "final_product.h5ad")
