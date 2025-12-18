read_h5ad_dataframe <- function(h5_file, group_name = 'obs') {
  if (!h5_file$exists(group_name)) return(NULL)
  grp <- h5_file[[group_name]]; attrs <- h5_node_attributes(grp)
  index_col_name <- if('_index' %in% names(attrs)) attrs[['_index']] else '_index'
  df_list <- list(); dataset_names <- grp$ls(recursive = FALSE)$name
  dataset_names <- dataset_names[!grepl('^__', dataset_names)]
  row_names_vals <- NULL
  for (name in dataset_names) {
    obj <- grp[[name]]
    if (inherits(obj, 'H5D')) {
      val <- obj$read()
      if (name == index_col_name) row_names_vals <- val
      else if (is.null(dim(val)) || length(dim(val)) == 1) df_list[[name]] <- val
    }
  }
  if (length(df_list) == 0) {
    if (!is.null(row_names_vals)) df <- data.frame(row.names = row_names_vals) else return(NULL)
  } else {
    df <- data.frame(df_list, stringsAsFactors = FALSE)
    if (!is.null(row_names_vals)) {
       if (anyDuplicated(row_names_vals)) row_names_vals <- make.unique(as.character(row_names_vals))
       if (length(row_names_vals) == nrow(df)) rownames(df) <- row_names_vals
    }
  }
  return(df)
}

add_reductions <- function(seu, h5_file, transpose_mode) {
  emb_group_name <- if(transpose_mode) 'varm' else 'obsm'
  cell_names <- colnames(seu)
  if (h5_file$exists(emb_group_name)) {
    grp <- h5_file[[emb_group_name]]
    for (k in grp$ls(recursive=FALSE)$name) {
      if (!inherits(grp[[k]], 'H5D') && !inherits(grp[[k]], 'H5Group')) next
      raw_emb <- tryCatch(grp[[k]]$read(), error=function(e) NULL)
      if (is.null(raw_emb)) next
      if (!is.matrix(raw_emb)) raw_emb <- as.matrix(raw_emb)
      if (ncol(raw_emb) == length(cell_names)) raw_emb <- t(raw_emb)
      if (nrow(raw_emb) != length(cell_names)) next
      clean_name <- tolower(gsub('^X_', '', k))
      seurat_key <- paste0(clean_name, '_')
      rownames(raw_emb) <- cell_names
      colnames(raw_emb) <- paste0(seurat_key, 1:ncol(raw_emb))
      dr_obj <- Seurat::CreateDimReducObject(embeddings = raw_emb, key = seurat_key, assay = Seurat::DefaultAssay(seu))
      seu[[clean_name]] <- dr_obj
    }
  }
  return(seu)
}

add_layers <- function(seu, h5_file, transpose_mode) {
  if (!h5_file$exists('layers')) return(seu)
  grp <- h5_file[['layers']]
  for (lname in grp$ls(recursive=FALSE)$name) {
    mat <- read_matrix_node(grp[[lname]], transpose_output = transpose_mode)
    if (is.null(mat)) next
    if (ncol(mat) != ncol(seu)) mat <- Matrix::t(mat)
    if (ncol(mat) == ncol(seu) && nrow(mat) == nrow(seu)) {
      target_slot <- if(lname %in% c('logcounts','norm_data')) 'data' else if(grepl('scale', lname)) 'scale.data' else lname
      da <- Seurat::DefaultAssay(seu)
      if (target_slot == 'scale.data') {
        Seurat::SetAssayData(seu, slot = 'scale.data', new.data = as.matrix(mat), assay = da)
      } else if (target_slot == 'data') {
        Seurat::SetAssayData(seu, slot = 'data', new.data = mat, assay = da)
      } else {
        Seurat::LayerData(seu, assay = da, layer = lname) <- mat
      }
    }
  }
  return(seu)
}
