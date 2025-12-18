write_matrix_h5 <- function(h5_group, mat, compression_level = 4) {
  mat_t <- Matrix::t(mat)
  h5_group$create_dataset('data', mat_t@x, gzip_level = compression_level)
  h5_group$create_dataset('indices', mat_t@i, gzip_level = compression_level)
  h5_group$create_dataset('indptr', mat_t@p, gzip_level = compression_level)
  h5_group$create_attr('encoding-type', 'csc_matrix')
  h5_group$create_attr('encoding-version', '0.1.0')
  h5_group$create_attr('shape', dim(mat_t))
}

write_dataframe_h5 <- function(h5_file, group_name, df, index_name = '_index') {
  if (is.null(df)) return()
  if (h5_file$exists(group_name)) h5_file$link_delete(group_name)
  grp <- h5_file$create_group(group_name)
  row_names <- rownames(df)
  if (is.null(row_names)) row_names <- as.character(1:nrow(df))
  grp$create_dataset(index_name, row_names)
  grp$create_attr('_index', index_name)
  grp$create_attr('encoding-type', 'dataframe')
  grp$create_attr('encoding-version', '0.2.0')
  grp$create_attr('column-order', c(index_name, colnames(df)))
  for (col in colnames(df)) {
    val <- df[[col]]
    if (is.factor(val)) val <- as.character(val)
    if (is.character(val)) val[is.na(val)] <- ''
    grp$create_dataset(col, val)
  }
}

write_obsm_h5 <- function(h5_file, seurat_obj) {
  reductions <- names(seurat_obj@reductions)
  if (length(reductions) == 0) return()
  if (h5_file$exists('obsm')) h5_file$link_delete('obsm')
  obsm_grp <- h5_file$create_group('obsm')
  for (red in reductions) {
    emb <- Seurat::Embeddings(seurat_obj[[red]])
    key_name <- paste0('X_', red)
    obsm_grp$create_dataset(key_name, t(emb))
  }
}

write_layers_h5 <- function(h5_file, seurat_obj) {
  assay_name <- Seurat::DefaultAssay(seurat_obj)
  layer_names <- names(seurat_obj[[assay_name]]@layers)
  if (length(layer_names) == 0) return()
  if (!h5_file$exists('layers')) layers_grp <- h5_file$create_group('layers') else layers_grp <- h5_file[['layers']]
  for (lname in layer_names) {
    if (lname == 'data') next
    mat <- Seurat::GetAssayData(seurat_obj, layer = lname, assay = assay_name)
    save_name <- if (lname == 'scale.data') 'scale_data' else lname
    if(layers_grp$exists(save_name)) layers_grp$link_delete(save_name)
    if (inherits(mat, 'matrix')) {
       layers_grp$create_dataset(save_name, t(mat))
    } else {
       sub_grp <- layers_grp$create_group(save_name); write_matrix_h5(sub_grp, mat)
    }
  }
}
