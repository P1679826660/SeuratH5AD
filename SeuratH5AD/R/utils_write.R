#' @importFrom Matrix t
#' @importFrom methods as
NULL

write_matrix_h5 <- function(h5_group, mat, compression_level = 4) {
  # Ensure input is a compressed sparse column matrix (dgCMatrix)
  # This guarantees the existence of @i, @p, and @x slots
  if (!inherits(mat, "CsparseMatrix")) {
    mat <- methods::as(mat, "CsparseMatrix")
  }

  # H5AD 'X' is typically stored as (n_obs x n_vars) -> (Cells x Genes)
  # Seurat stores data as (n_features x n_cells) -> (Genes x Cells)
  # We transpose here to match H5AD standard
  mat_t <- Matrix::t(mat)

  # Ensure the transposed matrix is also sparse
  if (!inherits(mat_t, "CsparseMatrix")) {
    mat_t <- methods::as(mat_t, "CsparseMatrix")
  }

  h5_group$create_dataset('data', mat_t@x, gzip_level = compression_level)
  h5_group$create_dataset('indices', mat_t@i, gzip_level = compression_level)
  h5_group$create_dataset('indptr', mat_t@p, gzip_level = compression_level)

  # Explicitly define attributes for CSC matrix
  # While CSR is common in Python, CSC is valid and native to R's internal representation
  h5_group$create_attr('encoding-type', 'csc_matrix')
  h5_group$create_attr('encoding-version', '0.1.0')
  h5_group$create_attr('shape', dim(mat_t))
}

write_dataframe_h5 <- function(h5_file, group_name, df, index_name = '_index') {
  if (is.null(df)) return()
  if (!is.data.frame(df)) df <- as.data.frame(df)

  if (h5_file$exists(group_name)) h5_file$link_delete(group_name)
  grp <- h5_file$create_group(group_name)

  row_names <- rownames(df)
  if (is.null(row_names)) row_names <- as.character(seq_len(nrow(df)))
  
  grp$create_dataset(index_name, row_names)
  grp$create_attr('_index', index_name)
  grp$create_attr('encoding-type', 'dataframe')
  grp$create_attr('encoding-version', '0.2.0')

  valid_cols <- c(index_name)

  for (col in colnames(df)) {
    val <- df[[col]]

    # Flatten list columns to avoid HDF5 write errors
    if (is.list(val)) {
      # Collapse list elements into semi-colon separated strings
      val <- vapply(val, function(x) paste(as.character(x), collapse = ";"), character(1))
    }

    if (is.factor(val)) val <- as.character(val)
    if (is.character(val)) val[is.na(val)] <- ''

    # Only write atomic vectors
    if (is.atomic(val)) {
      # Convert logicals to integers as HDF5 does not support boolean natively in some contexts
      if (is.logical(val)) val <- as.integer(val)

      tryCatch({
        grp$create_dataset(col, val)
        valid_cols <- c(valid_cols, col)
      }, error = function(e) {
        warning(sprintf("Skipping column '%s' in %s: %s", col, group_name, e$message))
      })
    }
  }
  grp$create_attr('column-order', valid_cols)
}

write_obsm_h5 <- function(h5_file, seurat_obj) {
  reductions <- names(seurat_obj@reductions)
  if (length(reductions) == 0) return()

  if (h5_file$exists('obsm')) h5_file$link_delete('obsm')
  obsm_grp <- h5_file$create_group('obsm')

  for (red in reductions) {
    emb <- tryCatch({
      Seurat::Embeddings(seurat_obj[[red]])
    }, error = function(e) NULL)

    if (is.null(emb)) next

    # Standardize naming convention to X_name
    key_name <- paste0('X_', red)

    if (!is.matrix(emb)) emb <- as.matrix(emb)

    # Embeddings are (Cells x PCs), which matches H5AD requirements. No transpose needed.
    obsm_grp$create_dataset(key_name, emb)
  }
}

write_layers_h5 <- function(h5_file, seurat_obj) {
  assay_name <- Seurat::DefaultAssay(seurat_obj)
  
  # Support Seurat V5 layers
  layer_names <- tryCatch({
    SeuratObject::Layers(seurat_obj, assay = assay_name)
  }, error = function(e) {
    # Fallback for older Seurat objects or specific assay types
    names(seurat_obj[[assay_name]]@layers)
  })

  if (length(layer_names) == 0) return()

  if (!h5_file$exists('layers')) layers_grp <- h5_file$create_group('layers') else layers_grp <- h5_file[['layers']]

  for (lname in layer_names) {
    # Skip standard slots as they are handled in X or raw
    if (lname %in% c('data', 'counts')) next

    mat <- tryCatch({
      Seurat::GetAssayData(seurat_obj, layer = lname, assay = assay_name)
    }, error = function(e) NULL)

    if (is.null(mat)) next

    save_name <- if (lname == 'scale.data') 'scale_data' else lname
    if (layers_grp$exists(save_name)) layers_grp$link_delete(save_name)

    sub_grp <- layers_grp$create_group(save_name)
    write_matrix_h5(sub_grp, mat)
  }
}
