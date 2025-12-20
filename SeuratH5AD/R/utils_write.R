#' @importFrom Matrix t
#' @importFrom methods as
NULL

write_matrix_h5 <- function(h5_group, mat, compression_level = 4) {
  # Robustness: Force conversion to dgCMatrix (CsparseMatrix)
  # This guarantees slots @i, @p, @x exist, preventing crashes with dense matrices.
  if (!inherits(mat, "CsparseMatrix")) {
    mat <- methods::as(mat, "CsparseMatrix")
  }

  # Transpose: Seurat (Genes x Cells) -> H5AD (Cells x Genes)
  mat_t <- Matrix::t(mat)

  # Re-verify sparseness after transpose
  if (!inherits(mat_t, "CsparseMatrix")) {
    mat_t <- methods::as(mat_t, "CsparseMatrix")
  }

  # Compatibility: Handle Logical/Boolean matrices
  # HDF5 does not natively support boolean arrays in this context. Convert to 0/1.
  data_values <- mat_t@x
  if (is.logical(data_values)) {
    data_values <- as.numeric(data_values)
  }

  h5_group$create_dataset('data', data_values, gzip_level = compression_level)
  h5_group$create_dataset('indices', mat_t@i, gzip_level = compression_level)
  h5_group$create_dataset('indptr', mat_t@p, gzip_level = compression_level)

  # Write Attributes
  h5_group$create_attr('encoding-type', 'csc_matrix')
  h5_group$create_attr('encoding-version', '0.1.0')
  h5_group$create_attr('shape', dim(mat_t))
}

write_dataframe_h5 <- function(h5_file, group_name, df, index_name = '_index') {
  if (is.null(df)) return()
  if (!is.data.frame(df)) df <- as.data.frame(df)

  # Clean existing group to prevent conflicts
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

    # Anti-Crash: Flatten list columns to semi-colon separated strings
    if (is.list(val)) {
      val <- vapply(val, function(x) paste(as.character(x), collapse = ";"), character(1))
    }

    if (is.factor(val)) val <- as.character(val)
    if (is.character(val)) val[is.na(val)] <- ''

    # Check for Atomic types before writing
    if (is.atomic(val)) {
      # Convert logicals to integers
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

    key_name <- paste0('X_', red)
    if (!is.matrix(emb)) emb <- as.matrix(emb)

    # Embeddings are (Cells x PCs). Direct write is correct.
    tryCatch({
      obsm_grp$create_dataset(key_name, emb)
    }, error = function(e) NULL)
  }
}

write_layers_h5 <- function(h5_file, seurat_obj) {
  assay_name <- Seurat::DefaultAssay(seurat_obj)
  
  # Seurat V5 Layer detection
  layer_names <- tryCatch({
    SeuratObject::Layers(seurat_obj, assay = assay_name)
  }, error = function(e) {
    names(seurat_obj[[assay_name]]@layers)
  })

  if (length(layer_names) == 0) return()

  if (!h5_file$exists('layers')) layers_grp <- h5_file$create_group('layers') else layers_grp <- h5_file[['layers']]

  for (lname in layer_names) {
    # Skip standard slots as they are handled in X or raw
    if (lname %in% c('data', 'counts')) next

    # Use 'layer' argument for Seurat V5 compatibility
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
