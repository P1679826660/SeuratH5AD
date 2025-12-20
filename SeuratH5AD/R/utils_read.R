read_h5ad_dataframe <- function(h5_file, group_name = 'obs') {
  if (!h5_file$exists(group_name)) return(NULL)

  grp <- h5_file[[group_name]]
  attrs <- h5_node_attributes(grp)
  index_col_name <- if ('_index' %in% names(attrs)) attrs[['_index']] else '_index'

  df_list <- list()
  
  # Retrieve dataset names, filtering out internal HDF5 references
  dataset_names <- grp$ls(recursive = FALSE)
  dataset_names <- dataset_names[dataset_names$link.type == 'H5L_TYPE_HARD', 'name']
  dataset_names <- dataset_names[!grepl('^__', dataset_names)]

  row_names_vals <- NULL

  for (name in dataset_names) {
    obj <- grp[[name]]
    if (inherits(obj, 'H5D')) {
      val <- tryCatch(obj$read(), error = function(e) NULL)
      if (is.null(val)) next

      if (name == index_col_name) {
        row_names_vals <- val
      } else if (is.null(dim(val)) || length(dim(val)) == 1) {
        df_list[[name]] <- val
      }
    }
  }

  if (length(df_list) == 0 && is.null(row_names_vals)) return(NULL)

  # Construct DataFrame safely
  if (length(df_list) > 0) {
    max_len <- if (!is.null(row_names_vals)) length(row_names_vals) else max(sapply(df_list, length))

    safe_list <- lapply(df_list, function(v) {
      if (length(v) == max_len) return(v)
      if (length(v) < max_len) return(c(v, rep(NA, max_len - length(v))))
      return(v[1:max_len])
    })

    df <- data.frame(safe_list, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(row.names = seq_along(row_names_vals))
  }

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

add_reductions <- function(seu, h5_file) {
  emb_group_name <- 'obsm'
  cell_names <- colnames(seu)

  if (h5_file$exists(emb_group_name)) {
    grp <- h5_file[[emb_group_name]]
    for (k in grp$ls(recursive = FALSE)$name) {
      obj <- grp[[k]]
      if (!inherits(obj, 'H5D')) next

      raw_emb <- tryCatch(obj$read(), error = function(e) NULL)
      if (is.null(raw_emb)) next
      if (!is.matrix(raw_emb)) raw_emb <- as.matrix(raw_emb)

      # Standard H5AD obsm is (Cells x PCs)
      # Check dimension compatibility
      if (nrow(raw_emb) == length(cell_names)) {
        # Correct orientation
      } else if (ncol(raw_emb) == length(cell_names)) {
        # Transposed orientation, correct it
        raw_emb <- t(raw_emb)
      } else {
        # Dimension mismatch, skip
        next
      }

      clean_name <- tolower(gsub('^X_', '', k))
      if (nchar(clean_name) == 0) clean_name <- "dimred"

      seurat_key <- paste0(clean_name, '_')
      rownames(raw_emb) <- cell_names
      colnames(raw_emb) <- paste0(seurat_key, seq_len(ncol(raw_emb)))

      dr_obj <- Seurat::CreateDimReducObject(embeddings = raw_emb, key = seurat_key, assay = Seurat::DefaultAssay(seu))
      seu[[clean_name]] <- dr_obj
    }
  }
  return(seu)
}

add_layers <- function(seu, h5_file) {
  if (!h5_file$exists('layers')) return(seu)
  grp <- h5_file[['layers']]

  for (lname in grp$ls(recursive = FALSE)$name) {
    # Layers follow X dimensions (Cells x Genes)
    mat <- read_matrix_node(grp[[lname]], transpose_output = FALSE)
    if (is.null(mat)) next

    # Transpose to (Genes x Cells) for Seurat
    mat <- Matrix::t(mat)

    if (ncol(mat) == ncol(seu) && nrow(mat) == nrow(seu)) {
      target_slot <- if (lname %in% c('logcounts', 'norm_data')) 'data' else if (grepl('scale', lname)) 'scale.data' else 'data'
      da <- Seurat::DefaultAssay(seu)

      tryCatch({
        if (target_slot == 'scale.data') {
          Seurat::SetAssayData(seu, slot = 'scale.data', new.data = as.matrix(mat), assay = da)
        } else {
          Seurat::SetAssayData(seu, slot = 'data', new.data = mat, assay = da)
        }
      }, error = function(e) {
        warning(sprintf("Could not add layer %s: %s", lname, e$message))
      })
    }
  }
  return(seu)
}
