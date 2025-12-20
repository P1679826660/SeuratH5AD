read_h5ad_dataframe <- function(h5_file, group_name = 'obs') {
  if (!h5_file$exists(group_name)) return(NULL)

  grp <- h5_file[[group_name]]
  attrs <- h5_node_attributes(grp)
  index_col_name <- if ('_index' %in% names(attrs)) attrs[['_index']] else '_index'

  df_list <- list()
  
  # Read datasets only (HARD links), skipping Groups (complex categoricals) for stability
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

  # Robust DataFrame Construction
  df <- tryCatch({
    if (length(df_list) > 0) {
      max_len <- if (!is.null(row_names_vals)) length(row_names_vals) else max(sapply(df_list, length))
      
      # Pad or truncate columns to ensure equal length
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
    message("    Warning: Could not construct dataframe for ", group_name)
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)

  if (!is.null(row_names_vals)) {
    row_names_vals <- as.character(row_names_vals)
    if (anyDuplicated(row_names_vals)) row_names_vals <- make.unique(row_names_vals)
    if (length(row_names_vals) == nrow(df)) rownames(df) <- row_names_vals
  }
  return(df)
}

add_reductions_standard <- function(seu, h5_file) {
  # Standard: obsm contains (Cells x PCs) matrices.
  if (!h5_file$exists('obsm')) return(seu)
  
  grp <- h5_file[['obsm']]
  cell_names <- colnames(seu)
  n_cells <- length(cell_names)

  items <- grp$ls(recursive = FALSE)
  for (k in items$name) {
    obj <- grp[[k]]
    # Skip non-Datasets
    if (!inherits(obj, 'H5D')) next

    # Read and convert to matrix
    raw_emb <- tryCatch(obj$read(), error = function(e) NULL)
    if (is.null(raw_emb)) next
    if (!is.matrix(raw_emb)) raw_emb <- as.matrix(raw_emb)

    # 1. Strict validation: Rows must equal Cell count
    if (nrow(raw_emb) != n_cells) {
       # 2. Transpose check (handling rare non-standard files)
       if (ncol(raw_emb) == n_cells) {
         raw_emb <- t(raw_emb)
       } else {
         # Dimension mismatch, skip silently to prevent errors
         next 
       }
    }

    # Name cleaning: X_pca -> pca
    clean_name <- tolower(gsub('^X_', '', k)) 
    if (nchar(clean_name) == 0) clean_name <- "dimred"
    
    seurat_key <- paste0(clean_name, '_')
    
    rownames(raw_emb) <- cell_names
    colnames(raw_emb) <- paste0(seurat_key, seq_len(ncol(raw_emb)))

    tryCatch({
      dr_obj <- Seurat::CreateDimReducObject(embeddings = raw_emb, key = seurat_key, assay = Seurat::DefaultAssay(seu))
      seu[[clean_name]] <- dr_obj
      message("    Added reduction: ", clean_name)
    }, error = function(e) {
      # Ignore failures
    })
  }
  return(seu)
}
