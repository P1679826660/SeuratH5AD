read_h5ad_dataframe <- function(h5_file, group_name = 'obs') {
  if (!h5_file$exists(group_name)) return(NULL)

  grp <- h5_file[[group_name]]
  attrs <- h5_node_attributes(grp)
  index_col_name <- if ('_index' %in% names(attrs)) attrs[['_index']] else '_index'

  df_list <- list()
  
  dataset_names <- grp$ls(recursive = FALSE)
  dataset_names <- dataset_names[dataset_names$link.type == 'H5L_TYPE_HARD', 'name']
  dataset_names <- dataset_names[!grepl('^__', dataset_names)]

  row_names_vals <- NULL

  for (name in dataset_names) {
    obj <- grp[[name]]
    
    if (inherits(obj, 'H5D')) {
      val <- tryCatch(obj$read(), error = function(e) NULL)
      if (is.null(val)) next


      obj_attrs <- h5_node_attributes(obj)
      if ('categories' %in% names(obj_attrs)) {
        cat_ref <- obj_attrs[['categories']]
        
        categories <- NULL
        tryCatch({
          if (inherits(cat_ref, "H5Ref")) {
             cat_dset <- h5_file[[cat_ref]]
             if (inherits(cat_dset, 'H5D')) {
               categories <- cat_dset$read()
             }
          }
        }, error = function(e) {
          message("Warning: Failed to resolve categories for ", name)
        })

        if (!is.null(categories)) {
          val_indices <- as.integer(val)
          
          # R is 1-based, HDF5 (Python) is 0-based
          val_r_indices <- val_indices + 1
          
          val_mapped <- categories[val_r_indices]
          
          val <- val_mapped
          
          # Optional: Convert to factor if needed
          # val <- factor(val, levels = categories)
        }
      }

      if (name == index_col_name) {
        row_names_vals <- val
      } else if (is.null(dim(val)) || length(dim(val)) == 1) {
        df_list[[name]] <- val
      }
    }
  }

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

#' @importFrom Seurat CreateDimReducObject SetAssayData
#' @importFrom Matrix t
NULL

add_reductions <- function(seu, file_h5) {
  if (!file_h5$exists('obsm')) return(seu)
  
  obsm_grp <- file_h5[['obsm']]
  red_names <- names(obsm_grp)
  
  for (name in red_names) {
    # Remove 'X_' prefix commonly found in H5AD (e.g., X_pca -> pca)
    clean_name <- sub("^X_", "", name)
    key_name <- paste0(toupper(clean_name), "_")
    
    mat <- tryCatch({
      obsm_grp[[name]]$read()
    }, error = function(e) NULL)
    
    if (is.null(mat)) next
    
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    
    # Check dimensions (rows in H5AD obsm must match columns in Seurat object)
    if (nrow(mat) != ncol(seu)) {
      message("Warning: Skipping reduction ", name, " due to dimension mismatch.")
      next
    }
    
    rownames(mat) <- colnames(seu)
    colnames(mat) <- paste0(key_name, seq_len(ncol(mat)))
    
    tryCatch({
      seu[[clean_name]] <- Seurat::CreateDimReducObject(
        embeddings = mat,
        key = key_name,
        assay = Seurat::DefaultAssay(seu)
      )
    }, error = function(e) {
      message("Warning: Could not create reduction ", clean_name)
    })
  }
  
  return(seu)
}

add_layers <- function(seu, file_h5) {
  if (!file_h5$exists('layers')) return(seu)
  
  layers_grp <- file_h5[['layers']]
  layer_names <- names(layers_grp)
  assay_name <- Seurat::DefaultAssay(seu)
  
  for (lname in layer_names) {
    node <- layers_grp[[lname]]
    
    # Read matrix using core utility
    mat <- read_matrix_node(node)
    
    if (is.null(mat)) next
    
    # Transpose to match Seurat format (Genes x Cells)
    mat <- Matrix::t(mat)
    
    # Validate dimensions
    if (ncol(mat) == ncol(seu) && nrow(mat) == nrow(seu)) {
      colnames(mat) <- colnames(seu)
      rownames(mat) <- rownames(seu)
      
      # Store in Seurat V5 layer
      tryCatch({
        seu <- Seurat::SetAssayData(seu, layer = lname, new.data = mat, assay = assay_name)
      }, error = function(e) {
        message("Warning: Could not set layer ", lname)
      })
    } else {
       message("Warning: Skipping layer ", lname, " due to dimension mismatch.")
    }
  }
  
  return(seu)
}
