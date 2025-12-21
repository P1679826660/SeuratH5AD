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


      categories <- NULL
      is_categorical <- FALSE
      
      if (obj$attr_exists('categories')) {
        tryCatch({
          cat_attr <- obj$h5attr('categories')
          
          if (inherits(cat_attr, "H5Ref")) {
             cat_dset <- h5_file[[cat_attr]] # 解引用
             categories <- cat_dset$read()
             is_categorical <- TRUE
          } 
          else if (is.character(cat_attr) && length(cat_attr) == 1) {
             if (h5_file$exists(cat_attr)) {
               categories <- h5_file[[cat_attr]]$read()
               is_categorical <- TRUE
             }
          }
        }, error = function(e) NULL)
      }
      

      if (!is_categorical) {
        fallback_path <- paste0(group_name, "/__categories/", name)
        if (h5_file$exists(fallback_path)) {
           categories <- h5_file[[fallback_path]]$read()
           is_categorical <- TRUE
        }
      }

      if (is_categorical && !is.null(categories)) {
        val_indices <- as.integer(val)
        
        # Python 是 0-based，R 是 1-based
        val_r_indices <- val_indices + 1
        
        val_r_indices[val_r_indices <= 0] <- NA
        
        val_mapped <- categories[val_r_indices]
        
       
        val <- factor(val_mapped, levels = categories)
      }

      if (name == index_col_name) {
        row_names_vals <- as.character(val)
      } else {
        if (!is.null(dim(val)) && length(dim(val)) == 1) val <- as.vector(val)
        
        if (is.null(dim(val))) {
           df_list[[name]] <- val
        }
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
        return(v[1:max_len]) # 截断
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
#' @importFrom hdf5r H5Ref
NULL

add_reductions <- function(seu, file_h5) {
  if (!file_h5$exists('obsm')) return(seu)
  
  obsm_grp <- file_h5[['obsm']]
  red_names <- names(obsm_grp)
  
  for (name in red_names) {
    clean_name <- sub("^X_", "", name)
    key_name <- paste0(toupper(clean_name), "_")
    
    mat <- tryCatch({
      obsm_grp[[name]]$read()
    }, error = function(e) NULL)
    
    if (is.null(mat)) next
    
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    
    if (nrow(mat) != ncol(seu)) {
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
      # ignore
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
    mat <- read_matrix_node(node)
    
    if (is.null(mat)) next
    
    mat <- Matrix::t(mat)
    
    if (ncol(mat) == ncol(seu) && nrow(mat) == nrow(seu)) {
      colnames(mat) <- colnames(seu)
      rownames(mat) <- rownames(seu)
      
      tryCatch({
        seu <- Seurat::SetAssayData(seu, layer = lname, new.data = mat, assay = assay_name)
      }, error = function(e) {
         # ignore
      })
    }
  }
  
  return(seu)
}
