#' @importFrom hdf5r h5attr
NULL

read_h5ad_dataframe <- function(h5_file, group_name = 'obs') {
  if (!h5_file$exists(group_name)) return(NULL)

  grp <- h5_file[[group_name]]
  attrs <- h5_node_attributes(grp)
  index_col_name <- if ('_index' %in% names(attrs)) attrs[['_index']] else '_index'

  df_list <- list()
  
  all_names <- grp$ls(recursive = FALSE)$name
  all_names <- all_names[!grepl('^__', all_names)]

  row_names_vals <- NULL

  for (name in all_names) {
    obj <- grp[[name]]
    

    if (inherits(obj, 'H5D')) {
      val <- tryCatch(obj$read(), error = function(e) NULL)
      if (is.null(val)) next

      if (name == index_col_name) {
        row_names_vals <- val
      } else if (is.null(dim(val)) || length(dim(val)) == 1) {
        df_list[[name]] <- val
      }
      

    # -----------------------------------------------------------
    } else if (inherits(obj, 'H5Group')) {
      if (obj$exists('categories') && obj$exists('codes')) {
        tryCatch({
          categories <- obj[['categories']]$read()
          codes <- obj[['codes']]$read()
          
          
          
          final_val <- rep(NA, length(codes))
          valid_mask <- codes >= 0
          
          if (any(valid_mask)) {
            r_indices <- codes[valid_mask] + 1
            
            safe_indices_mask <- r_indices <= length(categories)
            
            
            final_indices <- r_indices[safe_indices_mask]
            
            
            valid_codes_subset <- codes[valid_mask]
            mapped_values <- categories[valid_codes_subset + 1]
            
            final_val[valid_mask] <- mapped_values
          }
          
          df_list[[name]] <- final_val
          
        }, error = function(e) {
          
          message("    Warning: Failed to decode categorical '", name, "'. Skipping.")
        })
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
  if (!h5_file$exists('obsm')) return(seu)
  
  grp <- h5_file[['obsm']]
  cell_names <- colnames(seu)
  n_cells <- length(cell_names)

  items <- grp$ls(recursive = FALSE)
  for (k in items$name) {
    obj <- grp[[k]]
    if (!inherits(obj, 'H5D')) next

    raw_emb <- tryCatch(obj$read(), error = function(e) NULL)
    if (is.null(raw_emb)) next
    if (!is.matrix(raw_emb)) raw_emb <- as.matrix(raw_emb)

    if (nrow(raw_emb) != n_cells) {
       if (ncol(raw_emb) == n_cells) {
         raw_emb <- t(raw_emb)
       } else {
         next 
       }
    }

    clean_name <- tolower(gsub('^X_', '', k)) 
    if (nchar(clean_name) == 0) clean_name <- "dimred"
    
    seurat_key <- paste0(clean_name, '_')
    
    rownames(raw_emb) <- cell_names
    colnames(raw_emb) <- paste0(seurat_key, seq_len(ncol(raw_emb)))

    tryCatch({
      dr_obj <- Seurat::CreateDimReducObject(embeddings = raw_emb, key = seurat_key, assay = Seurat::DefaultAssay(seu))
      seu[[clean_name]] <- dr_obj
      message("    Added reduction: ", clean_name)
    }, error = function(e) NULL)
  }
  return(seu)
}
