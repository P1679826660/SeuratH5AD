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
          
          val_r_indices <- val_indices + 1
          
          val_mapped <- categories[val_r_indices]
          
          val <- val_mapped
          
          # val <- factor(val, levels = categories)
        }
      }
      # --- [FIX END] ---

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

