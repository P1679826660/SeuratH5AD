#' @importFrom Seurat CreateDimReducObject SetAssayData
#' @importFrom Matrix t
#' @importFrom hdf5r H5Ref
NULL

read_h5ad_dataframe <- function(h5_file, group_name = 'obs') {
  if (!h5_file$exists(group_name)) return(NULL)

  grp <- h5_file[[group_name]]
  attrs <- h5_node_attributes(grp)
  index_col_name <- if ('_index' %in% names(attrs)) attrs[['_index']] else '_index'

  df_list <- list()
  
  # [CRITICAL FIX 1] Do not filter by Link Type, keep all columns
  # Previous code filtered out Soft Links, causing some columns to be missing
  dataset_names <- names(grp)
  dataset_names <- dataset_names[!grepl('^__', dataset_names)]

  row_names_vals <- NULL

  for (name in dataset_names) {
    # Wrap the entire reading process in tryCatch to prevent single column errors from crashing the function
    tryCatch({
      obj <- grp[[name]]
      
      # Ensure it is a Dataset (H5D)
      if (inherits(obj, 'H5D')) {
        val <- obj$read()
        if (is.null(val)) next

        # ----------------------------------------------------------------
        # Ultimate Categories Resolution Logic (Schard + Legacy Support)
        # ----------------------------------------------------------------
        categories <- NULL
        is_categorical <- FALSE
        
        # 1. Priority: HDF5 Object Reference (Standard AnnData)
        if (obj$attr_exists('categories')) {
          tryCatch({
            cat_ref <- obj$h5attr('categories')
            
            # If it is a reference (H5Ref), dereference and read
            if (inherits(cat_ref, "H5Ref")) {
               categories <- h5_file[[cat_ref]]$read()
               is_categorical <- TRUE
            } 
            # If it is a path string (compatibility for some older files)
            else if (is.character(cat_ref) && length(cat_ref) == 1) {
               if (h5_file$exists(cat_ref)) {
                 categories <- h5_file[[cat_ref]]$read()
                 is_categorical <- TRUE
               }
            }
          }, error = function(e) NULL)
        }
        
        # 2. Fallback A: Scanpy internal storage location (obs/__categories/col_name)
        if (!is_categorical) {
          fallback_path <- paste0(group_name, "/__categories/", name)
          if (h5_file$exists(fallback_path)) {
             categories <- h5_file[[fallback_path]]$read()
             is_categorical <- TRUE
          }
        }
        
        # 3. [CRITICAL FIX 2] Fallback B: Legacy 'uns' storage (uns/col_name_categories)
        # Many legacy or converted H5AD files store dictionaries in 'uns'
        if (!is_categorical && h5_file$exists('uns')) {
           legacy_path <- paste0("uns/", name, "_categories")
           if (h5_file$exists(legacy_path)) {
             categories <- h5_file[[legacy_path]]$read()
             is_categorical <- TRUE
           }
        }

        # ----------------------------------------------------------------
        # Apply Mapping: Integer -> Text
        # ----------------------------------------------------------------
        if (is_categorical && !is.null(categories)) {
          # Force conversion to integer (handle potential float storage)
          val_indices <- as.integer(val)
          
          # Python (0-based) -> R (1-based)
          val_r_indices <- val_indices + 1
          
          # Handle -1 (Python NA) -> Map to R NA
          # Any index <= 0 is invalid in R, set to NA
          val_r_indices[val_r_indices <= 0] <- NA
          
          # Safe mapping
          if (length(categories) > 0) {
             val_mapped <- categories[val_r_indices]
             # Convert to factor to preserve order
             val <- factor(val_mapped, levels = categories)
          }
        }
        
        # ----------------------------------------------------------------
        # Store results
        # ----------------------------------------------------------------
        if (name == index_col_name) {
          row_names_vals <- as.character(val)
        } else {
          # Dimensionality reduction: If 1D array (N x 1), convert to vector
          if (!is.null(dim(val)) && length(dim(val)) == 1) val <- as.vector(val)
          
          # Only add pure vectors to DataFrame
          if (is.null(dim(val))) {
             df_list[[name]] <- val
          }
        }
      }
    }, error = function(e) {
      # Failure to read a single column should not crash the function, but warn
      warning(paste("Skipping column due to error:", name, e$message))
    })
  }

  if (length(df_list) == 0 && is.null(row_names_vals)) return(NULL)

  # Construct DataFrame
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

  # Handle Row Names
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

add_reductions <- function(seu, file_h5) {
  if (!file_h5$exists('obsm')) return(seu)
  
  obsm_grp <- file_h5[['obsm']]
  # Use names() instead of ls() to avoid link type filtering issues
  red_names <- names(obsm_grp)
  
  for (name in red_names) {
    clean_name <- sub("^X_", "", name)
    key_name <- paste0(toupper(clean_name), "_")
    
    mat <- tryCatch({
      obsm_grp[[name]]$read()
    }, error = function(e) NULL)
    
    if (is.null(mat)) next
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    
    if (nrow(mat) != ncol(seu)) next
    
    rownames(mat) <- colnames(seu)
    colnames(mat) <- paste0(key_name, seq_len(ncol(mat)))
    
    tryCatch({
      seu[[clean_name]] <- Seurat::CreateDimReducObject(
        embeddings = mat,
        key = key_name,
        assay = Seurat::DefaultAssay(seu)
      )
    }, error = function(e) NULL)
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
      }, error = function(e) NULL)
    }
  }
  return(seu)
}
