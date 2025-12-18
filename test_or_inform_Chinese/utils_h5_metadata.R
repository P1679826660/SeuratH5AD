# ==============================================================================
# Module Name: H5AD Metadata Reader (V2 Fixed)
# Function Description: Reads obs or var groups and converts them to DataFrames, enhanced with handling for empty column scenarios.
# ==============================================================================

#' Read obs or var groups and convert to DataFrame
read_h5ad_dataframe <- function(h5_file, group_name = "obs") {
  if (!h5_file$exists(group_name)) return(NULL)
  
  grp <- h5_file[[group_name]]
  
  # 1. Find Index
  # AnnData usually stores the index column name in the group attribute '_index'
  attrs <- h5_node_attributes(grp)
  index_col_name <- if("_index" %in% names(attrs)) attrs[["_index"]] else "_index"
  
  # 2. Iterate and read columns
  df_list <- list()
  dataset_names <- grp$ls(recursive = FALSE)$name
  
  # Filter out special __categories (if any)
  dataset_names <- dataset_names[!grepl("^__", dataset_names)]
  
  row_names_vals <- NULL
  
  for (name in dataset_names) {
    obj <- grp[[name]]
    
    # Process Datasets only
    if (inherits(obj, "H5D")) {
      val <- obj$read()
      
      if (name == index_col_name) {
        row_names_vals <- val
      } else {
        # Ensure reading a vector; if it is a multi-dimensional array (other than 1D), it might need handling
        if (is.null(dim(val)) || length(dim(val)) == 1) {
          df_list[[name]] <- val
        } else {
          # If a multi-dimensional array is encountered as a metadata column (rare), skip it or take the first column
          # For engineering stability, choosing to warn and skip here to prevent breaking the data.frame structure
          warning(paste("Skipping metadata column with complex structure:", name))
        }
      }
    }
  }
  
  # 3. Construct DataFrame (Core Fix Section)
  
  # Case A: Only index, no other data columns (df_list is empty)
  if (length(df_list) == 0) {
    if (!is.null(row_names_vals)) {
      # Create a DataFrame with only row names and no columns
      # Must explicitly specify row.names to establish row count
      df <- data.frame(row.names = row_names_vals)
    } else {
      # Neither data nor index exists, return NULL or empty object
      return(NULL)
    }
  } else {
    # Case B: Has data columns
    # Construct data frame first, do not set row names yet to avoid immediate error due to length mismatch
    df <- data.frame(df_list, stringsAsFactors = FALSE)
    
    # Safely set row names
    if (!is.null(row_names_vals)) {
      # 1. Check for duplicates
      if (anyDuplicated(row_names_vals)) {
        warning(paste0(group_name, " index contains duplicates, deduplicating via make.unique..."))
        row_names_vals <- make.unique(as.character(row_names_vals))
      }
      
      # 2. Check length consistency
      if (length(row_names_vals) == nrow(df)) {
        rownames(df) <- row_names_vals
      } else {
        # Critical Warning: Index length does not match data length
        warning(sprintf(
          "Warning: In %s, index length (%d) does not match data column length (%d). Ignoring original index, using numeric index.",
          group_name, length(row_names_vals), nrow(df)
        ))
        # Keep numeric index at this point, do not force assignment to prevent crash
      }
    }
  }
  
  return(df)
}

# Ensure helper functions exist (If your convert_h5ad.R does not source this independently, keep it here or ensure the matrix module is loaded)
# For safety, duplicate definitions are avoided here assuming it is defined in utils_h5_matrix.R.
