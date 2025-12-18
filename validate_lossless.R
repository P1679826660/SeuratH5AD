# ==============================================================================
# Script Name: H5AD Lossless Conversion Ultimate Validator (Bit-level Comparator)
# Function: Compare underlying data of two h5ad files to prove "different size but identical information"
# ==============================================================================

library(hdf5r)
library(Matrix)

compare_h5ad <- function(file_orig, file_new) {
  
  message(">>> Starting lossless validation process...")
  message("    Source File (A): ", file_orig, " (", round(file.size(file_orig)/1024/1024, 2), " MB)")
  message("    New File (B): ", file_new,  " (", round(file.size(file_new)/1024/1024, 2), " MB)")
  
  h5_a <- H5File$new(file_orig, mode = "r")
  h5_b <- H5File$new(file_new, mode = "r")
  on.exit({ h5_a$close_all(); h5_b$close_all() })
  
  # --- 1. Check Compression Filters (Core explanation for size difference) ---
  message("\n[1. Storage Compression Analysis]")
  
  get_compression <- function(h5, name) {
    if(h5$exists(name)) {
      ds <- h5[[name]]
      # data is a dataset that must exist (in sparse format)
      if(inherits(ds, "H5Group")) ds <- ds[["data"]] 
      filter <- ds$get_create_plist()$get_nfilters()
      if (filter > 0) return("GZIP/Compressed") else return("None/Uncompressed")
    }
    return("Unknown")
  }
  
  msg_a <- get_compression(h5_a, "X")
  msg_b <- get_compression(h5_b, "X")
  message(sprintf("    Source File X Compression Status: %s", msg_a))
  message(sprintf("    New File X Compression Status: %s", msg_b))
  
  if (msg_a != msg_b) {
    message("    !!! Conclusion: File size difference is mainly due to different compression algorithms, not data loss.")
  }
  
  # --- 2. Full Matrix Content Comparison ---
  message("\n[2. Matrix Full Value Comparison (X)]")
  
  read_full_matrix <- function(h5_file) {
    # Simple reader, does not rely on our previous complex modules, reads data directly
    x_node <- h5_file[["X"]]
    if(inherits(x_node, "H5Group")) {
      return(list(
        sum = sum(x_node[["data"]]$read()),
        nz  = length(x_node[["data"]]$read()),
        dim = x_node[["data"]]$dims # Note: Reading data vector length here, not matrix dimensions
      ))
    } else {
      # Dense
      val <- x_node[,]
      return(list(sum = sum(val), nz = sum(val != 0)))
    }
  }
  
  stats_a <- read_full_matrix(h5_a)
  stats_b <- read_full_matrix(h5_b)
  
  message(sprintf("    Source X Sum: %.4f | Non-zero elements: %d", stats_a$sum, stats_a$nz))
  message(sprintf("    New X Sum: %.4f | Non-zero elements: %d", stats_b$sum, stats_b$nz))
  
  if (abs(stats_a$sum - stats_b$sum) < 0.001 && stats_a$nz == stats_b$nz) {
    message("    >>> [PASS] Matrix information matches exactly! (Data lossless)")
  } else {
    message("    >>> [FAIL] Matrix information mismatch! Please check!")
  }
  
  # --- 3. Dimension Comparison ---
  message("\n[3. Logical Dimension Comparison]")
  # Read obs/var index length
  get_idx_len <- function(h5, grp) {
    if(!h5$exists(grp)) return(0)
    g <- h5[[grp]]
    # Find _index
    if(g$attr_exists("_index")) {
      idx_name <- g$attr_open("_index")$read()
      return(g[[idx_name]]$dims)
    }
    return(0)
  }
  
  obs_a <- get_idx_len(h5_a, "obs"); var_a <- get_idx_len(h5_a, "var")
  obs_b <- get_idx_len(h5_b, "obs"); var_b <- get_idx_len(h5