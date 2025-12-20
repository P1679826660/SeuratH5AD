#' Convert H5AD to Seurat V5
#' @param file_path Path to .h5ad file
#' @return Seurat object
#' @export
h5ad_to_seurat <- function(file_path) {
  message('>>> [Import] Loading: ', file_path)

  if (!file.exists(file_path)) stop("File not found: ", file_path)

  file_h5 <- tryCatch(hdf5r::H5File$new(file_path, mode = 'r'), error = function(e) stop("Cannot open HDF5 file: ", e$message))
  on.exit(file_h5$close_all(), add = TRUE)

  # 1. Load Metadata
  obs <- read_h5ad_dataframe(file_h5, 'obs')
  var <- read_h5ad_dataframe(file_h5, 'var')

  # 2. Locate Expression Matrix
  counts_node <- NULL
  data_node <- NULL
  
  # Check for raw counts (raw/X)
  has_raw <- file_h5$exists('raw') && file_h5[['raw']]$exists('X')
  
  if (has_raw) {
    message("    Found 'raw/X', using as counts.")
    counts_node <- file_h5[['raw/X']]
    
    # Use raw/var if available
    raw_var <- read_h5ad_dataframe(file_h5, 'raw/var')
    if (!is.null(raw_var)) {
      var <- raw_var
    }
    
    # Use root X as data (normalized) if raw is used for counts
    if (file_h5$exists('X')) data_node <- file_h5[['X']]
    
  } else {
    # Default to root X
    if (file_h5$exists('X')) counts_node <- file_h5[['X']]
  }

  if (is.null(counts_node)) stop('No expression matrix found (X or raw/X)')

  # 3. Read Matrix
  # H5AD Standard: X is (Cells x Genes)
  # read_matrix_node returns the matrix as stored in H5
  counts_mat <- read_matrix_node(counts_node, transpose_output = FALSE)
  if (is.null(counts_mat)) stop("Failed to read counts matrix.")

  # 4. Align Dimensions
  # Seurat requires (Genes x Cells). H5AD is (Cells x Genes).
  # Default action: Transpose.
  counts_mat <- Matrix::t(counts_mat)

  # Verify dimensions against metadata
  # counts_mat is now (Genes x Cells)
  
  dim_valid <- TRUE
  if (!is.null(obs)) {
    if (ncol(counts_mat) != nrow(obs)) dim_valid <- FALSE
  }
  if (!is.null(var)) {
    if (nrow(counts_mat) != nrow(var)) dim_valid <- FALSE
  }

  if (!dim_valid) {
    # If dimensions do not match, check if the matrix was stored as (Genes x Cells) in H5AD (Non-standard)
    # Revert transpose and check
    reverted_mat <- Matrix::t(counts_mat)
    
    reverted_valid <- TRUE
    if (!is.null(obs) && ncol(reverted_mat) != nrow(obs)) reverted_valid <- FALSE
    if (!is.null(var) && nrow(reverted_mat) != nrow(var)) reverted_valid <- FALSE
    
    if (reverted_valid) {
      message("    Detected non-standard matrix orientation. Adjusting.")
      counts_mat <- reverted_mat
    } else {
      # If still invalid, warn but proceed with standard expectation
      warning("Dimension mismatch detected between matrix and metadata. Proceeding with standard orientation.")
    }
  }

  # 5. Assign Names
  if (!is.null(obs)) colnames(counts_mat) <- rownames(obs)
  if (!is.null(var)) rownames(counts_mat) <- rownames(var)

  if (is.null(colnames(counts_mat))) colnames(counts_mat) <- paste0("Cell_", seq_len(ncol(counts_mat)))
  if (is.null(rownames(counts_mat))) rownames(counts_mat) <- paste0("Gene_", seq_len(nrow(counts_mat)))

  # 6. Create Object
  seu <- Seurat::CreateSeuratObject(counts = counts_mat, meta.data = obs, project = 'H5AD')

  # 7. Add Data Slot
  if (!is.null(data_node)) {
    data_mat <- read_matrix_node(data_node, transpose_output = FALSE)
    if (!is.null(data_mat)) {
      data_mat <- Matrix::t(data_mat)
      
      if (ncol(data_mat) == ncol(seu) && nrow(data_mat) == nrow(seu)) {
        colnames(data_mat) <- colnames(seu)
        rownames(data_mat) <- rownames(seu)
        seu <- Seurat::SetAssayData(seu, slot = 'data', new.data = data_mat)
      }
    }
  }

  # 8. Add Feature Metadata
  if (!is.null(var)) {
    da <- Seurat::DefaultAssay(seu)
    for (i in colnames(var)) {
      tryCatch({
        seu[[da]][[i]] <- var[[i]]
      }, error = function(e) NULL)
    }
    
    if ('highly_variable' %in% colnames(var)) {
      hv <- rownames(var)[which(var$highly_variable == TRUE)]
      if (length(hv) > 0) Seurat::VariableFeatures(seu) <- hv
    }
  }

  # 9. Add Reductions and Layers
  seu <- add_layers(seu, file_h5)
  seu <- add_reductions(seu, file_h5)

  message('>>> Conversion Complete!')
  return(seu)
}
