#' Convert H5AD to Seurat V5
#' @param file_path Path to .h5ad file
#' @return Seurat object
#' @export
h5ad_to_seurat <- function(file_path) {
  message('>>> [Import] Loading: ', file_path)

  if (!file.exists(file_path)) stop("File not found: ", file_path)

  file_h5 <- tryCatch(hdf5r::H5File$new(file_path, mode = 'r'), error = function(e) stop("Cannot open HDF5 file: ", e$message))
  on.exit(file_h5$close_all(), add = TRUE)

  obs <- read_h5ad_dataframe(file_h5, 'obs')
  var <- read_h5ad_dataframe(file_h5, 'var')

  counts_node <- NULL
  data_node <- NULL
  
  has_raw <- file_h5$exists('raw') && file_h5[['raw']]$exists('X')
  
  if (has_raw) {
    message("    Found 'raw/X', using as counts.")
    counts_node <- file_h5[['raw/X']]
    raw_var <- read_h5ad_dataframe(file_h5, 'raw/var')
    if (!is.null(raw_var)) var <- raw_var
    if (file_h5$exists('X')) data_node <- file_h5[['X']]
  } else {
    if (file_h5$exists('X')) counts_node <- file_h5[['X']]
  }

  if (is.null(counts_node)) stop('No expression matrix found (X or raw/X)')

  # Read as stored (H5AD usually Cells x Genes)
  counts_mat <- read_matrix_node(counts_node, transpose_output = FALSE)
  if (is.null(counts_mat)) stop("Failed to read counts matrix.")

  # Transpose to Seurat format (Genes x Cells)
  counts_mat <- Matrix::t(counts_mat)

  # Validate Dimensions
  dim_valid <- TRUE
  if (!is.null(obs) && ncol(counts_mat) != nrow(obs)) dim_valid <- FALSE
  if (!is.null(var) && nrow(counts_mat) != nrow(var)) dim_valid <- FALSE

  if (!dim_valid) {
    # Try reverting transpose (if H5AD was non-standard Genes x Cells)
    reverted_mat <- Matrix::t(counts_mat)
    reverted_valid <- TRUE
    if (!is.null(obs) && ncol(reverted_mat) != nrow(obs)) reverted_valid <- FALSE
    if (!is.null(var) && nrow(reverted_mat) != nrow(var)) reverted_valid <- FALSE
    
    if (reverted_valid) {
      message("    Detected non-standard matrix orientation. Adjusting.")
      counts_mat <- reverted_mat
    } else {
      warning("Dimension mismatch detected. Proceeding with standard orientation.")
    }
  }

  # Ensure names
  if (!is.null(obs)) colnames(counts_mat) <- rownames(obs)
  if (!is.null(var)) rownames(counts_mat) <- rownames(var)
  
  if (is.null(colnames(counts_mat))) colnames(counts_mat) <- paste0("Cell_", seq_len(ncol(counts_mat)))
  if (is.null(rownames(counts_mat))) rownames(counts_mat) <- paste0("Gene_", seq_len(nrow(counts_mat)))

  # [FIX] Ensure CsparseMatrix for Seurat Compatibility
  if (!inherits(counts_mat, "CsparseMatrix")) {
    counts_mat <- methods::as(counts_mat, "CsparseMatrix")
  }

  seu <- Seurat::CreateSeuratObject(counts = counts_mat, meta.data = obs, project = 'H5AD')

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

  if (!is.null(var)) {
    da <- Seurat::DefaultAssay(seu)
    for (i in colnames(var)) {
      tryCatch({ seu[[da]][[i]] <- var[[i]] }, error = function(e) NULL)
    }
    if ('highly_variable' %in% colnames(var)) {
      hv <- rownames(var)[which(var$highly_variable == TRUE)]
      if (length(hv) > 0) Seurat::VariableFeatures(seu) <- hv
    }
  }

  seu <- add_layers(seu, file_h5)
  seu <- add_reductions(seu, file_h5)

  message('>>> Conversion Complete!')
  return(seu)
}
