#' Convert H5AD to Seurat V5 (Robust Standard Mode)
#' @param file_path Path to .h5ad file
#' @return Seurat object
#' @export
h5ad_to_seurat <- function(file_path) {
  message('>>> [Import] Loading: ', file_path)

  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  file_h5 <- tryCatch(hdf5r::H5File$new(file_path, mode = 'r'), error = function(e) stop("Cannot open HDF5 file: ", e$message))
  on.exit(file_h5$close_all(), add = TRUE)

  # 1. Metadata
  obs <- read_h5ad_dataframe(file_h5, 'obs')
  var <- read_h5ad_dataframe(file_h5, 'var')

  # ------------------------------------------------------------------
  # 2. Locate Counts (Raw/Original Data)
  # ------------------------------------------------------------------
  counts_mat <- NULL
  counts_source <- ""
  var_for_counts <- NULL 
  
  # Priority 1: raw/X (Standard Scanpy location for full counts)
  if (file_h5$exists('raw') && file_h5[['raw']]$exists('X')) {
    message("    Found 'raw/X', reading as Counts (Full Gene Set).")
    counts_mat <- read_matrix_node(file_h5[['raw/X']])
    counts_source <- "raw"
    
    # Use raw/var if available
    raw_var <- read_h5ad_dataframe(file_h5, 'raw/var')
    var_for_counts <- if (!is.null(raw_var)) raw_var else var

  } else if (file_h5$exists('layers') && file_h5[['layers']]$exists('counts')) {
    # Priority 2: layers/counts
    message("    Found 'layers/counts', reading as Counts.")
    counts_mat <- read_matrix_node(file_h5[['layers/counts']])
    counts_source <- "layer"
    var_for_counts <- var
    
  } else if (file_h5$exists('X')) {
    # Priority 3: Root X
    message("    Using 'X' as Counts (No raw or layers found).")
    counts_mat <- read_matrix_node(file_h5[['X']])
    counts_source <- "X"
    var_for_counts <- var
  }

  if (is.null(counts_mat)) stop("Critical Error: No expression matrix found (checked raw/X, layers/counts, X).")
  
  # Transpose: H5AD (Cells x Genes) -> Seurat (Genes x Cells)
  counts_mat <- Matrix::t(counts_mat)
  
  # ------------------------------------------------------------------
  # 3. Safe Dimension Naming
  # ------------------------------------------------------------------
  
  # Cells
  if (!is.null(obs) && ncol(counts_mat) == nrow(obs)) {
    colnames(counts_mat) <- rownames(obs)
  } else {
    if (!is.null(obs)) message("    Warning: Cell count mismatch (Matrix: ", ncol(counts_mat), ", Obs: ", nrow(obs), "). Using generic cell names.")
    colnames(counts_mat) <- paste0("Cell_", seq_len(ncol(counts_mat)))
  }

  # Genes
  if (!is.null(var_for_counts) && nrow(counts_mat) == nrow(var_for_counts)) {
    rownames(counts_mat) <- rownames(var_for_counts)
  } else {
    if (!is.null(var_for_counts)) message("    Warning: Gene count mismatch (Matrix: ", nrow(counts_mat), ", Var: ", nrow(var_for_counts), "). Using generic gene names.")
    rownames(counts_mat) <- paste0("Gene_", seq_len(nrow(counts_mat)))
  }

  # Ensure CsparseMatrix for Seurat V5
  if (!inherits(counts_mat, "CsparseMatrix")) {
    counts_mat <- methods::as(counts_mat, "CsparseMatrix")
  }

  # Create Seurat Object
  seu <- Seurat::CreateSeuratObject(counts = counts_mat, meta.data = obs, project = 'H5AD')
  message("    Seurat Object created: ", nrow(seu), " features x ", ncol(seu), " cells.")

  # ------------------------------------------------------------------
  # 4. Attempt to Load Normalized Data (X)
  # ------------------------------------------------------------------
  # Only if we used 'raw' for counts, 'X' might contain normalized data
  if (counts_source == "raw" && file_h5$exists('X')) {
    message("    Attempting to load Normalized Data from 'X'...")
    data_mat <- read_matrix_node(file_h5[['X']])
    
    if (!is.null(data_mat)) {
      data_mat <- Matrix::t(data_mat) 
      
      # [Dimension Check] Only add if dimensions match exactly.
      if (ncol(data_mat) == ncol(seu) && nrow(data_mat) == nrow(seu)) {
        colnames(data_mat) <- colnames(seu)
        rownames(data_mat) <- rownames(seu)
        seu <- Seurat::SetAssayData(seu, layer = 'data', new.data = data_mat)
        message("    Success: Normalized data added.")
      } else {
        # Graceful skip for gene subsets (e.g. X contains only HVGs)
        message("    Skipped: 'X' dimensions (", nrow(data_mat), " genes) do not match Counts (", nrow(seu), " genes).")
        message("    Reason: Scanpy often stores a subset of genes in X. Seurat requires full alignment.")
      }
    }
  }

  # ------------------------------------------------------------------
  # 5. Metadata & Variable Features
  # ------------------------------------------------------------------
  if (!is.null(var)) {
    # Only add metadata for genes present in the Seurat object
    common_genes <- intersect(rownames(var), rownames(seu))
    
    if (length(common_genes) > 0) {
      da <- Seurat::DefaultAssay(seu)
      
      # Identify Variable Features (Scanpy standard: 'highly_variable')
      hv_cols <- c('highly_variable', 'highly_variable_genes')
      found_hv <- intersect(hv_cols, colnames(var))
      
      if (length(found_hv) > 0) {
        hv_col <- found_hv[1]
        hv_genes <- rownames(var)[which(var[[hv_col]] == TRUE)]
        hv_genes <- intersect(hv_genes, rownames(seu)) 
        if (length(hv_genes) > 0) {
          Seurat::VariableFeatures(seu) <- hv_genes
          message("    Set ", length(hv_genes), " variable features.")
        }
      }
    }
  }

  # ------------------------------------------------------------------
  # 6. Reductions
  # ------------------------------------------------------------------
  seu <- add_reductions_standard(seu, file_h5)

  message('>>> Conversion Complete!')
  return(seu)
}
