# ==============================================================================
# Module Name: H5AD Components Loader
# Function Description: Load Layers, Embeddings (obsm), Loadings (varm), Graphs (obsp)
# ==============================================================================

#' Read and inject dimensionality reduction information (Reductions)
#' @param seu Seurat Object
#' @param h5_file H5File handle
#' @param transpose_mode The "transpose mode" flag from convert_h5ad
add_reductions <- function(seu, h5_file, transpose_mode) {
  
  # === 1. Determine who stores Cell Embeddings (Coordinates) ===
  # Standard mode: obsm stores coordinates (Cells x Dims)
  # Transpose mode: varm stores coordinates (because var became Cells)
  emb_group_name <- if(transpose_mode) "varm" else "obsm"
  loadings_group_name <- if(transpose_mode) "obsm" else "varm"
  
  cell_names <- colnames(seu)
  
  if (h5_file$exists(emb_group_name)) {
    grp <- h5_file[[emb_group_name]]
    keys <- grp$ls(recursive=FALSE)$name
    
    for (k in keys) {
      # Filter out non-matrix data
      if (!inherits(grp[[k]], "H5D") && !inherits(grp[[k]], "H5Group")) next 
      
      # Read coordinates
      # Note: obsm reading is usually (N_cells, N_dims).
      # If it is a Dense Dataset, reading into R might result in (N_dims, N_cells) because R reads HDF5 column-major.
      raw_emb <- tryCatch(grp[[k]]$read(), error=function(e) NULL)
      if (is.null(raw_emb)) next
      
      # Ensure it is a matrix
      if (!is.matrix(raw_emb)) raw_emb <- as.matrix(raw_emb)
      
      # Dimension alignment check
      # We need rows = number of cells
      if (ncol(raw_emb) == length(cell_names)) {
        raw_emb <- t(raw_emb)
      }
      
      if (nrow(raw_emb) != length(cell_names)) {
        # Dimensions still don't match, skip
        next
      }
      
      # Name normalization X_pca -> pca
      clean_name <- tolower(gsub("^X_", "", k))
      seurat_key <- paste0(clean_name, "_")
      
      rownames(raw_emb) <- cell_names
      colnames(raw_emb) <- paste0(seurat_key, 1:ncol(raw_emb))
      
      # Create DR object
      dr_obj <- Seurat::CreateDimReducObject(
        embeddings = raw_emb,
        key = seurat_key,
        assay = Seurat::DefaultAssay(seu)
      )
      
      # === 2. Attempt to read corresponding Feature Loadings ===
      # If there is a corresponding varm (gene weights)
      if (h5_file$exists(loadings_group_name)) {
        # Guess name: PCA corresponds to PCs, usually keys match or are similar
        # Scanpy convention: obsm['X_pca'] <-> varm['PCs']
        loadings_key <- NULL
        if (clean_name == "pca" && h5_file[[loadings_group_name]]$exists("PCs")) loadings_key <- "PCs"
        else if (h5_file[[loadings_group_name]]$exists(k)) loadings_key <- k
        
        if (!is.null(loadings_key)) {
          raw_load <- h5_file[[loadings_group_name]][[loadings_key]]$read()
          if (!is.matrix(raw_load)) raw_load <- as.matrix(raw_load)
          
          # Loadings need to be (Features x Dims)
          # Check dimensions
          n_features <- nrow(seu)
          if (ncol(raw_load) == n_features) raw_load <- t(raw_load)
          
          if (nrow(raw_load) == n_features) {
            rownames(raw_load) <- rownames(seu)
            colnames(raw_load) <- paste0(seurat_key, 1:ncol(raw_load))
            
            # Inject Loadings
            Seurat::Loadings(dr_obj) <- raw_load
          }
        }
      }
      
      seu[[clean_name]] <- dr_obj
    }
  }
  return(seu)
}

#' Read and inject additional Layers (data, scale.data)
add_layers <- function(seu, h5_file, transpose_mode) {
  if (!h5_file$exists("layers")) return(seu)
  
  grp <- h5_file[["layers"]]
  layer_names <- grp$ls(recursive=FALSE)$name
  
  default_assay <- Seurat::DefaultAssay(seu)
  
  for (lname in layer_names) {
    # Decide if transposition is needed
    # If Standard Mode (obs=cells), layers are usually Cells x Genes. Seurat Layer needs Genes x Cells.
    # So in Standard Mode, read_matrix_node needs transpose_output=FALSE (because internal logic handles CSR transposition).
    # Or rely on read_matrix_node's smart logic?
    # Simply put: We want Genes x Cells.
    
    # The logic here must remain consistent with the Main Matrix reading in convert_h5ad.R.
    # If Main Matrix was transposed, Layer must also be transposed.
    
    # The transpose_output parameter here means: do we need to Transpose again on top of the read result?
    # If the file itself is Standard (Cells x Genes), reading CSR automatically becomes (Genes x Cells), we don't need to T again.
    # If the file is Transposed (Genes x Cells), reading CSR automatically becomes (Cells x Genes), we need T to get Genes x Cells.
    
    do_transpose <- transpose_mode 
    
    mat <- read_matrix_node(grp[[lname]], transpose_output = do_transpose)
    
    # Secondary dimension check
    if (ncol(mat) != ncol(seu) || nrow(mat) != nrow(seu)) {
      # If dimensions don't match, try transposing
      mat <- Matrix::t(mat)
    }
    
    if (ncol(mat) == ncol(seu) && nrow(mat) == nrow(seu)) {
      # Naming mapping
      # Scanpy: 'logcounts' -> Seurat: 'data'
      # Scanpy: 'scale' -> Seurat: 'scale.data'
      
      target_slot <- lname
      if (lname %in% c("logcounts", "norm_data", "normalize")) target_slot <- "data"
      if (grepl("scale", lname)) target_slot <- "scale.data"
      
      # Inject
      if (target_slot == "scale.data") {
        # scale.data must be a dense matrix
        seu <- Seurat::SetAssayData(seu, slot = "scale.data", new.data = as.matrix(mat), assay = default_assay)
      } else if (target_slot == "data") {
        seu <- Seurat::SetAssayData(seu, slot = "data", new.data = mat, assay = default_assay)
      } else {
        # Other layers, e.g., 'spliced'
        Seurat::LayerData(seu, assay = default_assay, layer = lname) <- mat
      }
      message(paste("    Layer added:", lname, "->", target_slot))
    }
  }
  return(seu)
}
