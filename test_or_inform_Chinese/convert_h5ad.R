# ==============================================================================
# Module Name: H5AD Converter Main (V6 Full Stack)
# Function Description: Final version of the scheduler. Handle Raw/X relationship, assemble Layers, Reductions, Meta.features
# ==============================================================================

#' Main function: convert H5AD to Seurat V5 Object (Full)
#' @export
h5ad_to_seurat <- function(file_path) {
  message(">>> [Phase 1] Initialization of HDF5: ", file_path)
  file_h5 <- hdf5r::H5File$new(file_path, mode = "r")
  on.exit(file_h5$close_all(), add = TRUE)
  
  # === Phase 2: Metadata reading and direction detection ===
  message(">>> [Phase 2] Detecting data direction (DNA Check)...")
  obs <- read_h5ad_dataframe(file_h5, "obs")
  var <- read_h5ad_dataframe(file_h5, "var")
  
  # Simplified DNA detection logic (based on previous V5)
  obs_names <- if(!is.null(obs)) rownames(obs) else character(0)
  var_names <- if(!is.null(var)) rownames(var) else character(0)
  
  # If obs names look like genes (non-DNA) and var looks like DNA -> Transpose Mode
  obs_dna <- calculate_dna_score(obs_names)
  var_dna <- calculate_dna_score(var_names)
  
  transpose_mode <- FALSE
  if (var_dna > 0.8 && obs_dna < 0.5) {
    message("!!! Transpose mode detected (obs=Genes, var=Cells)")
    transpose_mode <- TRUE
  }
  
  # Determine final metadata
  final_cells_meta <- if(transpose_mode) var else obs
  final_features_meta <- if(transpose_mode) obs else var
  
  # === Phase 3: Determine Counts Source (Raw vs X) ===
  # Logic: If file[['raw/X']] exists, it is usually counts. Then file[['X']] is usually data.
  # If no raw, file[['X']] is counts.
  
  counts_node <- NULL
  data_node <- NULL # For storing normalized data
  
  has_raw <- file_h5$exists("raw") && file_h5[["raw"]]$exists("X")
  
  if (has_raw) {
    message(">>> 'raw' node found, using raw/X as Counts, root X as Data...")
    counts_node <- file_h5[["raw/X"]]
    # Note: If raw is used, features metadata is usually in raw/var
    raw_var <- read_h5ad_dataframe(file_h5, "raw/var")
    if (!is.null(raw_var) && !transpose_mode) final_features_meta <- raw_var
    
    # Root X saved as data (if it exists)
    if (file_h5$exists("X")) data_node <- file_h5[["X"]]
    
  } else {
    message(">>> 'raw' not found, using root X as Counts...")
    if (file_h5$exists("X")) counts_node <- file_h5[["X"]]
  }
  
  if (is.null(counts_node)) stop("Could not find expression matrix (X or raw/X)")
  
  # === Phase 4: Read Main Matrix and Construct Object ===
  message(">>> [Phase 4] Reading Counts Matrix...")
  # Logic for transpose_output here:
  # If transpose_mode=TRUE (obs is genes), it means the file contains (Genes x Cells).
  # read_matrix_node transposes CSR by default. If reading CSR -> becomes (Cells x Genes).
  # We need (Genes x Cells). So we need transpose_output = TRUE.
  # The logic here is tricky; the best way is: read it, check dimensions, transpose if incorrect.
  
  counts_mat <- read_matrix_node(counts_node, transpose_output = FALSE)
  
  # Force dimension alignment
  target_cells <- nrow(final_cells_meta)
  target_feats <- nrow(final_features_meta)
  
  if (ncol(counts_mat) != target_cells) {
    counts_mat <- Matrix::t(counts_mat)
  }
  
  # Assign names
  if(!is.null(final_cells_meta)) colnames(counts_mat) <- rownames(final_cells_meta)
  if(!is.null(final_features_meta)) rownames(counts_mat) <- rownames(final_features_meta)
  
  message(">>> Constructing base Seurat Object...")
  seu <- Seurat::CreateSeuratObject(counts = counts_mat, meta.data = final_cells_meta, project = "H5AD")
  
  # === Phase 5: Inject Data (if available) ===
  if (!is.null(data_node)) {
    message(">>> [Phase 5] Injecting Normalized Data...")
    data_mat <- read_matrix_node(data_node, transpose_output = FALSE)
    if (ncol(data_mat) != target_cells) data_mat <- Matrix::t(data_mat)
    
    # Assign names for safety
    colnames(data_mat) <- colnames(seu)
    rownames(data_mat) <- rownames(seu)
    
    seu <- Seurat::SetAssayData(seu, slot = "data", new.data = data_mat)
  }
  
  # === Phase 6: Inject Feature Metadata ===
  if (!is.null(final_features_meta)) {
    message(">>> [Phase 6] Injecting Feature Metadata...")
    default_assay <- Seurat::DefaultAssay(seu)
    for (i in colnames(final_features_meta)) {
      seu[[default_assay]][[i]] <- final_features_meta[[i]]
    }
    
    # Detect highly variable genes
    if ("highly_variable" %in% colnames(final_features_meta)) {
      hv_genes <- rownames(final_features_meta)[which(final_features_meta$highly_variable == TRUE)]
      if (length(hv_genes) > 0) {
        Seurat::VariableFeatures(seu) <- hv_genes
        message(paste("    Highly variable genes marked:", length(hv_genes)))
      }
    }
  }
  
  # === Phase 7: Inject other Layers, Reductions ===
  # Must source utils_h5_components.R first
  if (exists("add_reductions")) {
    message(">>> [Phase 7] Injecting Reductions & Layers...")
    seu <- add_layers(seu, file_h5, transpose_mode)
    seu <- add_reductions(seu, file_h5, transpose_mode)
  }
  
  message(">>> Conversion fully completed!")
  return(seu)
}
