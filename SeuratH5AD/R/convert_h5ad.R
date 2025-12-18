#' Convert H5AD to Seurat V5
#' @param file_path Path to .h5ad file
#' @return Seurat object
#' @export
h5ad_to_seurat <- function(file_path) {
  message('>>> [Import] Loading: ', file_path)
  file_h5 <- hdf5r::H5File$new(file_path, mode = 'r')
  on.exit(file_h5$close_all(), add = TRUE)
  obs <- read_h5ad_dataframe(file_h5, 'obs')
  var <- read_h5ad_dataframe(file_h5, 'var')
  obs_dna <- calculate_dna_score(if(!is.null(obs)) rownames(obs) else character(0))
  var_dna <- calculate_dna_score(if(!is.null(var)) rownames(var) else character(0))
  transpose_mode <- FALSE
  if (var_dna > 0.8 && obs_dna < 0.5) {
    message('!!! Transpose detected (obs=Genes, var=Cells)')
    transpose_mode <- TRUE
  }
  final_cells_meta <- if(transpose_mode) var else obs
  final_features_meta <- if(transpose_mode) obs else var
  counts_node <- NULL; data_node <- NULL
  has_raw <- file_h5$exists('raw') && file_h5[['raw']]$exists('X')
  if (has_raw) {
    counts_node <- file_h5[['raw/X']]
    raw_var <- read_h5ad_dataframe(file_h5, 'raw/var')
    if (!is.null(raw_var) && !transpose_mode) final_features_meta <- raw_var
    if (file_h5$exists('X')) data_node <- file_h5[['X']]
  } else { if (file_h5$exists('X')) counts_node <- file_h5[['X']] }
  if (is.null(counts_node)) stop('No expression matrix found (X or raw/X)')
  counts_mat <- read_matrix_node(counts_node, transpose_output = FALSE)
  target_cells <- nrow(final_cells_meta)
  if (ncol(counts_mat) != target_cells) counts_mat <- Matrix::t(counts_mat)
  if(!is.null(final_cells_meta)) colnames(counts_mat) <- rownames(final_cells_meta)
  if(!is.null(final_features_meta)) rownames(counts_mat) <- rownames(final_features_meta)
  seu <- Seurat::CreateSeuratObject(counts = counts_mat, meta.data = final_cells_meta, project = 'H5AD')
  if (!is.null(data_node)) {
    data_mat <- read_matrix_node(data_node, transpose_output = FALSE)
    if (ncol(data_mat) != target_cells) data_mat <- Matrix::t(data_mat)
    colnames(data_mat) <- colnames(seu); rownames(data_mat) <- rownames(seu)
    seu <- Seurat::SetAssayData(seu, slot = 'data', new.data = data_mat)
  }
  if (!is.null(final_features_meta)) {
    da <- Seurat::DefaultAssay(seu)
    for (i in colnames(final_features_meta)) seu[[da]][[i]] <- final_features_meta[[i]]
    if ('highly_variable' %in% colnames(final_features_meta)) {
      hv <- rownames(final_features_meta)[which(final_features_meta$highly_variable == TRUE)]
      if(length(hv)>0) Seurat::VariableFeatures(seu) <- hv
    }
  }
  seu <- add_layers(seu, file_h5, transpose_mode)
  seu <- add_reductions(seu, file_h5, transpose_mode)
  message('>>> Conversion Complete!')
  return(seu)
}
