#' @importFrom Matrix sparseMatrix t
#' @importFrom hdf5r h5attr_names h5attr
NULL

read_matrix_node <- function(node, transpose_output = FALSE) {
  mat <- NULL
  if (inherits(node, 'H5Group')) {
    if(!exists('h5attr_names', where = asNamespace('hdf5r'))) return(NULL)
    attrs <- hdf5r::h5attr_names(node)
    encoding_attr <- if('encoding-type' %in% attrs) hdf5r::h5attr(node, 'encoding-type') else 'csr_matrix'
    shape_attr <- if('shape' %in% attrs) hdf5r::h5attr(node, 'shape') else NULL
    data <- node[['data']]$read()
    indices <- node[['indices']]$read()
    indptr <- node[['indptr']]$read()
    if (encoding_attr == 'csr_matrix') {
      mat <- Matrix::sparseMatrix(i = indices + 1, p = indptr, x = data, dims = c(shape_attr[2], shape_attr[1]), index1 = TRUE)
      if (!transpose_output) mat <- Matrix::t(mat)
    } else if (encoding_attr == 'csc_matrix') {
      mat <- Matrix::sparseMatrix(i = indices + 1, p = indptr, x = data, dims = c(shape_attr[1], shape_attr[2]), index1 = TRUE)
      if (transpose_output) mat <- Matrix::t(mat)
    }
  } else {
    raw_data <- node[,]
    mat <- methods::as(raw_data, 'CsparseMatrix')
    if (transpose_output) mat <- Matrix::t(mat)
  }
  return(mat)
}

calculate_dna_score <- function(identifiers) {
  if (length(identifiers) == 0) return(0)
  n_sample <- min(length(identifiers), 100)
  sample_ids <- identifiers[1:n_sample]
  clean_ids <- toupper(gsub('[^a-zA-Z]', '', sample_ids))
  total_chars <- sum(nchar(clean_ids))
  if (total_chars == 0) return(0)
  non_acgt_chars <- gsub('[ACGT]', '', clean_ids)
  return(1 - (sum(nchar(non_acgt_chars)) / total_chars))
}

h5_node_attributes <- function(node) {
  n <- hdf5r::h5attr_names(node)
  res <- list(); for (i in n) res[[i]] <- hdf5r::h5attr(node, i); return(res)
}
