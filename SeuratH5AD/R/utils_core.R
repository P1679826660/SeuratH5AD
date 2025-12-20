#' @importFrom Matrix sparseMatrix t
#' @importFrom hdf5r h5attr_names h5attr
NULL

read_matrix_node <- function(node, transpose_output = FALSE) {
  mat <- NULL

  if (is.null(node) || !node$is_valid) return(NULL)

  if (inherits(node, 'H5Group')) {
    if (!exists('h5attr_names', where = asNamespace('hdf5r'))) return(NULL)
    attrs <- hdf5r::h5attr_names(node)

    encoding_attr <- if ('encoding-type' %in% attrs) hdf5r::h5attr(node, 'encoding-type') else 'csr_matrix'
    shape_attr <- if ('shape' %in% attrs) hdf5r::h5attr(node, 'shape') else NULL

    data <- node[['data']]$read()
    indices <- node[['indices']]$read()
    indptr <- node[['indptr']]$read()

    # Convert 64-bit integers to standard numeric for Matrix package compatibility
    if (inherits(indices, "integer64")) indices <- as.numeric(indices)
    if (inherits(indptr, "integer64")) indptr <- as.numeric(indptr)

    if (encoding_attr == 'csr_matrix') {
      # CSR logic: Standard H5AD format
      # CSR stores row pointers (indptr) and column indices (indices).
      # Matrix::sparseMatrix constructs columns from 'p'.
      # Passing CSR 'indptr' to 'p' effectively creates a Transposed matrix (CSC).
      # We construct the transposed matrix directly, then transpose it back to original orientation.
      mat <- Matrix::sparseMatrix(j = indices, p = indptr, x = data,
                                  dims = c(shape_attr[2], shape_attr[1]),
                                  index1 = FALSE)
      mat <- Matrix::t(mat) 
      if (transpose_output) mat <- Matrix::t(mat)

    } else if (encoding_attr == 'csc_matrix') {
      # CSC logic: Direct mapping
      mat <- Matrix::sparseMatrix(i = indices, p = indptr, x = data,
                                  dims = c(shape_attr[1], shape_attr[2]),
                                  index1 = FALSE)
      if (transpose_output) mat <- Matrix::t(mat)
    }

  } else if (inherits(node, 'H5D')) {
    # Dense matrix logic
    raw_data <- node$read()
    if (is.null(dim(raw_data))) raw_data <- as.matrix(raw_data)
    
    mat <- methods::as(raw_data, 'CsparseMatrix')
    if (transpose_output) mat <- Matrix::t(mat)
  }

  return(mat)
}

h5_node_attributes <- function(node) {
  if (!node$is_valid) return(list())
  n <- hdf5r::h5attr_names(node)
  res <- list()
  for (i in n) {
    res[[i]] <- tryCatch(hdf5r::h5attr(node, i), error = function(e) NULL)
  }
  return(res)
}
