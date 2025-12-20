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

    # Read raw vectors
    # [Performance Note] For extremely large matrices, we might need chunked reading,
    # but for now we read fully to construct the sparse matrix.
    data <- tryCatch(node[['data']]$read(), error = function(e) numeric(0))
    indices <- tryCatch(node[['indices']]$read(), error = function(e) integer(0))
    indptr <- tryCatch(node[['indptr']]$read(), error = function(e) integer(0))

    # Convert 64-bit integers to standard numeric for Matrix package compatibility
    # HDF5 often returns integer64 which Matrix doesn't like
    if (inherits(indices, "integer64")) indices <- as.numeric(indices)
    if (inherits(indptr, "integer64")) indptr <- as.numeric(indptr)

    # Validate inputs to avoid cryptic Matrix errors
    if (length(data) == 0 || length(indices) == 0 || length(indptr) == 0) {
      return(NULL) 
    }

    if (encoding_attr == 'csr_matrix') {
      # [CRITICAL FIX] CSR Logic
      # CSR stores (rows, cols). 'indptr' points to row starts. 'indices' are column indices.
      # To load into R (CSC default), we simulate loading the TRANSPOSE of this matrix.
      # The Transpose of CSR is a CSC matrix where:
      #   - Original Row Pointers (indptr) -> Become Column Pointers (p)
      #   - Original Col Indices (indices) -> Become Row Indices (i)  <-- Was 'j' in error
      #   - Dims are swapped: (Cols, Rows)
      
      mat <- Matrix::sparseMatrix(i = indices, p = indptr, x = data,
                                  dims = c(shape_attr[2], shape_attr[1]),
                                  index1 = FALSE)
      
      # Now mat is A_transpose. We need A.
      # Transposing it back gives us the original matrix in correct R format.
      mat <- Matrix::t(mat) 
      
      if (transpose_output) mat <- Matrix::t(mat)

    } else if (encoding_attr == 'csc_matrix') {
      # CSC Logic: Direct mapping
      # 'indptr' is column pointers (p)
      # 'indices' is row indices (i)
      mat <- Matrix::sparseMatrix(i = indices, p = indptr, x = data,
                                  dims = c(shape_attr[1], shape_attr[2]),
                                  index1 = FALSE)
      if (transpose_output) mat <- Matrix::t(mat)
    }

  } else if (inherits(node, 'H5D')) {
    # Dense matrix logic
    raw_data <- node$read()
    if (is.null(dim(raw_data))) raw_data <- as.matrix(raw_data)
    
    # Cast to sparse to save memory in R
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
