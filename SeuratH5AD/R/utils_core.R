#' @importFrom Matrix sparseMatrix t
#' @importFrom hdf5r h5attr_names h5attr
#' @importFrom methods as
NULL

read_matrix_node <- function(node, transpose_output = FALSE) {
  mat <- NULL

  # Validation: Check if node exists and is valid
  if (is.null(node) || !node$is_valid) return(NULL)

  if (inherits(node, 'H5Group')) {
    if (!exists('h5attr_names', where = asNamespace('hdf5r'))) return(NULL)
    attrs <- hdf5r::h5attr_names(node)

    encoding_attr <- if ('encoding-type' %in% attrs) hdf5r::h5attr(node, 'encoding-type') else 'csr_matrix'
    shape_attr <- if ('shape' %in% attrs) hdf5r::h5attr(node, 'shape') else NULL

    # Safe read of sparse components
    data <- tryCatch(node[['data']]$read(), error = function(e) numeric(0))
    indices <- tryCatch(node[['indices']]$read(), error = function(e) integer(0))
    indptr <- tryCatch(node[['indptr']]$read(), error = function(e) integer(0))

    # Compatibility: Convert 64-bit integers to standard numeric
    if (inherits(indices, "integer64")) indices <- as.numeric(indices)
    if (inherits(indptr, "integer64")) indptr <- as.numeric(indptr)

    # Handle empty matrices
    if (length(data) == 0) {
       if (!is.null(shape_attr) && prod(shape_attr) > 0) {
         return(Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(shape_attr[1], shape_attr[2])))
       }
       return(NULL)
    }

    if (encoding_attr == 'csr_matrix') {
      # CSR Logic:
      # H5AD CSR stores (rows, cols). 'indptr' points to row starts. 'indices' are column indices.
      # To load efficiently into R (CSC default), we construct the TRANSPOSE.
      # The transpose of a CSR matrix is a CSC matrix where:
      #   - Original Row Pointers (indptr) -> Become Column Pointers (p)
      #   - Original Col Indices (indices) -> Become Row Indices (i)
      #   - Dimensions are swapped
      
      mat <- Matrix::sparseMatrix(i = indices, p = indptr, x = data,
                                  dims = c(shape_attr[2], shape_attr[1]),
                                  index1 = FALSE)
      
      # Transpose back to original orientation (Rows x Cols)
      mat <- Matrix::t(mat)
      
      if (transpose_output) mat <- Matrix::t(mat)

    } else if (encoding_attr == 'csc_matrix') {
      # CSC Logic: Direct mapping
      mat <- Matrix::sparseMatrix(i = indices, p = indptr, x = data,
                                  dims = c(shape_attr[1], shape_attr[2]),
                                  index1 = FALSE)
      if (transpose_output) mat <- Matrix::t(mat)
    }

  } else if (inherits(node, 'H5D')) {
    # Dense Matrix Logic
    raw_data <- tryCatch(node$read(), error = function(e) NULL)
    if (is.null(raw_data)) return(NULL)
    
    if (is.null(dim(raw_data))) raw_data <- as.matrix(raw_data)
    
    # Force conversion to sparse to optimize memory usage in R
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
