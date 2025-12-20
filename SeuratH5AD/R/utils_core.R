#' @importFrom Matrix sparseMatrix t
#' @importFrom hdf5r h5attr_names h5attr
NULL

read_matrix_node <- function(node, transpose_output = FALSE) {
  mat <- NULL

  # Validity check
  if (is.null(node) || !node$is_valid) return(NULL)

  if (inherits(node, 'H5Group')) {
    if (!exists('h5attr_names', where = asNamespace('hdf5r'))) return(NULL)
    attrs <- hdf5r::h5attr_names(node)

    encoding_attr <- if ('encoding-type' %in% attrs) hdf5r::h5attr(node, 'encoding-type') else 'csr_matrix'
    shape_attr <- if ('shape' %in% attrs) hdf5r::h5attr(node, 'shape') else NULL

    # Safe read with tryCatch
    data <- tryCatch(node[['data']]$read(), error = function(e) numeric(0))
    indices <- tryCatch(node[['indices']]$read(), error = function(e) integer(0))
    indptr <- tryCatch(node[['indptr']]$read(), error = function(e) integer(0))

    # Convert integer64 to numeric (standard R integer/numeric)
    if (inherits(indices, "integer64")) indices <- as.numeric(indices)
    if (inherits(indptr, "integer64")) indptr <- as.numeric(indptr)

    # If empty data, return NULL or empty matrix logic
    if (length(data) == 0 && !is.null(shape_attr) && prod(shape_attr) > 0) {
      # Fallback for empty sparse matrix but with shape
      return(Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(shape_attr[1], shape_attr[2])))
    }
    if (length(data) == 0) return(NULL)

    if (encoding_attr == 'csr_matrix') {
      # [CRITICAL FIX] CSR Loading Logic
      # CSR: indptr=row_start, indices=col_idx.
      # To load into R (CSC) efficiently, we construct the TRANSPOSE.
      # Transpose of CSR (M x N) is CSC (N x M) where:
      #   - CSR indptr (Row Ptrs) -> CSC p (Col Ptrs)
      #   - CSR indices (Col Idx) -> CSC i (Row Idx)
      #   - Dimensions are swapped: c(N, M)
      
      mat <- Matrix::sparseMatrix(i = indices, p = indptr, x = data,
                                  dims = c(shape_attr[2], shape_attr[1]),
                                  index1 = FALSE)
      
      # Now transpose back to get original orientation (M x N)
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
    # Dense matrix logic
    raw_data <- tryCatch(node$read(), error = function(e) NULL)
    if (is.null(raw_data)) return(NULL)
    
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
