# ==============================================================================
# Module Name: H5AD Matrix Loader (V2 Generic)
# Function Description: Generic HDF5 matrix loader, supports X, raw/X, layers/*
# ==============================================================================

#' Generic matrix reading function
#' @param node H5Group or H5Dataset (e.g., file[["X"]] or file[["layers"]][["norm_data"]])
#' @param transpose_output Logical. If TRUE, transpose the result (to fit Seurat Genes x Cells)
#' @return Matrix object (CsparseMatrix or Dense)
read_matrix_node <- function(node, transpose_output = FALSE) {
  
  mat <- NULL
  
  # 1. Identify storage type
  if (inherits(node, "H5Group")) {
    # === Sparse Matrix (CSR/CSC) ===
    attrs <- hdf5r::h5attr_names(node)
    encoding_attr <- if("encoding-type" %in% attrs) hdf5r::h5attr(node, "encoding-type") else "csr_matrix"
    shape_attr <- if("shape" %in% attrs) hdf5r::h5attr(node, "shape") else NULL
    
    data <- node[["data"]]$read()
    indices <- node[["indices"]]$read()
    indptr <- node[["indptr"]]$read()
    
    # Assume it is CSR (Python default)
    # Python CSR: (Cells, Genes). indices=cols, indptr=rows
    # R sparseMatrix default construction is CSC.
    # If we feed CSR data into the CSC constructor:
    #   i (indices) -> becomes row indices
    #   p (indptr) -> becomes column pointers
    #   Result: The matrix is physically transposed. i.e., (Cells x Genes) -> (Genes x Cells)
    
    if (encoding_attr == "csr_matrix") {
      # At this point, mat is already Genes x Cells (if original data was Cells x Genes)
      mat <- Matrix::sparseMatrix(
        i = indices + 1,
        p = indptr,
        x = data,
        dims = c(shape_attr[2], shape_attr[1]), # Reverse dimensions
        index1 = TRUE
      )
      
      # If it is already transposed here (became Genes x Cells), and the user also requests transpose_output
      # The logic here is a bit tricky.
      # Current status: CSR read into Seurat is "naturally transposed".
      # If transpose_output = TRUE (meaning we need Genes x Cells), then no action is needed here.
      # If transpose_output = FALSE (meaning we need to keep Cells x Genes), we need to manually t() it back.
      
      if (!transpose_output) {
        mat <- Matrix::t(mat)
      }
      
    } else if (encoding_attr == "csc_matrix") {
      # CSC is read in as-is (Cells x Genes)
      mat <- Matrix::sparseMatrix(
        i = indices + 1,
        p = indptr,
        x = data,
        dims = c(shape_attr[1], shape_attr[2]),
        index1 = TRUE
      )
      if (transpose_output) {
        mat <- Matrix::t(mat)
      }
    }
    
  } else {
    # === Dense Matrix (Dataset) ===
    # Reading HDF5 Dataset into R preserves dimension order, but R is column-major.
    # Usually reads out as (Genes, Cells) if originally stored as (Cells, Genes) with no special attributes.
    # To be safe, read as standard matrix first
    raw_data <- node[,]
    
    # Convert to sparse for unified handling (unless dense is specifically required, convert here first)
    mat <- as(raw_data, "CsparseMatrix")
    
    # If user requests transposition
    if (transpose_output) {
      mat <- Matrix::t(mat)
    }
  }
  
  return(mat)
}

# Helper: Get node attributes
h5_node_attributes <- function(node) {
  if(!exists("h5attr_names", where = asNamespace("hdf5r"))) return(list())
  
  n <- hdf5r::h5attr_names(node)
  res <- list()
  for (i in n) res[[i]] <- hdf5r::h5attr(node, i)
  return(res)
}
