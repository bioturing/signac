#' Initialize and setup the Signac object
#'
#' Initializes the Signac object and some optional filtering
#' @param raw.data Raw input data
#' @param type Type of the input file "mtx", "h5", "csv", "tsv", "sparse",
#' "matrix"
#' @param project Project name (string)
#' @param ... parameters for h5 file
#' @return Returns a Signac object with the raw data stored in object@@raw.data.
#' object@@data, object@@meta.data, object@@ident, also initialized.
#'
#' @importFrom methods new
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @examples
#' pbmc_small <- CreateSignacObject(raw.data = path2pbmc, type = "mtx")
#' pbmc_small
#'
CreateSignacObject <- function(
  raw.data,
  type = "mtx",
  project = "UntitleProject",
  ...
) {

  if (type %in% c("mtx", "h5", "csv", "tsv")) {
    CheckFileExist(raw.data)
  }

  if (type == "matrix") {
    stopifnot("matrix" %in% class(raw.data))
  }

  if (type == "sparse") {
    stopifnot("sparseMatrix" %in% class(raw.data))
  }

  ReadFunction <- switch(type,
                         `mtx` = Read10X,
                         `csv` = function(path, ...) ReadMatrix(path, sep = ","),
                         `tsv` = function(path, ...) ReadMatrix(path, sep = "\t"),
                         `h5`  = NULL)
  data <- ReadFunction(raw.data, ...)
  signac.version <- packageVersion("Signac")
  object <- new(
    Class = "Signac",
    raw.data = data,
    project.name = project,
    version = signac.version
  )
  return(object)
}

#' This function was imported from Seurat R package
#' citation: Butler et al., Nature Biotechnology 2018
#'
#' Load in data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load
#' several data directories. If a named vector is given, the cell barcode names
#' will be prefixed with the name.
#'
#' @return Returns a sparse matrix with rows and columns labeled
#'
#' @importFrom Matrix readMM
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Read10X(data.dir = data_dir)
#' signac_object = CreateSignacObject(raw.data = expression_matrix)
#' }
#'
Read10X <- function(data.dir = NULL){
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (! dir.exists(run)) {
      stop("Directory provided does not exist")
    }
    if(!grepl("\\/$", run)) {
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv")
    gene.loc <- paste0(run, "genes.tsv")
    matrix.loc <- paste0(run, "matrix.mtx")
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc)) {
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc)) {
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(
        x = as.character(
          x = sapply(
            X = cell.names,
            FUN = ExtractField,
            field = 1,
            delim = "-"
          )
        )
      )
    }
    rownames(x = data) <- make.unique(
      names = as.character(
        x = sapply(
          X = gene.names,
          FUN = ExtractField,
          field = 2,
          delim = "\\t"
        )
      )
    )
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}

#' Read expresison matrix from CSV or TSV
#' @param mat.path Path to expression matrix (can be zipped)
#' @param sep Separator. Default is <code>,</code>
#' @param sep Data has a header. Default is TRUE
ReadMatrix <- function(mat.path, sep = ",", header = TRUE) {

  # Get header. Handle a case when the 1st element is missing
  getHeader <- function(mat.path, sep) {
    lines <- readr::read_lines(mat.path, n_max = 2)
    lines <- lapply(lines, getLine, sep = sep)
    if (length(lines[[1]]) == length(lines[[2]])) {
      lines[[1]] <- line[[1]][-1]
    }
    lines[[1]] <- c('gene', lines[[1]])
    return(lines[[1]])
  }

  # Get line's components
  getLine <- function(line, sep) {
    return(strsplit(line, sep)[[1]])
  }

  # Fake header
  fakeHeader <- function(mat.path, sep) {
    line1 <- getLine(readr::read_lines(mat.path, n_max = 1), sep)
    return(paste0('c', seq_along(line1)))
  }

  # convert tibble table to a sparse matrix
  convertTibbleToSparseMatrix <- function(mat) {
    rn <- mat$gene
    mat <- Matrix(as.matrix(mat[, -1]), sparse = TRUE)
    rownames(mat) <- rn
    return(mat)
  }

  header.arr <- if (header) getHeader(mat.path, sep) else fakeHeader(mat.path, sep) 
  n.skip <- if (header) 1 else 0
  mat <- suppressWarnings(readr::read_delim(
    mat.path,
    delim = sep,
    skip = n.skip,
    col_names = header.arr,
    progress = FALSE
  ))
  mat <- convertTibbleToSparseMatrix(mat)
  return(mat)
}
