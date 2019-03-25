#' Initialize and setup the Signac object
#'
#' Initializes the Signac object and some optional filtering
#' @param raw.data Raw input data
#' @param type Type of the input file "mtx", "h5", "csv", "tsv", "sparse",
#' "matrix"
#' @param project Project name (string)
#' @param ... parameters for h5 file, tsv file,
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
                         `csv` = function(path, ...) ReadDelim(path, sep = ",", ...),
                         `tsv` = function(path, ...) ReadDelim(path, sep = "\t", ...),
                         `h5`  = function(path, ...) Read10XH5(path, ...)
                         )

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
#' @importFrom Matrix readMM
#' @importFrom Matrix Matrix
#' @importFrom readr read_lines read_delim
#'
ReadDelim <- function(mat.path, sep = ",", header = TRUE) {

  # Get header. Handle a case when the 1st element is missing
  getHeader <- function(mat.path, sep) {
    lines <- read_lines(mat.path, n_max = 2)
    lines <- lapply(lines, getLine, sep = sep)
    if (length(lines[[1]]) == length(lines[[2]])) {
      lines[[1]] <- lines[[1]][-1]
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
    line1 <- getLine(read_lines(mat.path, n_max = 1), sep)
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
  mat <- suppressWarnings(read_delim(
    mat.path,
    delim = sep,
    skip = n.skip,
    col_names = header.arr,
    progress = FALSE
  ))
  mat <- convertTibbleToSparseMatrix(mat)
  return(mat)
}

#' This function was initialized by the idea
#' from Seurat R package
#' citation: Butler et al., Nature Biotechnology 2018
#'
#' Read 10X hdf5 file
#'
#' Read count matrix from 10X CellRanger hdf5 file.
#' This can be used to read both scATAC-seq and scRNA-seq matrices.
#'
#' @param filename Path to h5 file
#' @param use.names Label row names with feature names rather than ID numbers.
#' @param unique.features Make feature names unique (default TRUE)
#'
#' @return Returns a sparse matrix with rows and columns labeled. If multiple
#' genomes are present, returns a list of sparse matrices (one per genome).
#'
#' @importFrom hdf5r H5File
#' @importFrom Matrix sparseMatrix
#'
Read10XH5 <- function(filename, use.names = TRUE, unique.features = TRUE) {

  # Get version of HDF5 file
  GetVersion <- function(infile) {
    return(if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) 3 else 2)
  }

  # Get row names (for genes)
  GetFeatureSlot <- function(version, use.names) {
    slot <- if (version == 3) {
      if (use.names) "features/name" else "features/id"
    } else {
      if (use.names) "gene_names" else "genes"
    }
  }

  # Get features
  GetFeatures <- function(infile, base.slot, feature.slot, unique.features) {
    features <- infile[[paste0(base.slot, '/', feature.slot)]][]
    if (unique.features) {
      features <- make.unique(names = features)
    }
    return(features)
  }

  # Get data from v2 H5 file
  GetOutputV2 <- function(infile, feature.slot, unique.features) {
    output <- list()
    genomes <- names(x = infile)
    for (genome in genomes) {
      counts <- infile[[paste0(genome, '/data')]]
      indices <- infile[[paste0(genome, '/indices')]]
      indptr <- infile[[paste0(genome, '/indptr')]]
      shp <- infile[[paste0(genome, '/shape')]]
      features <- GetFeatures(infile, genome, feature.slot, unique.features)
      barcodes <- infile[[paste0(genome, '/barcodes')]]
      sparse.mat <- sparseMatrix(
        i = indices[] + 1,
        p = indptr[],
        x = as.numeric(x = counts[]),
        dims = shp[],
        giveCsparse = FALSE
      )
      rownames(x = sparse.mat) <- features
      colnames(x = sparse.mat) <- barcodes[]
      sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
      output[[genome]] <- sparse.mat
    }
    infile$close_all()
    return(output)
  }

  # Get data from v3 H5 file
  GetOutputV3 <- function(infile, feature.slot, unique.features) {

    # Get genomes and indices (of features)
    GetGenomes <- function(infile, base.slot) {
      genome.arr <- infile[[paste0(base.slot, "/features/genome")]][]
      name <- unique(genome.arr)
      genome.arr <- lapply(name, function(x) which(genome.arr == x))
      names(genome.arr) <- name
      return(genome.arr)
    }

    # Get feature type
    GetFeatureType <- function(infile, base.slot) {
      return (if (infile$exists(name = paste0(base.slot, '/features/feature_type'))) {
        infile[[paste0(base.slot, '/features/feature_type')]][]
      } else {
        NULL
      })
    }

    output <- list()
    base.slot <- names(x = infile)
    genomes <- GetGenomes(infile, base.slot)
    counts <- infile[[paste0(base.slot, '/data')]]
    indices <- infile[[paste0(base.slot, '/indices')]]
    indptr <- infile[[paste0(base.slot, '/indptr')]]
    shp <- infile[[paste0(base.slot, '/shape')]]
    features <- GetFeatures(infile, base.slot, feature.slot, unique.features)
    types <- GetFeatureType(infile, base.slot)
    barcodes <- infile[[paste0(base.slot, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      giveCsparse = FALSE
    )
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    types.unique <- unique(types)
    output <- lapply(names(genomes), function(genome) {
      ft.idx <- genomes[[genome]]
      sparse.mat2 <- sparse.mat[ft.idx, ]
      types2 <- types[ft.idx]
      sparse.mat2 <- lapply(types.unique, function(x) sparse.mat2[which(x = types2 == x), ])
      names(sparse.mat2) <- types.unique
      return(sparse.mat2)
    })
    names(output) <- names(genomes)
    infile$close_all()
    return(output)
  }

  infile <- H5File$new(filename = filename, mode = 'r')
  version <- GetVersion(infile)
  feature.slot <- GetFeatureSlot(version, use.names)
  GetOutputFunc <- if (version == 2) GetOuputV2 else GetOutputV3
  return(GetOutputFunc(infile, feature.slot, unique.features))
}
