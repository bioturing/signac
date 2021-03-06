#' Check whether a file is exists, then stop with a message
#' @param path path to  a file on your system
CheckFileExist <- function(path) {
  if (!file.exists(path)) {
    stop(paste("File", path, "does not exists!"))
  }
}

#' This function was imported from Seurat R package
#' citation: Butler et al., Nature Biotechnology 2018
#' Extract delimiter information from a string.
#'
#' Parses a string (usually a cell name) and extracts fields based on a delimiter
#'
#' @param string String to parse.
#' @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
#' @param delim Delimiter to use, set to underscore by default.
#'
#' @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#'
#' @export
#'
#' @examples
#' ExtractField(string = 'Hello World', field = 1, delim = '_')
#'
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

#' Colourise text for display in the terminal.
#'
#' If R is not currently running in a system that supports terminal colours
#' the text will be returned unchanged.
#' This function was imported from testthat package
#' https://github.com/r-lib/testthat
#' Allowed colours are: black, blue, brown, cyan, dark gray, green, light
#' blue, light cyan, light gray, light green, light purple, light red,
#' purple, red, white, yellow
#'
#' @param text character vector
#' @param fg foreground colour, defaults to white
#' @param bg background colour, defaults to transparent
#' @examples
#' print(Colourise("Red", "red"))
#' cat(Colourise("Red", "red"), "\n")
#' cat(Colourise("White on red", "white", "red"), "\n")
Colourise <- function(text, fg = "black", bg = NULL) {
  term <- Sys.getenv()["TERM"]
  colour_terms <- c("xterm-color","xterm-256color", "screen", "screen-256color")

  if(rcmd_running() || !any(term %in% colour_terms, na.rm = TRUE)) {
    return(text)
  }

  col_escape <- function(col) {
    paste0("\033[", col, "m")
  }

  col <- .fg_colours[tolower(fg)]
  if (!is.null(bg)) {
    col <- paste0(col, .bg_colours[tolower(bg)], sep = ";")
  }

  init <- col_escape(col)
  reset <- col_escape("0")
  paste0(init, text, reset)
}

.fg_colours <- c(
  "black" = "0;30",
  "blue" = "0;34",
  "green" = "0;32",
  "cyan" = "0;36",
  "red" = "0;31",
  "purple" = "0;35",
  "brown" = "0;33",
  "light gray" = "0;37",
  "dark gray" = "1;30",
  "light blue" = "1;34",
  "light green" = "1;32",
  "light cyan" = "1;36",
  "light red" = "1;31",
  "light purple" = "1;35",
  "yellow" = "1;33",
  "white" = "1;37"
)

.bg_colours <- c(
  "black" = "40",
  "red" = "41",
  "green" = "42",
  "brown" = "43",
  "blue" = "44",
  "purple" = "45",
  "cyan" = "46",
  "light gray" = "47"
)

rcmd_running <- function() {
  nchar(Sys.getenv('R_TESTS')) != 0
}

#' Flexible check file exists in an 10X folder
#' @param dir to the a 10X folder
#' @return a list of three file names: barcode, genes, and matrix
#'
ReturnTripleFiles <- function(dir) {
  barcodes <- list("barcodes.tsv", "barcodes.tsv.gz")
  genes <- list("genes.tsv", "features.tsv.gz")
  matrix <- list("matrix.mtx", "matrix.mtx.gz")

  bx.exists <- sapply(file.path(dir, barcodes), file.exists)
  genes.exists <- sapply(file.path(dir, genes), file.exists)
  matrix.exists <- sapply(file.path(dir, matrix), file.exists)

  if (!any(bx.exists)) {
    stop(paste("Both", barcodes[[1]], "and", barcodes[[2]], "do not exists in",
               dir, sep = " "))
  }

  if (!any(genes.exists)) {
    stop(paste("Both", genes[[1]], "and", genes[[2]], "do not exists in",
               dir, sep = " "))
  }

  if (!any(matrix.exists)) {
    stop(paste("Both", matrix[[1]], "and", matrix[[2]], "do not exists in",
               dir, sep = " "))
  }

  return(list(genes    = genes[genes.exists][[1]],
              barcodes = barcodes[bx.exists][[1]],
              matrix   = matrix[matrix.exists][[1]]))
}


#' Ask for a choice. Use to select one genome from multi-species data
#'
#' @param choices : multi genomes selection
#'
#' @return selected genome
#'
AskForChoices <- function(choices) {
  n <- length(choices)
  choices.txt <- paste(seq(1, n), choices, sep = ": ")
  choices.txt <- Colourise(paste(c(choices.txt, ""), collapse = "\n"), "green")

  opts <- as.character(seq(1, n))
  select <- NA
  while (! select %in% opts) {
    message("We only support single-specie analysis. Please select from ", 1, " to ", n)
    select <- readline(choices.txt)
  }
  return(choices[as.numeric(select)])
}

#' Write Sparse Matrix using rhdf5 api
#'
#' @export
#'
WriteSpMtAsS4 <- function(filePath, groupName, mat) {
  if (!file.exists(filePath)) {
    rhdf5::h5createFile(filePath)
  }
  rhdf5::h5createGroup(filePath, groupName)
  rhdf5::h5createDataset(filePath, paste0(groupName, "/indices"), storage.mode = "integer",
                         dims = length(mat@i), chunk=min(10000, length(mat@i)), level=1)
  rhdf5::h5write(mat@i, filePath, paste0(groupName, "/indices"))
  rhdf5::h5createDataset(filePath, paste0(groupName, "/data"), storage.mode = "double",
                         dims = length(mat@x), chunk=min(10000, length(mat@x)), level=1)
  rhdf5::h5write(mat@x, filePath, paste0(groupName, "/data"))
  rhdf5::h5write(mat@p, filePath, paste0(groupName, "/indptr"))
  rhdf5::h5write(rownames(mat), filePath, paste0(groupName, "/features"))
  rhdf5::h5write(colnames(mat), filePath, paste0(groupName, "/barcodes"))
  rhdf5::h5write(dim(mat), filePath, paste0(groupName, "/shape"))
}
