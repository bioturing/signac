

#' Filter epxression matrix
#'
#' @param object Signac object
#' @param slot Data for filtering. Default is raw.data
#' @param min.cells The minimum number of cells for a gene. Can be a numeric variable.
#'                  If NULL, it will be the at least 10 or the 2% of the total cells.
#' @param min.genes The minimum number of genes for a cell. Default is 200.
#' @param not.expressed Genes that have expression value equal or less than
#'                      this number is considered not expressed. Default is 0.
#' @param verbose Talkative or not
#' @importFrom Matrix rowSums colSums
FilterData <- function(
    object,
    slot = "raw.data",
    min.cells = NULL,
    min.genes = 200,
    not.expressed = 0,
    verbose = TRUE
) {

    # Check if the object is valid for filtering
    CheckInput <- function(object, slot) {
        stopifnot(class(object)[1] == "Signac")
        stopifnot(class(attr(object, slot))[1] == "dgCMatrix")
    }

    Cat <- function(meow, ...) {
        if (meow) {
            cat(..., '\n')
        }
    }

    FilterGenes <- function(data, min.cells, not.expressed) {
        if (is.null(min.cells)) {
            min.cells <- max(10, ncol(data) * 0.02)
        }
        Cat(verbose, "[Signac] min.cells:", min.cells)
        data <- data[Matrix::rowSums(data > not.expressed) >= min.cells, ]
        return(data)
    }

    FilterCells <- function(data, min.genes, not.expressed) {
        Cat(verbose, "[Signac] min.genes:", min.genes)
        return(data <- data[, Matrix::colSums(data > not.expressed) >= min.genes])
    }

    CheckInput(object, slot)
    data <- attr(object, slot)
    data <- FilterGenes(data, min.cells, not.expressed)
    data <- FilterCells(data, min.genes, not.expressed)
    Cat(verbose, "[Signac] After filtering:", nrow(data), "genes X", ncol(data), "cells")
    object@filtered.data <- data
    return(object)
}