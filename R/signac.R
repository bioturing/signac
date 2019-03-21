#' Signac class
#'
#' All functions will operate on this object
#' @slot raw.data The original expression matrix
#' @slot corrected.data Data after batch effect removal
#' @slot scaled.data Data after normalization
#' @slot log.data Data after log
#' @slot cell.names Id of cells
#' @slot batch Batches for data come from multiple runs/sources
#' @slot metadata Labels and annotations
#' @slot vdj TCR information
#' @slot ssn Shared nearest neighbor
#' @slot dimred A named list of dimensionality reduction result
#' @slot assay A named list of multimodal analysis
#' @slot trajectory Pseudotime result
signac <- setClass(
    "signac",

    # Data
    raw.data = "ANY",
    corrected.data = "ANY",
    scaled.data = "ANY",
    log.data = "ANY",

    # Metadata
    cell.names = "character"
    batch = "character",
    metadata = "data.frame",
    vdj = "data.frame"

    # Supported data
    snn = "list",
    dimred = "list",
    assay = "list",
    trajectory = "list"
)
