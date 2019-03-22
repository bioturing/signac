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
    "Signac",
    slots = c(
        # Data
        raw.data = "ANY",
        corrected.data = "ANY",
        scaled.data = "ANY",
        log.data = "ANY",

        # Metadata
        cell.names = "character",
        batch = "character",
        metadata = "data.frame",
        vdj = "data.frame",

        # Supported data
        snn = "list",
        dimred = "list",
        assay = "list",
        trajectory = "ANY",

        # Version
        version = "ANY",

        # Project name
        project.name = "ANY"
    )
)

# How to show signac
setMethod("show", "Signac", function(object) {
    cat(
        # Intro
        "Class:", Colourise(class(object), "green"), "\n",

        # Detail
        "Original data:", nrow(object@raw.data), "genes X",
        ncol(object@raw.data), "cells\n",

        # Project name
        "Project name:", Colourise(object@project.name, "cyan"), "\n",

        # Signac version
        "Signac version:", Colourise(object@version, "cyan"), "\n"
    )
})
