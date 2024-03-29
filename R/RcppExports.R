# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' FastGetCurrentDate
#'
#' This function returns a current date (YYYY-MM-DD)
#'
#' @export
FastGetCurrentDate <- function() {
    .Call(`_Signac_FastGetCurrentDate`)
}

#' FastDiffVector
#'
#' This function is used to diff 2 vector
#'
#' @param a An integer vector
#' @param b An integer vector
#' @export
FastDiffVector <- function(a, b) {
    .Call(`_Signac_FastDiffVector`, a, b)
}

#' FastRandVector
#'
#' This function create a random vector
#'
#' @param num An integer number
#' @export
FastRandVector <- function(num) {
    .Call(`_Signac_FastRandVector`, num)
}

#' WriteSpMtAsS4
#'
#' This function is used to write a sparse S4 matrix
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @param mat A sparse matrix
NULL

#' WriteSpMtAsSpMat
#'
#' This function is used to write a sparse ARMA matrix
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @param mat A sparse matrix
#' @export
WriteSpMtAsSpMat <- function(filePath, groupName, mat) {
    invisible(.Call(`_Signac_WriteSpMtAsSpMat`, filePath, groupName, mat))
}

#' WriteSpMtAsSpMatFromS4
#'
#' This function is used to write a sparse ARMA matrix from S4
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @param mat A sparse matrix
#' @export
WriteSpMtAsSpMatFromS4 <- function(filePath, groupName, mat) {
    invisible(.Call(`_Signac_WriteSpMtAsSpMatFromS4`, filePath, groupName, mat))
}

#' ReadSpMtAsSPMat
#'
#' This function is used to read a sparse matrix from HDF5 file
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @export
ReadSpMtAsSPMat <- function(filePath, groupName) {
    .Call(`_Signac_ReadSpMtAsSPMat`, filePath, groupName)
}

#' ReadSpMtAsS4
#'
#' This function is used to read a sparse matrix from HDF5 file
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @export
ReadSpMtAsS4 <- function(filePath, groupName) {
    .Call(`_Signac_ReadSpMtAsS4`, filePath, groupName)
}

#' ReadRowSumSpMt
#'
#' Read rows sums
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @export
ReadRowSumSpMt <- function(filePath, groupName) {
    .Call(`_Signac_ReadRowSumSpMt`, filePath, groupName)
}

#' ReadColSumSpMt
#'
#' Read cols sums
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @export
ReadColSumSpMt <- function(filePath, groupName) {
    .Call(`_Signac_ReadColSumSpMt`, filePath, groupName)
}

#' GetListAttributes
#'
#' Get list attribute of a dataset
#'
#' @param filePath A HDF5 path
#' @param groupName A string (HDF5 dataset)
#' @param datasetName A dataset name
#' @export
GetListAttributes <- function(filePath, groupName, datasetName) {
    .Call(`_Signac_GetListAttributes`, filePath, groupName, datasetName)
}

#' GetListObjectNames
#'
#' Get list object of a group
#'
#' @param filePath A HDF5 path
#' @param groupName A string (HDF5 dataset)
#' @export
GetListObjectNames <- function(filePath, groupName) {
    .Call(`_Signac_GetListObjectNames`, filePath, groupName)
}

#' GetListRootObjectNames
#'
#' Get list groups
#'
#' @param filePath An HDF5 path
#' @export
GetListRootObjectNames <- function(filePath) {
    .Call(`_Signac_GetListRootObjectNames`, filePath)
}

#' Read10XH5
#'
#' Get list triplet of SEXP
#'
#' @param filePath A fiel path
#' @param use_names Use names flag
#' @param unique_features Unique features flag
#' @export
Read10XH5Content <- function(filePath, use_names, unique_features) {
    .Call(`_Signac_Read10XH5Content`, filePath, use_names, unique_features)
}

#' WriteRootDataset
#'
#' Write a string vector to root group
#'
#' @param filePath A fiel path
#' @param datasetName A dataset name
#' @param datasetVal A string vector
#' @export
WriteRootDataset <- function(filePath, datasetName, datasetVal) {
    invisible(.Call(`_Signac_WriteRootDataset`, filePath, datasetName, datasetVal))
}

#' ReadIntegerVector
#'
#' This function is used to read a integer vector from hdf5 file
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @export
ReadIntegerVector <- function(filePath, groupName, datasetName) {
    .Call(`_Signac_ReadIntegerVector`, filePath, groupName, datasetName)
}

#' ReadDoubleVector
#'
#' This function is used to read a double vector from hdf5 file
#'
#' @param filePath A string (HDF5 path)
#' @param groupName A string (HDF5 dataset)
#' @export
ReadDoubleVector <- function(filePath, groupName, datasetName) {
    .Call(`_Signac_ReadDoubleVector`, filePath, groupName, datasetName)
}

#' FastMatMult
#'
#' This function is used to add two matrix
#'
#' @param mat1 A matrix
#' @param mat2 A matrix
#' @export
FastMatMult <- function(mat1, mat2) {
    .Call(`_Signac_FastMatMult`, mat1, mat2)
}

#' FastGetRowsOfMat
#'
#' This function is used to get some rows of matrix
#'
#' @param mat A matrix
#' @param vec A row vector
#' @export
FastGetRowsOfMat <- function(mat, vec) {
    .Call(`_Signac_FastGetRowsOfMat`, mat, vec)
}

FastGetColsOfMat <- function(mat, vec) {
    .Call(`_Signac_FastGetColsOfMat`, mat, vec)
}

#' FastGetSubMat
#'
#' This function is used to get subview matrix
#'
#' @param mat A matrix
#' @param vec A row vector
#' @param vec A col vector
#' @export
FastGetSubMat <- function(mat, rvec, cvec) {
    .Call(`_Signac_FastGetSubMat`, mat, rvec, cvec)
}

#' FastCreateSparseMat
#'
#' Generate a sparse matrix with the elements along the main diagonal set to one and off-diagonal elements set to zero
#'
#' @param nrow An integer number
#' @param ncol An integer number
#' @export
FastCreateSparseMat <- function(nrow, ncol) {
    .Call(`_Signac_FastCreateSparseMat`, nrow, ncol)
}

#' FastStatsOfSparseMat
#'
#' Get statistic of sparse matrix
#'
#' @param mat An sparse matrix
#' @export
FastStatsOfSparseMat <- function(mat) {
    .Call(`_Signac_FastStatsOfSparseMat`, mat)
}

#' FastCreateFromTriplet
#'
#' Create sparse matrix from triplet
#'
#' @param vec1 An integer vector
#' @param vec2 An integer vector
#' @param vec_val An double vector
#' @export
FastCreateFromTriplet <- function(vec1, vec2, vec_val) {
    .Call(`_Signac_FastCreateFromTriplet`, vec1, vec2, vec_val)
}

#' FastConvertToSparseMat
#'
#' Convert SEXP to sparse matrix
#'
#' @param s A SEXP type
#' @export
FastConvertToSparseMat <- function(s) {
    .Call(`_Signac_FastConvertToSparseMat`, s)
}

#' FastConvertToTripletMat
#'
#' Get list triplet of SEXP
#'
#' @param s A SEXP type
#' @export
FastConvertToTripletMat <- function(s) {
    .Call(`_Signac_FastConvertToTripletMat`, s)
}

#' FastSparseMatSqrt
#'
#' SQRT element i,j in a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastSparseMatSqrt <- function(mat) {
    .Call(`_Signac_FastSparseMatSqrt`, mat)
}

#' FastSparseMatMult
#'
#' Multiply two sparse matrix
#'
#' @param mat1 A sparse matrix
#' @param mat2 A sparse matrix
#' @export
FastSparseMatMult <- function(mat1, mat2) {
    .Call(`_Signac_FastSparseMatMult`, mat1, mat2)
}

#' FastSparseMatAddition
#'
#' Add two sparse matrix
#'
#' @param mat1 A sparse matrix
#' @param mat2 A sparse matrix
#' @export
FastSparseMatAddition <- function(mat1, mat2) {
    .Call(`_Signac_FastSparseMatAddition`, mat1, mat2)
}

#' FastSparseMatMultWithNum
#'
#' Multiply sparse matrix with number
#'
#' @param mat A sparse matrix
#' @param k A integer number
#' @export
FastSparseMatMultWithNum <- function(mat, k) {
    .Call(`_Signac_FastSparseMatMultWithNum`, mat, k)
}

#' FastSparseMatSymmatl
#'
#' Symmatl sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastSparseMatSymmatl <- function(mat) {
    .Call(`_Signac_FastSparseMatSymmatl`, mat)
}

#' FastSparseMatTranspose
#'
#' Transpose sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastSparseMatTranspose <- function(mat) {
    .Call(`_Signac_FastSparseMatTranspose`, mat)
}

#' FastSparseMatTrimatu
#'
#' Trimatu sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastSparseMatTrimatu <- function(mat) {
    .Call(`_Signac_FastSparseMatTrimatu`, mat)
}

#' FastSparseMatTrace
#'
#' Trace sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastSparseMatTrace <- function(mat) {
    .Call(`_Signac_FastSparseMatTrace`, mat)
}

#' FastConvertToDiagonalSparseMat
#'
#' Get diag sparse matrix from a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastConvertToDiagonalSparseMat <- function(mat) {
    .Call(`_Signac_FastConvertToDiagonalSparseMat`, mat)
}

#' FastSparseMatSquare
#'
#' Square a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastSparseMatSquare <- function(mat) {
    .Call(`_Signac_FastSparseMatSquare`, mat)
}

#' FastSparseMatRepmat
#'
#' Repmat a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastSparseMatRepmat <- function(mat, i, j) {
    .Call(`_Signac_FastSparseMatRepmat`, mat, i, j)
}

#' FastSparseMatSign
#'
#' Sign a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastSparseMatSign <- function(mat) {
    .Call(`_Signac_FastSparseMatSign`, mat)
}

#' FastSparseMatMultSD
#'
#' Multiply a sparse matrix with a dense matrix
#'
#' @param mat1 A sparse matrix
#' @param mat2 A dense matrix
#' @export
FastSparseMatMultSD <- function(mat1, mat2) {
    .Call(`_Signac_FastSparseMatMultSD`, mat1, mat2)
}

#' FastSparseMatMultDD
#'
#' Multiply two dense matrix
#'
#' @param mat1 A dense matrix
#' @param mat2 A dense matrix
#' @export
FastSparseMatMultDD <- function(mat1, mat2) {
    .Call(`_Signac_FastSparseMatMultDD`, mat1, mat2)
}

#' FastGetRowOfSparseMat
#'
#' Get row of sparse matrix
#'
#' @param mat A dense matrix
#' @param i A integer number
#' @export
FastGetRowOfSparseMat <- function(mat, i) {
    .Call(`_Signac_FastGetRowOfSparseMat`, mat, i)
}

#' FastGetColOfSparseMat
#'
#' Get col of sparse matrix
#'
#' @param mat A dense matrix
#' @param i A integer number
#' @export
FastGetColOfSparseMat <- function(mat, j) {
    .Call(`_Signac_FastGetColOfSparseMat`, mat, j)
}

#' FastGetRowsOfSparseMat
#'
#' Get rows of sparse matrix
#'
#' @param mat A sparse matrix
#' @param start A integer number
#' @param end A integer number
#' @export
FastGetRowsOfSparseMat <- function(mat, start, end) {
    .Call(`_Signac_FastGetRowsOfSparseMat`, mat, start, end)
}

#' FastGetColsOfSparseMat
#'
#' Get cols of sparse matrix
#'
#' @param mat A sparse matrix
#' @param start A integer number
#' @param end A integer number
#' @export
FastGetColsOfSparseMat <- function(mat, start, end) {
    .Call(`_Signac_FastGetColsOfSparseMat`, mat, start, end)
}

#' FastGetSubSparseMat
#'
#' Subview sparse matrix
#'
#' @param mat A sparse matrix
#' @param rrvec A row vector
#' @param ccvec A col vector
#' @param need_perform_row A bool
#' @param need_perform_col A bool
#' @export
FastGetSubSparseMat <- function(mat, rrvec, ccvec, need_perform_row, need_perform_col) {
    .Call(`_Signac_FastGetSubSparseMat`, mat, rrvec, ccvec, need_perform_row, need_perform_col)
}

#' FastGetSubSparseMatByRows
#'
#' Get subview of a sparse matrix by rows
#'
#' @param mat A sparse matrix
#' @param rvec A row vector
#' @export
FastGetSubSparseMatByRows <- function(mat, rvec) {
    .Call(`_Signac_FastGetSubSparseMatByRows`, mat, rvec)
}

#' FastGetSubSparseMatByCols
#'
#' Get subview of a sparse matrix by cols
#'
#' @param mat A sparse matrix
#' @param cvec A col vector
#' @export
FastGetSubSparseMatByCols <- function(mat, cvec) {
    .Call(`_Signac_FastGetSubSparseMatByCols`, mat, cvec)
}

#' FastGetSumSparseMatByRows
#'
#' Sum of rows in a sparse matrix
#'
#' @param mat A sparse matrix
#' @param rvec A col vector
#' @export
FastGetSumSparseMatByRows <- function(mat, rvec) {
    .Call(`_Signac_FastGetSumSparseMatByRows`, mat, rvec)
}

#' FastGetSumSparseMatByCols
#'
#' Sum of cols in a sparse matrix
#'
#' @param mat A sparse matrix
#' @param cvec A col vector
#' @export
FastGetSumSparseMatByCols <- function(mat, cvec) {
    .Call(`_Signac_FastGetSumSparseMatByCols`, mat, cvec)
}

#' FastGetSumSparseMatByAllRows
#'
#' Sum all rows in a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastGetSumSparseMatByAllRows <- function(mat) {
    .Call(`_Signac_FastGetSumSparseMatByAllRows`, mat)
}

#' FastGetSumSparseMatByAllCols
#'
#' Sum all cols in a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastGetSumSparseMatByAllCols <- function(mat) {
    .Call(`_Signac_FastGetSumSparseMatByAllCols`, mat)
}

#' FastGetMedianSparseMatByAllRows
#'
#' Median all rows in a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastGetMedianSparseMatByAllRows <- function(mat) {
    .Call(`_Signac_FastGetMedianSparseMatByAllRows`, mat)
}

#' FastGetMedianSparseMatByAllCols
#'
#' Median all cols in a sparse matrix
#'
#' @param mat A sparse matrix
#' @export
FastGetMedianSparseMatByAllCols <- function(mat) {
    .Call(`_Signac_FastGetMedianSparseMatByAllCols`, mat)
}

#' VeniceMarker
#'
#' Find gene marker for a cluster in sparse matrix
#'
#' @param S4_mtx A sparse matrix
#' @param cluster A numeric vector
#' @export
VeniceMarker <- function(S4_mtx, cluster, threshold = 0L, threshold_pct = 0, perm = 0L, correct = TRUE, verbose = FALSE) {
    .Call(`_Signac_VeniceMarker`, S4_mtx, cluster, threshold, threshold_pct, perm, correct, verbose)
}

#' VeniceMarkerTransposedH5
#'
#' Find gene marker for a cluster in tranposed H5 file
#'
#' @param hdf5Path A string path
#' @param cluster A numeric vector
#' @export
VeniceMarkerTransposedH5 <- function(hdf5Path, group_name, cluster, threshold = 0L, threshold_pct = 0, perm = 0L, correct = TRUE, verbose = FALSE) {
    .Call(`_Signac_VeniceMarkerTransposedH5`, hdf5Path, group_name, cluster, threshold, threshold_pct, perm, correct, verbose)
}

