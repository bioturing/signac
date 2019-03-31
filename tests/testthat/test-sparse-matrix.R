library(Matrix)
context("test-sparse-matrix-util")

test_that("FastConvertToSparseMat", {
    mtxt <- c("11   0   0  14   0  16",
              " 0  22   0   0  25  26",
              " 0   0  33  34   0  36",
              "41   0  43  44   0  46")
    M <- as.matrix(read.table(text=mtxt))
    dimnames(M) <- NULL
    SM <- Matrix(M, sparse=TRUE)
    MAT <- FastConvertToSparseMat(SM)
    expect_equal(length(SM), length(MAT))
})

test_that("FastConvertToTripletMat", {
    i <- c(1,3:8)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    SM <- sparseMatrix(i, j, x = x)
    List <- FastConvertToTripletMat(SM)
    expect_equal(List$nrow, 8)
})

test_that("FastCreateSparseMat", {
    i <- c(1,3:8)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    SM <- sparseMatrix(i, j, x = x)
    MAT <- FastCreateSparseMat(length(i), length(j))
    expect_equal(length(MAT), 49)
})

test_that("FastCreateFromTriplet", {
    i <- c(1,4:9)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    SM <- sparseMatrix(i, j, x = x)
    MAT <- Signac::FastCreateFromTriplet(i, j, x)
    expect_equal(length(MAT), 110)
})

test_that("FastStatsOfSparseMat", {
    i <- c(1,4:9)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    SM <- sparseMatrix(i, j, x = x)
    List <- Signac::FastStatsOfSparseMat(SM)
    expect_equal(List[[1]], 9)
})

test_that("FastSparseMatTranspose", {
    i <- c(1,4:9)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    MAT <- sparseMatrix(i, j, x = x)
    TMAT <- Signac::FastSparseMatTranspose(MAT)
    expect_equal(length(TMAT), length(MAT))
})

test_that("FastSparseMatSqrt", {
    i <- c(1,4:9)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    MAT <- sparseMatrix(i, j, x = x)
    sqrtMAT <- Signac::FastSparseMatSqrt(MAT)
    expect_equal(length(sqrtMAT), length(MAT))
})

test_that("FastSparseMatMult", {
    set.seed(123)
    v1 <- sample(1e3)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.99, .01)),
                 length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)

    v2 <- sample(1e3)
    m2  <- Matrix(sample(c(0, 1), length(v2) ^ 2, T, c(.89, .01)),
                  length(v2), length(v2), sparse = F)
    MAT2 <- Matrix(m2, sparse = T)

    FinalMAT <- Signac::FastSparseMatMult(MAT1, MAT2)
    expect_equal(length(FinalMAT), length(MAT1))
})

test_that("FastSparseMatSymmatl", {
    MAT <- new("dtCMatrix", Dim = c(5L, 5L), uplo = "L",
               x = c(10, 1, 3, 10, 1, 10, 1, 10, 10),
               i = c(0L, 2L, 4L, 1L, 3L,2L, 4L, 3L, 4L),
               p = c(0L, 3L, 5L, 7:9))
    SMAT <- Signac::FastSparseMatSymmatl(MAT)
    expect_equal(length(SMAT), length(MAT))
})

test_that("FastSparseMatAddition", {
    set.seed(123)
    v1 <- sample(1e3)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.99, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)

    v2 <- sample(1e3)
    m2  <- Matrix(sample(c(0, 1), length(v2) ^ 2, T, c(.89, .01)),
                  length(v2), length(v2), sparse = F)
    MAT2 <- Matrix(m2, sparse = T)

    FinalMAT <- Signac::FastSparseMatAddition(MAT1, MAT2)
    expect_equal(length(FinalMAT), length(MAT1))
})


test_that("FastSparseMatMultWithNum", {
    set.seed(123)
    v1 <- sample(1e1)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.99, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)

    FinalMAT <- Signac::FastSparseMatMultWithNum(MAT1, 2)
    expect_equal(length(FinalMAT), length(MAT1))
})

test_that("FastSparseMatTrimatu", {
    i <- c(1,4:9)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    MAT <- sparseMatrix(i, j, x = x, symmetric = TRUE)
    TMAT <- Signac::FastSparseMatTrimatu(MAT)
    expect_equal(length(TMAT), length(MAT))
})

test_that("FastSparseMatTrace - sum of diagonal elements", {
    dgT <- new("dgTMatrix",
               i = c(1L,1L,0L,3L,3L),
               j = c(2L,2L,4L,0L,0L),
               x=10*1:5, Dim=4:5)
    dgT_t <- Signac::FastSparseMatTranspose(dgT)
    MAT <- Signac::FastSparseMatMult(dgT, dgT_t)
    iNum <- Signac::FastSparseMatTrace(MAT)
    expect_equal(iNum, 9900)
})

test_that("FastConvertToDiagonalSparseMat", {
    i <- c(1,4:9)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    MAT <- sparseMatrix(i, j, x = x)
    TMAT <- Signac::FastConvertToDiagonalSparseMat(MAT)
    expect_equal(length(TMAT), length(MAT))
})

test_that("FastSparseMatSquare", {
    i <- c(1,4:9)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    MAT <- sparseMatrix(i, j, x = x)
    TMAT <- Signac::FastSparseMatSquare(MAT)
    expect_equal(length(TMAT), length(MAT))
})

test_that("FastSparseMatRepmat", {
    MAT <- new("dtRMatrix", Dim = c(2L,2L),
               x = c(5, 1:2), p = c(0L,2:3), j= c(0:1,1L))
    TMAT <- Signac::FastSparseMatRepmat(MAT, 2, 2)
    expect_equal(length(TMAT), 4 * length(MAT))
})

test_that("FastSparseMatSign", {
    MAT <- new("dsRMatrix", Dim = c(2L,2L),
               x = c(-3,1), j = c(1L,1L), p = 0:2)
    TMAT <- Signac::FastSparseMatSign(MAT)
    expect_equal(length(TMAT), length(MAT))
})

test_that("FastSparseMatMultSD", {
    set.seed(123)
    v1 <- sample(5e1)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.99, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)

    n <- 5e1
    m2 <- rsparsematrix(n, n, 0.11, rand.x=function(n) rpois(n, 1) + 1)
    MAT2 <- as.matrix(m2)

    FinalMAT <- Signac::FastSparseMatMultSD(MAT1, MAT2)
    expect_equal(length(FinalMAT), 2500)
})

test_that("FastSparseMatMultDS", {
    set.seed(123)
    n <- 5e3
    m1 <- rsparsematrix(n, n, 0.11, rand.x=function(n) rpois(n, 1) + 1)
    MAT1 <- as.matrix(m1)

    v2 <- sample(5e3)
    m2  <- Matrix(sample(c(0, 1), length(v2) ^ 2, T, c(.89, .01)),
                  length(v2), length(v2), sparse = F)
    MAT2 <- Matrix(m2, sparse = T)

    FinalMAT <- Signac::FastSparseMatMultDS(MAT1, MAT2)
    expect_equal(length(FinalMAT), 25000000)
})

test_that("FastSparseMatMultDD", {
    set.seed(123)
    n <- 5e3
    a <- rsparsematrix(n, n, 0.01, rand.x=function(n) rpois(n, 1) + 1)
    b <- rsparsematrix(n, n, 0.01, rand.x=function(n) rpois(n, 1) + 1)

    #Dense matrix
    MAT1 <- as.matrix(a)
    MAT2 <- as.matrix(b)

    FinalMAT <- Signac::FastSparseMatMultDD(MAT1, MAT2)
    expect_equal(length(FinalMAT), length(MAT1))
})

test_that("FastGetRowOfSparseMat", {
    MAT <- new("dsRMatrix", Dim = c(2L,2L),
               x = c(-3,1), j = c(1L,1L), p = 0:2)
    row <- 1
    TMAT <- Signac::FastGetRowOfSparseMat(MAT, row)
    expect_equal(length(TMAT), 2)
})

test_that("FastGetColOfSparseMat", {
    MAT <- new("dsRMatrix", Dim = c(2L,2L),
               x = c(-3,1), j = c(1L,1L), p = 0:2)
    col <- 2
    TMAT <- Signac::FastGetColOfSparseMat(MAT, col)
    expect_equal(length(TMAT), col)
})

test_that("FastGetSubSparseMat", {
    set.seed(123)
    v1 <- sample(1e2)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.89, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)
    FinalMAT <- Signac::FastGetSubSparseMat(MAT1, c(3,4), c(2,5,7))
    expect_equal(length(FinalMAT), 6)
})

test_that("FastGetSubSparseMatByRows", {
    set.seed(123)
    v1 <- sample(1e1)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.89, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)
    FinalMAT <- Signac::FastGetSubSparseMatByRows(MAT1, c(4992,4997))
    expect_equal(length(FinalMAT), 20)
})

test_that("FastGetSubSparseMatByCols", {
    set.seed(123)
    v1 <- sample(1e1)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.89, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)
    FinalMAT <- Signac::FastGetSubSparseMatByCols(MAT1, c(5,16))
    expect_equal(length(FinalMAT), 20)
})

test_that("FastGetSumSparseMatByRows", {
    set.seed(123)
    v1 <- sample(1e1)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.89, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)
    ListSum <- Signac::FastGetSumSparseMatByRows(MAT1, c(4,16))
    expect_equal(length(ListSum), 2)
})

test_that("FastGetSumSparseMatByAllRows", {
    set.seed(123)
    v1 <- sample(1e3)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.89, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)
    ListSum <- Signac::FastGetSumSparseMatByAllRows(MAT1)
    expect_equal(length(FinalMAT), 20)
})

test_that("FastGetSumSparseMatByAllRowsV2", {
    set.seed(123)
    v1 <- sample(1e3)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.89, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)
    ListSum <- Signac::FastGetSumSparseMatByAllRowsV2(MAT1)
    expect_equal(length(FinalMAT), 20)
})

test_that("FastGetSumSparseMatByAllCols", {
    set.seed(123)
    v1 <- sample(1e3)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.89, .01)),
                  length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)
    ListSum <- Signac::FastGetSumSparseMatByAllCols(MAT1)
    expect_equal(length(FinalMAT), 20)
})
