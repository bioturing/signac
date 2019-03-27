library(Matrix)
context("test-sparse-matrix-util")

test_that("Sqrt sparse matrix", {
    i <- c(1,4:9)
    j <- c(2,9,6:10)
    x <- 7 * (1:7)
    MAT <- sparseMatrix(i, j, x = x)
    sqrtMAT <- Signac::FastSparseMatSqrt(MAT)
    expect_equal(length(sqrtMAT), length(MAT))
})

test_that("Multiply two sparse matrix", {
    T0 = Sys.time()
    set.seed(123)
    v1 <- sample(1e3)
    m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.99, .01)),
                 length(v1), length(v1), sparse = F)
    MAT1 <- Matrix(m1, sparse = T)
    #T1 = Sys.time()
    #cat("\nTime to init MAT1:", T1 - T0)

    v2 <- sample(1e3)
    m2  <- Matrix(sample(c(0, 1), length(v2) ^ 2, T, c(.99, .01)),
                  length(v2), length(v2), sparse = F)
    MAT2 <- Matrix(m2, sparse = T)
    #T2 = Sys.time()
    #cat("\nTime to init MAT2:", T2 - T1)

    FinalMAT <- Signac::FastSparseMatMult(MAT1, MAT2)
    #T3 = Sys.time()
    #cat("\nTime3:", T3 - T2)
    expect_equal(length(FinalMAT), length(MAT1))
})
