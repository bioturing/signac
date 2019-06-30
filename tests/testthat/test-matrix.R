library(Matrix)
context("test-matrix-util")

test_that("FastMatMult", {
    m1  <- matrix(1:9, 3, 3)
    m2  <- matrix(1:12, 3, 4)
    FinalMAT <- FastMatMult(m1, m2)
    expect_equal(length(FinalMAT), length(m2))
})

test_that("FastGetRowsOfMat", {
    m1  <- matrix(1:9, 3, 3)
    FinalMAT <- FastGetRowsOfMat(m1, c(1,3))
    expect_equal(length(FinalMAT), 6)
})

test_that("FastGetColsOfMat", {
    m1  <- matrix(1:9, 3, 3)
    FinalMAT <- FastGetColsOfMat(m1, c(1,3))
    expect_equal(length(FinalMAT), 6)
})

test_that("FastGetSubMat", {
    m1  <- matrix(1:9, 3, 3)
    FinalMAT <- FastGetSubMat(m1, c(1,3), c(2,3))
    expect_equal(length(FinalMAT), 4)
})
