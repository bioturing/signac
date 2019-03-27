library(Matrix)
context("test-matrix-util")

test_that("Multiply two matrix", {
    m1  <- matrix(1:9, 3, 3)
    m2  <- matrix(1:12, 3, 4)
    FinalMAT <- Signac::FastMatMult(m1, m2)
    expect_equal(length(FinalMAT), length(m2))
})
