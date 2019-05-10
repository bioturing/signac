library(Matrix)
library(numbers)
context("test-number-util")

test_that("FastGetCurrentDate", {
    x <- Signac::FastGetCurrentDate()
    expect_equal(length(x), 1)
})

test_that("FastDiffVector", {
    a <- c(1,3,4,7)
    b <- c(3,4,8)
    x <- Signac::FastDiffVector(a,b)
    expect_equal(length(x), 2)
})
