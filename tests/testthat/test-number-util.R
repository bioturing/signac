library(Matrix)
library(numbers)
context("test-number-util")

test_that("Compute GCD", {
    a <- 3214
    b <- 9876
    x <- FastComputeGCD(a,b)
    expect_equal(x, GCD(a,b))
})

test_that("Compute LCM", {
    a <- 3214
    b <- 9876
    x <- FastComputeLCM(a,b)
    expect_equal(x, LCM(a,b))
})
