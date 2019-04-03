context("test-io")

test_that("WriteSpMtV2", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "matrix"
    auto.update <- FALSE
    mat <- Signac::ReadSpMt( h5.path, group.name, auto.update)
    status <- Signac::WriteSpMtV2(h5.path, "matrixv2", mat)
    expect_equal(status, TRUE)
})

test_that("WriteSpMtV1", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "matrix"
    auto.update <- FALSE
    mat <- Signac::ReadSpMt( h5.path, group.name, auto.update)
    status <- Signac::WriteSpMtV1(h5.path, "matrixv1", mat)
    expect_equal(status, TRUE)
})

test_that("ReadSpMtV2", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "matrixv2"
    auto.update <- FALSE
    mat <- Signac::ReadSpMtV2(h5.path, group.name)
    expect_equal(TRUE, TRUE)
})

test_that("ReadSpMtV1", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "matrix"
    mat <- Signac::ReadSpMtV1( h5.path, group.name)
    expect_equal(status, TRUE)
})
