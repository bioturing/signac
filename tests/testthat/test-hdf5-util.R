context("test-io")

test_that("WriteSpMtAsSPMat", {
    h5.path <- system.file("extdata", "v3input.h5", package = "Signac")
    group.name <- "matrix"
    auto.update <- FALSE
    mat <- Signac::ReadSpMt( h5.path, group.name, auto.update)
    status <- Signac::WriteSpMtAsSPMat(h5.path, "matrixv2", mat)
    expect_equal(status, TRUE)
})

test_that("WriteSpMtAsS4", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "matrix"
    auto.update <- FALSE
    mat <- Signac::ReadSpMt( h5.path, group.name, auto.update)
    status <- Signac::WriteSpMtAsS4(h5.path, "matrixv1", mat)
    expect_equal(status, TRUE)
})

test_that("ReadSpMtAsSPMat", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "matrixv2"
    auto.update <- FALSE
    mat <- Signac::ReadSpMtAsSPMat(h5.path, group.name)
    expect_equal(TRUE, TRUE)
})

test_that("ReadSpMtAsS4", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "matrix"
    mat <- Signac::ReadSpMtAsS4( h5.path, group.name)
    expect_equal(status, TRUE)
})
