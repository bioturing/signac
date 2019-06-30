context("test-io")

test_that("WriteSpMtAsS4", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "bioturing"
    auto.update <- FALSE
    mat <- ReadSpMt( h5.path, group.name)
    status <- WriteSpMtAsS4(tempfile(), "matrixv1", mat)
    expect_equal(status, NULL)
})

test_that("ReadSpMtAsS4", {
    h5.path <- system.file("extdata", "v3inputb.h5", package = "Signac")
    group.name <- "bioturing"
    mat <- ReadSpMtAsS4(h5.path, group.name)
    expect_equal(TRUE, TRUE)
})
