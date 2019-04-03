context("test-io")

test_that("Read from zip tsv", {
    file.test <- system.file("extdata", "GSM2629435_AB2430.txt.gz", package = "Signac")
    obj <- CreateSignacObject(file.test, type = "tsv")
    expect_equal("Signac" %in%  class(obj), TRUE)
})

test_that("Read from mtx", {
    file.test <- system.file("extdata", "10xLite", package = "Signac")
    obj <- CreateSignacObject(file.test, type = "mtx")
    expect_equal("Signac" %in%  class(obj), TRUE)
})

test_that("Read from H5", {
    file.test <- system.file("extdata", "v3input.h5", package = "Signac")
    obj <- CreateSignacObject(file.test, type = "h5")
    expect_equal("Signac" %in%  class(obj), TRUE)
})

test_that("ReadSpMt", {
    h5.path <- system.file("extdata", "v3input.h5", package = "Signac")
    group.name = "matrix"
    auto.update = TRUE
    mat <- Signac::ReadSpMt( h5.path, group.name, auto.update)
    expect_equal("dgCMatrix" %in%  class(mat), TRUE)
})
