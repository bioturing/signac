context("test-preprocess")

test_that("Filter basic", {
    file.test <- system.file("extdata", "GSM2629435_AB2430.txt.gz", package = "Signac")
    obj <- CreateSignacObject(file.test, type = "tsv")
    obj <- FilterDataBasic(obj, verbose = FALSE)
    expect_equal(class(obj@filtered.data)[1], "dgCMatrix")
    expect_equal(nrow(obj@filtered.data), 2298)
    expect_equal(ncol(obj@filtered.data), 163)
})