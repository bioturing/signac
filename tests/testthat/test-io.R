context("test-io")
library(Signac)

test_that("Read from zip tsv works", {
    file_test <- system.file("extdata", "GSM2629435_AB2430.txt.gz", package = "Signac")
    obj <- CreateSignacObject(file_test, type = "tsv")
    expect_equal("Signac" %in%  class(obj), TRUE)
})
