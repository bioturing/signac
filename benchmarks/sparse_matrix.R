library(rbenchmark)
library(numbers)
library(Signac)
library(Matrix)
library(microbenchmark)

set.seed(123)
n <- 5e3
a <- rsparsematrix(n, n, 0.01, rand.x=function(n) rpois(n, 1) + 1)
m2 <- as.matrix(a)
MAT <- Matrix(m2, sparse = T)

stopifnot(identical(Signac::FastGetSumSparseMatByAllRows(MAT), Matrix::rowSums(MAT)))
stopifnot(identical(Signac::FastGetSumSparseMatByAllCols(MAT), Matrix::colSums(MAT)))

res <- microbenchmark::microbenchmark(rFast = Signac::FastGetSumSparseMatByAllRows(MAT),
                 rBultin = Matrix::rowSums(MAT),
                 cFast = Signac::FastGetSumSparseMatByAllCols(MAT),
                 cBultin = Matrix::colSums(MAT), times = 500)

#res <- microbenchmark(Matrix::rowSums(MAT),
#                 Signac::FastGetSumSparseMatByAllRows(MAT),
#                 order="relative")

print(res)
