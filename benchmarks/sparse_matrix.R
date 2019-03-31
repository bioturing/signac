library(rbenchmark)
library(numbers)
library(Signac)
library(Matrix)

set.seed(123)
v1 <- sample(1e3)
m1  <- Matrix(sample(c(0, 1), length(v1) ^ 2, T, c(.89, .01)),
              length(v1), length(v1), sparse = F)
MAT <- Matrix(m1, sparse = T)

res <- benchmark(r1 = c(Signac::FastGetSumSparseMatByAllRows(MAT), Signac::FastGetSumSparseMatByAllCols(MAT)),
                 r2 = c(Matrix::rowSums(MAT), Matrix::colSums(MAT)),
                 replications = 500)
print(res[,1:4])
