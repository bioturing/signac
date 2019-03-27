library(rbenchmark)
library(numbers)
library(Signac)

a <- 3214
b <- 9876

res <- benchmark(r1 = c(FastComputeGCD(a,b), FastComputeLCM(a,b)),
                 r2 = c(GCD(a,b), LCM(a,b)),
                 replications = 50000)
print(res[,1:4])
