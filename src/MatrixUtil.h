#ifndef MATRIX_UTIL
#define MATRIX_UTIL

#include <RcppArmadillo.h>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;

arma::mat FastMatMult(const arma::mat &mat1, const arma::mat &mat2);

#endif //MATRIX_UTIL
