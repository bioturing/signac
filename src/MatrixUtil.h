#ifndef MATRIX_UTIL
#define MATRIX_UTIL

#include <RcppArmadillo.h>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;

arma::mat FastMatMult(const arma::mat &mat1, const arma::mat &mat2);
arma::mat FastGetRowsOfMat(const arma::mat &mat, arma::uvec vec);
arma::mat FastGetColsOfMat(const arma::mat &mat, arma::uvec vec);
arma::mat FastGetSubMat(const arma::mat &mat, arma::uvec rvec, arma::uvec cvec);

#endif //MATRIX_UTIL
