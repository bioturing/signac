#ifndef COMMON_MATRIX_UTIL
#define COMMON_MATRIX_UTIL

#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG
#define ARMA_USE_HDF5

#include <RcppArmadillo.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;
using namespace arma;

int FastComputeGCD(const int &num1, const int &num2);
int FastComputeLCM(const int &num1, const int &num2);
Rcpp::Date FastGetCurrentDate();
void PerformRVector(const arma::uvec &vec, const int &total, arma::uvec &new_vec);
void PerformRVector(const arma::urowvec &vec, const int &total, arma::urowvec &new_vec);
void PerformRIndex(const int &i, const int &total, int &new_i);
void PerformRMultiIndex(const int &start, const int &total_start, int &new_start, const int &end, const int &total_end, int &new_end);
arma::uvec FastDiffVector(const arma::uvec& a, const arma::uvec& b);
arma::uvec FastRandVector(int num);

#endif //COMMON_MATRIX_UTIL
