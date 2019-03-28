#ifndef COMMON_MATRIX_UTIL
#define COMMON_MATRIX_UTIL

#include <RcppArmadillo.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <boost/math/common_factor.hpp>
#include <boost/date_time/gregorian/gregorian_types.hpp>
using namespace Rcpp;

int FastComputeGCD(const int &num1, const int &num2);
int FastComputeLCM(const int &num1, const int &num2);
Rcpp::Date FastGetCurrentDate();
void PerformRVector(const arma::uvec &vec, const int &total, arma::uvec &new_vec);
void PerformRVector(const arma::urowvec &vec, const int &total, arma::urowvec &new_vec);
void PerformRIndex(const int &i, const int &total, int &new_i);
void PerformRMultiIndex(const int &start, const int &total_start, int &new_start, const int &end, const int &total_end, int &new_end);

#endif //COMMON_MATRIX_UTIL
