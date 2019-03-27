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

#endif //COMMON_MATRIX_UTIL
