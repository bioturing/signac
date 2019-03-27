#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <boost/math/common_factor.hpp>
#include <boost/date_time/gregorian/gregorian_types.hpp>

using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]

// [[Rcpp::export]]
int FastComputeGCD(const int &num1, const int &num2) {
    return boost::math::gcd(num1, num2);
}

// [[Rcpp::export]]
int FastComputeLCM(const int &num1, const int &num2) {
    return boost::math::lcm(num1, num2);
}

// [[Rcpp::export]]
Rcpp::Date FastGetCurrentDate() {
    boost::gregorian::date current_date(boost::gregorian::day_clock::local_day());
    return Rcpp::wrap(Rcpp::Date(current_date.year(), current_date.month(), current_date.day()));
}
