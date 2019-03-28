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

template <typename T>
void PerformRVector(const T &vec, const int &total, T &new_vec) {
    for(int i = 0; i < vec.size(); i++) {
        int ival = vec.at(i);
        if(ival <= 0) {
            ival = 1;
        }

        if(ival > total) {
            ival = total;
        }
        new_vec(i) = (ival - 1);
    }
}

void PerformRVector(const arma::uvec &vec, const int &total, arma::uvec &new_vec) {
    PerformRVector<arma::uvec>(vec, total, new_vec);
}

void PerformRVector(const arma::urowvec &vec, const int &total, arma::urowvec &new_vec) {
    PerformRVector<arma::urowvec>(vec, total, new_vec);
}

void PerformRIndex(const int &i, const int &total, int &new_i) {
    new_i = i;
    if(i <= 0) {
        new_i = 1;
    }

    if(i > total) {
        new_i = total;
    }

    new_i = new_i - 1;
}

void PerformRMultiIndex(const int &start, const int &total_start, int &new_start, const int &end, const int &total_end, int &new_end) {
    new_start = start;
    if(start <= 0) {
        new_start = 1;
    }

    if(start > total_start) {
        new_start = total_start;
    }

    new_end = end;
    if(end <= 0) {
        new_end = 1;
    }

    if(end > total_end) {
        new_end = total_end;
    }

    if(new_start > new_end) {
        int item = new_end;
        new_end = new_start;
        new_start = item;
    }

    new_start = new_start -1;
    new_end = new_end - 1;
}

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
