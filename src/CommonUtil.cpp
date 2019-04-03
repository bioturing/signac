#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#include "CommonUtil.h"
#include <boost/date_time/gregorian/gregorian_types.hpp>

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]

template <typename T>
void PerformRVector(const T &vec, const int &total, T &new_vec) {
    for(int i = 0; i < vec.size(); i++) {
        int ival = vec.at(i);
        if(ival <= 0) {
            std::stringstream ostr;
            ostr << "Invalid start index:" << ival;
            throw std::range_error(ostr.str());
        }

        if(ival > total) {
            std::stringstream ostr;
            ostr << "Invalid end index:" << ival;
            throw std::range_error(ostr.str());
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
        std::stringstream ostr;
        ostr << "Invalid index:" << i;
        throw std::range_error(ostr.str());
    }

    if(i > total) {
        std::stringstream ostr;
        ostr << "Invalid index:" << i;
        throw std::range_error(ostr.str());
    }

    new_i = new_i - 1;
}

void PerformRMultiIndex(const int &start, const int &total_start, int &new_start, const int &end, const int &total_end, int &new_end) {
    new_start = start;
    if(start <= 0) {
        std::stringstream ostr;
        ostr << "Invalid start index:" << start;
        throw std::range_error(ostr.str());
    }

    if(start > total_start) {
        std::stringstream ostr;
        ostr << "Invalid start index:" << start;
        throw std::range_error(ostr.str());
    }

    new_end = end;
    if(end <= 0) {
        std::stringstream ostr;
        ostr << "Invalid end index:" << start;
        throw std::range_error(ostr.str());
    }

    if(end > total_end) {
        std::stringstream ostr;
        ostr << "Invalid end index:" << start;
        throw std::range_error(ostr.str());
    }

    if(new_start > new_end) {
        std::stringstream ostr;
        ostr << "Invalid start/end index:" << start << ">" << end;
        throw std::range_error(ostr.str());
    }

    new_start = new_start -1;
    new_end = new_end - 1;
}

//' FastGetCurrentDate
//'
//' This function returns a current date (YYYY-MM-DD)
//'
//' @export
// [[Rcpp::export]]
Rcpp::Date FastGetCurrentDate() {
    boost::gregorian::date current_date(boost::gregorian::day_clock::local_day());
    return Rcpp::wrap(Rcpp::Date(current_date.year(), current_date.month(), current_date.day()));
}

//' FastDiffVector
//'
//' This function is used to diff 2 vector
//'
//' @param a An integer vector
//' @param b An integer vector
//' @export
// [[Rcpp::export]]
arma::uvec FastDiffVector(const arma::uvec& a, const arma::uvec& b) {
    std::vector<int> x = arma::conv_to< std::vector<int> >::from(arma::sort(a));
    std::vector<int> y = arma::conv_to< std::vector<int> >::from(arma::sort(b));
    std::vector<int> z;
    std::set_difference(x.begin(), x.end(), y.begin(), y.end(),
                        std::inserter(z, z.end()));
    return arma::conv_to<arma::uvec>::from(z);
}

//' FastRandVector
//'
//' This function create a random vector
//'
//' @param num An integer number
//' @export
// [[Rcpp::export]]
arma::uvec FastRandVector(int num) {
    arma::uvec result(num);
    for (int i=0; i < num; ++i) {
        result(i) = i;
    }
    std::random_shuffle(result.begin(), result.end());
    return result;
}
