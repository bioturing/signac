#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#include "CommonUtil.h"
//#include <boost/math/common_factor.hpp>
#include <boost/date_time/gregorian/gregorian_types.hpp>

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]

/*
Source:
http://www.unicode.org/versions/Unicode7.0.0/UnicodeStandard-7.0.pdf
*/
size_t IsUtf8Str(unsigned char *str, size_t len, int *faulty_bytes)
{
    size_t i = 0;
    *faulty_bytes = 0;

    while (i < len)
    {
        if (str[i] <= 0x7F) /* 00..7F */
        {
            i += 1;
        }
        else if (str[i] >= 0xC2 && str[i] <= 0xDF) /* C2..DF 80..BF */
        {
            if (i + 1 < len) /* Expect a 2nd byte */
            {
                if (str[i + 1] < 0x80 || str[i + 1] > 0xBF)
                {
                    *faulty_bytes = 2;
                    return i;
                }
            }
            else
            {
                *faulty_bytes = 1;
                return i;
            }
            i += 2;
        }
        else if (str[i] == 0xE0) /* E0 A0..BF 80..BF */
        {
            if (i + 2 < len) /* Expect a 2nd and 3rd byte */
            {
                if (str[i + 1] < 0xA0 || str[i + 1] > 0xBF)
                {
                    *faulty_bytes = 2;
                    return i;
                }
                if (str[i + 2] < 0x80 || str[i + 2] > 0xBF)
                {
                    *faulty_bytes = 3;
                    return i;
                }
            }
            else
            {
                *faulty_bytes = 1;
                return i;
            }
            i += 3;
        }
        else if (str[i] >= 0xE1 && str[i] <= 0xEC) /* E1..EC 80..BF 80..BF */
        {
            if (i + 2 < len) /* Expect a 2nd and 3rd byte */
            {
                if (str[i + 1] < 0x80 || str[i + 1] > 0xBF)
                {
                    *faulty_bytes = 2;
                    return i;
                }
                if (str[i + 2] < 0x80 || str[i + 2] > 0xBF)
                {
                    *faulty_bytes = 3;
                    return i;
                }
            }
            else
            {
                *faulty_bytes = 1;
                return i;
            }
            i += 3;
        }
        else if (str[i] == 0xED) /* ED 80..9F 80..BF */
        {
            if (i + 2 < len) /* Expect a 2nd and 3rd byte */
            {
                if (str[i + 1] < 0x80 || str[i + 1] > 0x9F)
                {
                    *faulty_bytes = 2;
                    return i;
                }
                if (str[i + 2] < 0x80 || str[i + 2] > 0xBF)
                {
                    *faulty_bytes = 3;
                    return i;
                }
            }
            else
            {
                *faulty_bytes = 1;
                return i;
            }
            i += 3;
        }
        else if (str[i] >= 0xEE && str[i] <= 0xEF) /* EE..EF 80..BF 80..BF */
        {
            if (i + 2 < len) /* Expect a 2nd and 3rd byte */
            {
                if (str[i + 1] < 0x80 || str[i + 1] > 0xBF)
                {
                    *faulty_bytes = 2;
                    return i;
                }
                if (str[i + 2] < 0x80 || str[i + 2] > 0xBF)
                {
                    *faulty_bytes = 3;
                    return i;
                }
            }
            else
            {
                *faulty_bytes = 1;
                return i;
            }
            i += 3;
        }
        else if (str[i] == 0xF0) /* F0 90..BF 80..BF 80..BF */
        {
            if (i + 3 < len) /* Expect a 2nd, 3rd 3th byte */
            {
                if (str[i + 1] < 0x90 || str[i + 1] > 0xBF)
                {
                    *faulty_bytes = 2;
                    return i;
                }
                if (str[i + 2] < 0x80 || str[i + 2] > 0xBF)
                {
                    *faulty_bytes = 3;
                    return i;
                }
                if (str[i + 3] < 0x80 || str[i + 3] > 0xBF)
                {
                    *faulty_bytes = 4;
                    return i;
                }
            }
            else
            {
                *faulty_bytes = 1;
                return i;
            }
            i += 4;
        }
        else if (str[i] >= 0xF1 && str[i] <= 0xF3) /* F1..F3 80..BF 80..BF 80..BF */
        {
            if (i + 3 < len) /* Expect a 2nd, 3rd 3th byte */
            {
                if (str[i + 1] < 0x80 || str[i + 1] > 0xBF)
                {
                    *faulty_bytes = 2;
                    return i;
                }
                if (str[i + 2] < 0x80 || str[i + 2] > 0xBF)
                {
                    *faulty_bytes = 3;
                    return i;
                }
                if (str[i + 3] < 0x80 || str[i + 3] > 0xBF)
                {
                    *faulty_bytes = 4;
                    return i;
                }
            }
            else
            {
                *faulty_bytes = 1;
                return i;
            }
            i += 4;
        }
        else if (str[i] == 0xF4) /* F4 80..8F 80..BF 80..BF */
        {
            if (i + 3 < len) /* Expect a 2nd, 3rd 3th byte */
            {
                if (str[i + 1] < 0x80 || str[i + 1] > 0x8F)
                {
                    *faulty_bytes = 2;
                    return i;
                }
                if (str[i + 2] < 0x80 || str[i + 2] > 0xBF)
                {
                    *faulty_bytes = 3;
                    return i;
                }
                if (str[i + 3] < 0x80 || str[i + 3] > 0xBF)
                {
                    *faulty_bytes = 4;
                    return i;
                }
            }
            else
            {
                *faulty_bytes = 1;
                return i;
            }
            i += 4;
        }
        else
        {
            *faulty_bytes = 1;
            return i;
        }
    }
    return 0;
}

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

/*
// [[Rcpp::export]]
int FastComputeGCD(const int &num1, const int &num2) {
    return boost::math::gcd(num1, num2);
}

// [[Rcpp::export]]
int FastComputeLCM(const int &num1, const int &num2) {
    return boost::math::lcm(num1, num2);
}
*/

// [[Rcpp::export]]
Rcpp::Date FastGetCurrentDate() {
    boost::gregorian::date current_date(boost::gregorian::day_clock::local_day());
    return Rcpp::wrap(Rcpp::Date(current_date.year(), current_date.month(), current_date.day()));
}

// [[Rcpp::export]]
arma::uvec FastDiffVector(const arma::uvec& a, const arma::uvec& b) {
    std::vector<int> x = arma::conv_to< std::vector<int> >::from(arma::sort(a));
    std::vector<int> y = arma::conv_to< std::vector<int> >::from(arma::sort(b));
    std::vector<int> z;
    std::set_difference(x.begin(), x.end(), y.begin(), y.end(),
                        std::inserter(z, z.end()));
    return arma::conv_to<arma::uvec>::from(z);
}

// [[Rcpp::export]]
arma::uvec FastRandVector(int num) {
    arma::uvec result(num);
    for (int i=0; i < num; ++i) {
        result(i) = i;
    }
    std::random_shuffle(result.begin(), result.end());
    return result;
}
