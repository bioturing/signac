#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#define HARMONY_RATIO 1e-5
#define HARMONY_EPS 1e-50

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>
#include <functional>
#include <assert.h>

#include "CommonUtil.h"
#include "incbeta.h"

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

// [[Rcpp::export]]
void HarmonyFit(const arma::sp_mat &mtx,
                               const std::vector<std::pair<std::string, int>> &cluster,
                               int cid,
                               std::vector<std::pair<double,double> > &res)
{
    res.resize(mtx.n_rows, std::make_pair(0.0,0.0));
    arma::sp_mat::const_col_iterator c_it;

    for (int i = 0; i < mtx.n_cols; ++i) {
        for (c_it = mtx.begin_col(i); c_it != mtx.end_col(i); ++c_it) {
            if (cluster[i].second != cid)
                continue;

            res[c_it.row()].first += (*c_it);
        }
    }

    int total = 0;
    for (int i = 0 ; i < mtx.n_cols; ++i)
        total += (cluster[i].second == cid);

    for (int i = 0; i < mtx.n_rows; ++i)
        res[i].first /= total;

    std::vector<int> zero;
    zero.resize(mtx.n_rows, total);

    for (int i = 0; i < mtx.n_cols; ++i) {
        for (c_it = mtx.begin_col(i); c_it != mtx.end_col(i); ++c_it) {
            if (cluster[i].second != cid)
                continue;
            res[c_it.row()].second += std::pow((*c_it) - res[c_it.row()].first, 2);
            --zero[c_it.row()];
        }
    }

    for (int i = 0; i < mtx.n_rows; ++i) {
        res[i].second += zero[i] * std::pow(res[i].first,2);
        res[i].second /= total - 1;

        if (res[i].second < res[i].first * (1 + HARMONY_RATIO))
            res[i].second = res[i].first * (1 + HARMONY_RATIO);
    }
}
