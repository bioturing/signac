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

double HarmonicMean(double a, double b) {
    if (fmin(a,b) < HARMONY_EPS)
        return 0;

    return 2 * a * b / (a + b);
}

bool Compare(int i, int j, const std::vector<double> &res) {
    return res[i] < res[j];
}

void SetCluster(const Rcpp::NumericVector &col_idx,
                  std::vector<int> &result)
{
    //result.fill(2); //default is 2
    for (int i = 0; i < col_idx.length(); ++i) {
        assert(i < (result.size() - 1) && "Index must be smaller than number of columns");
        result[col_idx[i] - 1] = 1; //Convert
    }
}

double LnPvalue(double score, int n1, int n2, int bin_cnt) {
        --bin_cnt;
        assert(bin_cnt >= 0);

        double mean = 0.25 * bin_cnt * (1.0/n1 + 1.0/n2);
        assert(mean >= 0);
        double var = 0.125 * bin_cnt * pow(1.0/n1 + 1.0/n2,2);
        double scale =  (mean * (1 - mean) / var - 1);
        double shape1 = (1 - mean) * scale;
        double shape2 = mean * scale;

        assert(shape1 > 0 && shape2 > 0);

        return logincbeta(shape1, shape2, score);
}

void HarmonyTest(
        const arma::sp_mat &mtx,
        const std::vector<int> &cluster,
        std::vector<std::tuple<double, double> > &res)
{
    std::vector<std::vector<std::pair<int,int>>> count(mtx.n_rows);

    res.resize(mtx.n_rows, std::make_tuple(0.0, 1.0));
    arma::sp_mat::const_col_iterator c_it;

    for(int i = 0; i < mtx.n_cols; ++i) {
        for (c_it = mtx.begin_col(i); c_it != mtx.end_col(i); ++c_it) {
            int cnt = (*c_it);

            std::vector<std::pair<int, int>> &r_cnt = count[c_it.row()];

            if (cnt >= r_cnt.size())
                r_cnt.resize(cnt + 1);

            if (cluster[i] == 1)
                ++r_cnt[cnt].first;
            else if (cluster[i] == 2)
                ++r_cnt[cnt].second;
        }
    }

    int total_in = 0;
    int total_out = 0;

    for (int i = 0; i < mtx.n_cols; ++i) {
        total_in += cluster[i] == 1;
        total_out += cluster[i] == 2;
    }

    int total = total_in + total_out;

    for (int i = 0; i < mtx.n_rows; ++i) {
        int zero_in = total_in;
        int zero_out = total_out;

        for(int j = 1; j < count[i].size(); ++j) {
            zero_in -= count[i][j].first;
            zero_out -= count[i][j].second;
        }

        if (count[i].size() < 1)
            count[i].resize(1);

        count[i][0] = std::make_pair(zero_in, zero_out);
    }

    double thres = std::max(std::min<double>(total,5.0),std::pow(total,0.6))/total;

    assert(thres < 0.5);

    for (int i = 0; i < mtx.n_rows; ++i) {

        double prob_in = 0;
        double prob_out = 0;
        double prob = 0;
        double prob_rem = 1;

        int bin_cnt = 1; //bin_count actually

        for(int j = 0; j < count[i].size(); ++j) {
            prob_in += (double)count[i][j].first/total_in;
            prob_out += (double)count[i][j].second/total_out;

            double p = (double)(count[i][j].first + count[i][j].second)/total;
            prob += p;
            prob_rem -= p;

            if (prob >= thres && prob_rem >= thres) {
                std::get<0>(res[i]) += HarmonicMean(prob_in, prob_out);
                prob_in = prob_out = prob = 0;
                ++bin_cnt;
            }
        }

        std::get<0>(res[i]) += HarmonicMean(prob_in, prob_out);
        std::get<1>(res[i]) = LnPvalue(std::get<0>(res[i]), total_in, total_out, bin_cnt);
    }
}

// [[Rcpp::export]]
DataFrame HarmonyMarker(const arma::sp_mat &mtx,
            const Rcpp::NumericVector &col_idx)
{
    std::vector<int> cluster(mtx.n_cols, 2);
    SetCluster(col_idx, cluster);

    std::vector<std::tuple<double, double>> res;
    HarmonyTest(mtx, cluster, res);
    Rcout << "Done calculate" << std::endl;

    std::vector<int> order(mtx.n_rows);
    std::vector<double> pvalue(mtx.n_rows);

    for (int i = 0; i < mtx.n_rows; ++i) {
        order[i] = i;
        pvalue[i] = std::get<1>(res[i]);
    }

    std::sort(order.begin(), order.end(), std::bind(Compare, std::placeholders::_1, std::placeholders::_2, pvalue));

    std::vector<double> p_value(mtx.n_rows), similatiry(mtx.n_rows);
    std::vector<double> p_adjusted(mtx.n_rows);
    std::vector<int> g_index(mtx.n_rows);

    //Adjust p value
    double prev = 0;
    for(int i = 0; i < mtx.n_rows; ++i) {
        int k = order[i];
        double p = std::exp(pvalue[k]) * mtx.n_rows / (i + 1);

        if (p > 1)
            p = 1;

        if (p >= prev)
            prev = p;

        p_adjusted[i] = prev;
    }

    for(int i = 0; i < mtx.n_rows; ++i) {
        int k = order[i];
        g_index[i] = k + 1;
        similatiry[i] = std::get<0>(res[k]);
        p_value[i] = std::get<1>(res[k]);
    }

    Rcout << "Done all" << std::endl;
    return DataFrame::create( Named("Gene Index") = wrap(g_index),
                              Named("Similarity") = wrap(similatiry),
                              Named("Log P value") = wrap(p_value),
                              Named("P-adjusted value") = wrap(p_adjusted));
}
