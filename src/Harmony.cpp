#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#define HARMONY_RATIO 1e-5
#define HARMONY_EPS 1e-50
#define MINIMAL_SAMPLE 5
#define GROUPING_RATE 0.6

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

#include "CommonUtil.h"
#include "incbeta.h"

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace std::placeholders;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

double HarmonicMean(double a, double b)
{
    if (fmin(a,b) < HARMONY_EPS)
        return 0;

    return 2 * a * b / (a + b);
}

bool Compare(int i, int j, const std::vector<double> &res)
{
    return res[i] < res[j];
}

void SetCluster(std::vector<int> &cluster,
                const Rcpp::NumericVector &col_idx,
                int id)
{
    for (int i = 0; i < col_idx.length(); ++i) {
        if(col_idx[i] <= 0 || col_idx[i] > cluster.size())
            throw std::runtime_error("Index must be smaller than "
                                    "number of columns");

        if (cluster[col_idx[i] - 1])
            throw std::runtime_error("Each cell should only "
                                    "be assigned to one cluster");

        cluster[col_idx[i] - 1] = id;
    }
}

void FillRest(std::vector<int> &cluster, int id) {

    for (int i = 0; i < cluster.size(); ++i) {
        if (cluster[i] == 0)
            cluster[i] = id;
    }
}

double LnPvalue(double score, int n1, int n2, int bin_cnt)
{
    if (bin_cnt <= 0)
        throw std::domain_error("Bin count should be positive");

    if (bin_cnt == 1)
        return std::numeric_limits<double>::infinity();

    double mean = 0.25 * (bin_cnt - 1) * (1.0 / n1 + 1.0 / n2);

    if (mean > 1 || mean < 0) {
        throw std::runtime_error("Mean estimation lies outside [0,1] range. "
                                 "Cluster sizes may be too small.");
    }

    double var = 0.125 * (bin_cnt - 1) * pow(1.0 / n1 + 1.0 / n2,2);
    double scale =  (mean * (1 - mean) / var - 1);
    double shape1 = (1 - mean) * scale;
    double shape2 = mean * scale;

    return logincbeta(shape1, shape2, score);
}

std::string classify(double down, double mid, double up) {
    if (down > 0 && up < 0)
        return "Down";
    else if (up > 0 && down < 0)
        return "Up";
    else if (up > 0 && down > 0)
        return "Mid-down";
    else if (up < 0 && down < 0)
        return "Mid-up";

    return "Other";
}

void HarmonyTest(
        const arma::sp_mat &mtx,
        const std::vector<int> &cluster,
        std::vector<std::tuple<double, double, std::string> > &res)
{
    std::vector<std::vector<std::pair<int,int>>> count(mtx.n_rows);

    res.resize(mtx.n_rows);
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

    Rcout << total_in << " " << total_out << " " << total << std::endl;

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

    int thres = std::max(std::min(total, MINIMAL_SAMPLE),
                         (int) std::pow(total,GROUPING_RATE));

    Rcout << thres << std::endl;

    if(thres * 2 > total)
        throw std::runtime_error("Not enough bins to compare. "
                                "Please choose larger clusters to compare");

    for (int i = 0; i < mtx.n_rows; ++i) {

        int cnt_in = 0;
        int cnt_out = 0;
        int cnt_both = 0;

        int l = 0;
        int bin_cnt = 0;

        double up_cnt = 0;
        double down_cnt = 0;
        double mid_cnt = 0;

        double sim = 0;

        double mean_cnt = (double) total_in * total_in * total_out / 3;

        for(int j = 0; j < count[i].size(); ++j) {
            int a = count[i][j].first;
            int b = count[i][j].second;

            cnt_in += a;
            cnt_out += b;
            cnt_both += a + b;

            int r = total_in - l - a;
            double up  = 0.5 * (l + 1) * l + 0.5 * l * a + 1.0/6 * a * (a + 1);
            double down = 0.5 * (r + 1) * r + 0.5 * r * a + 1.0/6 * a * (a + 1);
            double mid = (total_in) * (total_in + 1) / 2 - up - down;
            l += a;

            up_cnt += up * b;
            down_cnt += down * b;
            mid_cnt += mid * b;

            if (cnt_both >= thres) {
                sim += HarmonicMean((double)cnt_in/total_in,
                                    (double)cnt_out/total_out);
                cnt_in = cnt_out = cnt_both = 0;
                ++bin_cnt;
            }
        }

        if (cnt_both >= MINIMAL_SAMPLE) {
            sim += HarmonicMean((double)cnt_in/total_in,
                                (double)cnt_out/total_out);
            ++bin_cnt;
        } else {
            sim += ((double) cnt_in / total_in
                  + (double) cnt_out / total_out) / 2.0;
        }

        res[i] = std::make_tuple(
                    sim,
                    LnPvalue(sim, total_in, total_out, bin_cnt),
                    classify(down_cnt - mean_cnt,
                            mid_cnt - mean_cnt,
                            up_cnt - mean_cnt
                    )
                );
    }
}

// [[Rcpp::export]]
DataFrame HarmonyMarker(const arma::sp_mat &mtx,
            const Rcpp::NumericVector &in_idx,
            const Rcpp::Nullable<Rcpp::NumericVector> &out_idx = R_NilValue)
{
    std::vector<int> cluster(mtx.n_cols);

    SetCluster(cluster, in_idx, 1);

    if (out_idx.isNull())
        FillRest(cluster, 2);
    else
        SetCluster(cluster, out_idx.get(), 2);


    std::vector<std::tuple<double, double, std::string>> res;
    HarmonyTest(mtx, cluster, res);
    Rcout << "Done calculate" << std::endl;

    std::vector<int> order(mtx.n_rows);
    std::vector<double> pvalue(mtx.n_rows);

    for (int i = 0; i < mtx.n_rows; ++i) {
        order[i] = i;
        pvalue[i] = std::get<1>(res[i]);
    }

    std::sort(order.begin(), order.end(),
              std::bind(Compare, _1, _2, pvalue));

    std::vector<double> p_value(mtx.n_rows), similatiry(mtx.n_rows);
    std::vector<double> p_adjusted(mtx.n_rows);
    std::vector<int> g_index(mtx.n_rows);
    std::vector<std::string> diff_class(mtx.n_rows);


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
        diff_class[i] = std::get<2>(res[k]);
    }

    Rcout << "Done all" << std::endl;
    return DataFrame::create( Named("Gene Index") = wrap(g_index),
                              Named("Similarity") = wrap(similatiry),
                              Named("Log P value") = wrap(p_value),
                              Named("P-adjusted value") = wrap(p_adjusted),
                              Named("Type") = wrap(diff_class)
                            );
}
