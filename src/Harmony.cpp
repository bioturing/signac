#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#define HARMONY_RATIO 1e-5
#define HARMONY_EPS 1e-50
#define MINIMAL_SAMPLE 10
#define GROUPING_RATE 0.6
#define GROUP_NAME "bioturing"

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
#include "SparseMatrixUtil.h"
#include "incbeta.h"
#include "Hdf5Util.h"

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

void GetTotalCount(const Rcpp::NumericVector &cluster, std::array<int, 2> &total_cnt)
{
    total_cnt[0] = total_cnt[1] = 0;

    for (int i = 0; i < cluster.size(); ++i)
        if ((int)cluster[i])
            ++total_cnt[(int)cluster[i] - 1];
}

int GetThreshold(int total)
{
    int thres = std::max(std::min(total, MINIMAL_SAMPLE),
                        (int) std::pow(total,GROUPING_RATE));

    Rcout << thres << std::endl;

    if(thres * 2 > total)
        throw std::runtime_error("Not enough bins to compare. "
                                "Please choose larger clusters to compare");
    return thres;
}

std::tuple<double, double, double> ComputeSimilarity(
        const Rcpp::NumericVector &cluster,
        std::vector<std::pair<double, int>> &exp,
        std::array<int, 2> total_cnt,
        std::array<int, 2> zero_cnt,
        int thres)
{
    int total_in = total_cnt[0];
    int total_out = total_cnt[1];
    int zero_in = zero_cnt[0];
    int zero_out = zero_cnt[1];
    int cnt_in = zero_in;
    int cnt_out = zero_out;
    int cnt_both = cnt_in + cnt_out;

    int bin_cnt = 0;
    double sim = 0;



    std::sort(exp.begin(), exp.end());

    for(int i = 0, l = exp.size(); i < l; ++i) {
        if ((!i ||  std::fabs(exp[i].first - exp[i - 1].first) > HARMONY_EPS)) {
            if (cnt_both >= thres) {
                sim += HarmonicMean((double)cnt_in/total_in,
                                    (double)cnt_out/total_out);

                cnt_in = cnt_out = cnt_both = 0;
                ++bin_cnt;
            }
        }

        cnt_in += exp[i].second == 1;
        cnt_out += exp[i].second == 2;
        ++cnt_both;
    }

    if (cnt_both >= MINIMAL_SAMPLE) {
        sim += HarmonicMean((double)cnt_in/total_in,
                            (double)cnt_out/total_out);
        ++bin_cnt;
    } else {
        sim += ((double) cnt_in / total_in
                + (double) cnt_out / total_out) / 2.0;
    }

    int cnt_up = 0;
    int cnt_down = 0;

    int x = -zero_in - 1;
    int y = -zero_out - 1;
    int z = -zero_out - 1;
    int cnt_x = -1;
    int cnt_y = -1;
    int cnt_z = -1;

    for (int i = 0; i < total_in; ++i) {
        double j = (double) i / (total_in - 1) * (total_out - 1);
        int j_down =  std::floor(j);
        int j_up = std::ceil(j);


        while (cnt_x < i) {
            ++x;
            cnt_x += x < 0 || exp[x].second == 1;
        }


        while (cnt_y < j_down) {
            ++y;
            cnt_y += y < 0 || exp[y].second == 2;
        }

        while (cnt_z < j_up) {
            ++z;
            cnt_z += z < 0 || exp[z].second == 2;
        }

        double i_exp = (x < 0? 0 : exp[x].first);
        double j_down_exp = (y < 0? 0 : exp[y].first);
        double j_up_exp = (z < 0? 0 : exp[z].first);

        double j_exp = (j_up == j_down?j_down_exp: (j_up_exp *(j - j_down) + j_down_exp * (j_up - j)));

        cnt_up += i_exp > j_exp;
        cnt_down += i_exp < j_exp;
    }

    return std::make_tuple(
        sim,
        LnPvalue(sim, total_in, total_out, bin_cnt),
        1.0 * (cnt_up - cnt_down) / (cnt_up + cnt_down)
    );
}

void HarmonyTest(
        const arma::sp_mat &mtx,
        const Rcpp::NumericVector &cluster,
        std::vector<std::tuple<double, double, double> > &res,
        std::array<int, 2> total_cnt)
{
    int thres = GetThreshold(total_cnt[0] + total_cnt[1]);
    int n_genes = mtx.n_rows;

    if (cluster.size() != mtx.n_cols)
        throw std::domain_error("Input cluster size is not equal to the number of rows in matrix");

    std::vector<std::vector<std::pair<double, int>>> exp(n_genes);
    std::vector<std::array<int, 2>> zero_cnt(n_genes);

    res.resize(n_genes);
    arma::sp_mat::const_col_iterator c_it;

    for (int i = 0; i < mtx.n_cols; ++i) {
        if (!(int)cluster[i])
            continue;

        for (c_it = mtx.begin_col(i); c_it != mtx.end_col(i); ++c_it) {
            if (*c_it) {
                int r = c_it.row();
                exp[r].push_back({*c_it, (int)cluster[i]});
                ++zero_cnt[r][(int)cluster[i] - 1];
            }
        }
    }

    for (int i = 0; i < n_genes; ++i) {
        zero_cnt[i][0] = total_cnt[0] - zero_cnt[i][0];
        zero_cnt[i][1] = total_cnt[1] - zero_cnt[i][1];
        res[i] = ComputeSimilarity(cluster, exp[i], total_cnt, zero_cnt[i], thres);
    }
}

void HarmonyTest(
        com::bioturing::Hdf5Util &oHdf5Util,
        HighFive::File *file,
        const Rcpp::NumericVector &cluster,
        std::vector<std::tuple<double, double, double> > &res,
        std::array<int, 2> total_cnt)
{
    int thres = GetThreshold(total_cnt[0] + total_cnt[1]);
    std::vector<int> shape;
    oHdf5Util.ReadDatasetVector<int>(file, GROUP_NAME, "shape", shape);

    int n_genes = shape[1];

    if (cluster.size() != shape[0])
        throw std::domain_error("Input cluster size is not equal to the number of rows in matrix");

    res.resize(n_genes);

    for (int i = 0; i < n_genes; ++i) {
        std::vector<int> col_idx;
        std::vector<double> g_exp;
        oHdf5Util.ReadGeneExpH5(file, GROUP_NAME, i, col_idx,  g_exp);

        std::vector<std::pair<double, int>> exp;
        std::array<int, 2> zero_cnt = {total_cnt[0], total_cnt[1]};

        for (int k = 0; k < col_idx.size(); ++k) {
            int idx = (int)cluster[col_idx[k]];
            if (idx && g_exp[k]) {
                exp.push_back({g_exp[k], idx});
                --zero_cnt[idx - 1];
            }
        }

        res[i] = ComputeSimilarity(cluster, exp, total_cnt, zero_cnt, thres);
    }
}

DataFrame PostProcess(
        std::vector<std::tuple<double, double, double>> &res,
        std::vector<std::string> &rownames)
{
    int n_gene = res.size();
    std::vector<int> order(n_gene);
    std::vector<double> pvalue(n_gene);

    for (int i = 0; i < n_gene; ++i) {
        order[i] = i;
        pvalue[i] = std::get<1>(res[i]);
    }

    std::sort(order.begin(), order.end(),
              std::bind(Compare, _1, _2, pvalue));

    std::vector<double> p_value(n_gene), similarity(n_gene);
    std::vector<double> p_adjusted(n_gene);
    std::vector<std::string> g_names(n_gene);
    std::vector<double> diff_class(n_gene);


    //Adjust p value
    double prev = 0;
    for(int i = 0; i < n_gene; ++i) {
        int k = order[i];
        double p = std::exp(pvalue[k]) * n_gene / (i + 1);

        if (p > 1)
            p = 1;

        if (p >= prev)
            prev = p;

        p_adjusted[i] = prev;
    }

    for(int i = 0; i < n_gene; ++i) {
        int k = order[i];
        g_names[i] = rownames[k];
        similarity[i] = std::get<0>(res[k]);
        p_value[i] = std::get<1>(res[k]);
        diff_class[i] = std::get<2>(res[k]);
    }

    Rcout << "Done all" << std::endl;
    return DataFrame::create( Named("Gene Name") = wrap(g_names),
                              Named("Similarity") = wrap(similarity),
                              Named("Log P value") = wrap(p_value),
                              Named("P-adjusted value") = wrap(p_adjusted),
                              Named("Type") = wrap(diff_class)
                            );
}

// [[Rcpp::export]]
DataFrame HarmonyMarker(Rcpp::S4 &S4_mtx, const Rcpp::NumericVector &cluster)
{
    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    const arma::sp_mat mtx = Rcpp::as<arma::sp_mat>(S4_mtx);
    std::vector<std::tuple<double, double, double>> res;

    HarmonyTest(mtx, cluster, res, total_cnt);
    Rcout << "Done calculate" << std::endl;

    Rcpp::List dim_names = Rcpp::List(S4_mtx.attr("Dimnames"));
    std::vector<std::string> rownames = dim_names[0];
    return PostProcess(res, rownames);
}

// [[Rcpp::export]]
DataFrame HarmonyMarkerH5(const std::string &hdf5Path, const Rcpp::NumericVector &cluster)
{
    com::bioturing::Hdf5Util oHdf5Util(hdf5Path);
    HighFive::File *file = oHdf5Util.Open(1);

    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    Rcout << "Group1 " << total_cnt[0] << "Group2 " << total_cnt[1] << std::endl;

    std::vector<std::tuple<double, double, double>> res;

    HarmonyTest(oHdf5Util, file, cluster, res, total_cnt);
    Rcout << "Done calculate" << std::endl;
    std::vector<std::string> rownames;
    oHdf5Util.ReadDatasetVector<std::string>(file, GROUP_NAME, "features/name", rownames);

    oHdf5Util.Close(file);

    return PostProcess(res, rownames);
}
