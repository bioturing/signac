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

double LnPvalue(double score, int n1, int n2, int bin_cnt)
{
    if (bin_cnt <= 0)
        throw std::domain_error("Bin count should be positive");

    if (bin_cnt == 1)
        return 0.0;

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

void GetTotalCount(
    const Rcpp::NumericVector &cluster,
    std::array<int, 2> &total_cnt)
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

// Expression matrix needs to be sorted
double ComputeUDScore(
        const Rcpp::NumericVector &cluster,
        const std::vector<std::pair<double, int>> &exp,
        const std::array<int, 2> &cnt,
        const std::array<int, 2> &zero_cnt)
{
    double up = 0;
    double down = 0;

    const int in = cnt[0];
    const int out = cnt[1];

    int c_in = 0;
    int c_out = 0;
    int e_in = zero_cnt[0];
    int e_out = zero_cnt[1];

    double prev_exp = 0;
    int n = exp.size();
    for (int i = 0; i < n; ++i) {
        double e = exp[i].first;
        if (std::fabs(e - prev_exp) > HARMONY_EPS) {
            if (e_in > 0) {
                double l = (double) e_in/in;

                double in_l = (double) c_in/in;
                double in_r = (double) (c_in + e_in)/in;

                double out_l = (double) c_out/out;
                double out_r = (double) (c_out + e_out)/out;

                up += std::min(l, std::max(0.0, out_l - in_l));
                down += std::min(l, std::max(0.0, in_r - out_r));
            }

            prev_exp = e;

            c_in += e_in;
            c_out += e_out;
            e_in = e_out = 0;
        }

        int group = exp[i].second;

        e_in += group == 1;
        e_out += group == 2;
    }

    up += std::max(0.0, (double)c_out/out - (double)c_in/in);

    return (up - down) / (up + down);
}

// Expression matrix needs to be sorted
std::tuple<double, double> ComputeSimilarity(
        const Rcpp::NumericVector &cluster,
        const std::vector<std::pair<double, int>> &exp,
        const std::array<int, 2> &cnt,
        const std::array<int, 2> &zero_cnt,
        int thres)
{
    const int in = cnt[0];
    const int out = cnt[1];
    const int zero_in = zero_cnt[0];
    const int zero_out = zero_cnt[1];

    int bin_in = zero_in;
    int bin_out = zero_out;
    int bin_both = bin_in + bin_out;

    int bin_cnt = 0;
    double sim = 0;

    double prev_exp = 0;

    int n = exp.size();
    for (int i = 0; i < n; ++i) {
        double e = exp[i].first;
        if (std::fabs(e - prev_exp) > HARMONY_EPS) {
            if (bin_both >= thres) {
                sim += HarmonicMean((double)bin_in/in,
                                    (double)bin_out/out);

                bin_in = bin_out = bin_both = 0;
                ++bin_cnt;
            }
            prev_exp = e;
        }

        int group = exp[i].second;
        bin_in += group == 1;
        bin_out += group == 2;
        ++bin_both;
    }

    if (bin_both >= MINIMAL_SAMPLE) {
        sim += HarmonicMean((double)bin_in/in,
                            (double)bin_out/out);
        ++bin_cnt;
    } else {
        sim += ((double) bin_in / in
                + (double) bin_out / out) / 2;
    }

    return std::make_tuple(
        sim,
        LnPvalue(sim, in, out, bin_cnt)
    );
}

std::tuple<double, double, double> ProcessGene(
        const Rcpp::NumericVector &cluster,
        std::vector<std::pair<double, int>> exp,
        const std::array<int, 2> &cnt,
        const std::array<int, 2> &zero_cnt,
        int thres)
{
    std::sort(exp.begin(), exp.end());

    if (!exp.empty() && exp[0].first < HARMONY_EPS)
        throw std::domain_error("Zero expression should not "
                                "be included in exp matrix");

    return std::tuple_cat(
        ComputeSimilarity(cluster, exp, cnt, zero_cnt, thres),
        std::make_tuple(ComputeUDScore(cluster, exp, cnt, zero_cnt))
    );
}

std::vector<std::tuple<double, double, double>> HarmonyTest(
        const arma::sp_mat &mtx,
        const Rcpp::NumericVector &cluster,
        const std::array<int, 2> &total_cnt)
{
    int thres = GetThreshold(total_cnt[0] + total_cnt[1]);
    int n_genes = mtx.n_rows;

    if (cluster.size() != mtx.n_cols)
        throw std::domain_error("Input cluster size is not equal "
                                "to the number of columns in matrix");

    std::vector<std::vector<std::pair<double, int>>> exp(n_genes);
    std::vector<std::array<int, 2>> nz_cnt(n_genes);

    std::vector<std::tuple<double, double, double>> res(n_genes);
    arma::sp_mat::const_col_iterator c_it;

    for (int i = 0; i < mtx.n_cols; ++i) {
        if (!(int)cluster[i])
            continue;

        for (c_it = mtx.begin_col(i); c_it != mtx.end_col(i); ++c_it) {
            if (*c_it) {
                int r = c_it.row();
                exp[r].push_back({*c_it, (int)cluster[i]});
                ++nz_cnt[r][(int)cluster[i] - 1];
            }
        }
    }

    for (int i = 0; i < n_genes; ++i) {
        res[i] = ProcessGene(
            cluster,
            std::move(exp[i]),
            total_cnt,
            {total_cnt[0] - nz_cnt[i][0], total_cnt[1] - nz_cnt[i][1]},
            thres
        );

        if ((i + 1) % 1000 == 0)
        Rcout << "Processed " << i + 1 << " genes\r";
    }

    return res;
}

std::vector<std::tuple<double, double, double>> HarmonyTest(
        com::bioturing::Hdf5Util &oHdf5Util,
        HighFive::File *file,
        const Rcpp::NumericVector &cluster,
        const std::array<int, 2> &total_cnt)
{
    int thres = GetThreshold(total_cnt[0] + total_cnt[1]);
    std::vector<int> shape;
    oHdf5Util.ReadDatasetVector<int>(file, GROUP_NAME, "shape", shape);

    int n_genes = shape[1];

    if (cluster.size() != shape[0])
        throw std::domain_error("Input cluster size is not equal to "
                                "the number of columns in matrix");

    std::vector<std::tuple<double, double, double>> res(n_genes);

    for (int i = 0; i < n_genes; ++i) {
        std::vector<int> col_idx;
        std::vector<double> g_exp;
        oHdf5Util.ReadGeneExpH5(file, GROUP_NAME, i, col_idx,  g_exp);

        std::vector<std::pair<double, int>> exp(col_idx.size());
        std::array<int, 2> zero_cnt = {total_cnt[0], total_cnt[1]};

        for (int k = 0; k < col_idx.size(); ++k) {
            int idx = (int)cluster[col_idx[k]];
            if (idx && g_exp[k]) {
                exp[k] = {g_exp[k], idx};
                --zero_cnt[idx - 1];
            }
        }

        res[i] = ProcessGene(
            cluster,
            std::move(exp),
            total_cnt,
            zero_cnt,
            thres
        );
    }

    return res;
}

DataFrame PostProcess(
        std::vector<std::tuple<double, double, double>> &res,
        std::vector<std::string> &rownames)
{
    int n_gene = res.size();
    std::vector<std::pair<double,int>> order(n_gene);

    for (int i = 0; i < n_gene; ++i)
        order[i] = std::make_pair(std::get<1>(res[i]), i);

    std::sort(order.begin(), order.end());

    std::vector<double> p_value(n_gene), similarity(n_gene);
    std::vector<double> p_adjusted(n_gene);
    std::vector<std::string> g_names(n_gene);
    std::vector<double> change(n_gene);

    //Adjust p value
    double prev = 0;
    for(int i = 0; i < n_gene; ++i) {
        double p = std::exp(order[i].first) * n_gene / (i + 1);

        if (p > 1)
            p = 1;

        if (p >= prev)
            prev = p;

        p_adjusted[i] = prev;
    }

    for(int i = 0; i < n_gene; ++i) {
        int k = order[i].second;
        g_names[i] = rownames[k];
        std::tie(similarity[i], p_value[i], change[i]) = res[k];
    }

    Rcout << "Done all" << std::endl;
    return DataFrame::create( Named("Gene Name") = wrap(g_names),
                              Named("Similarity") = wrap(similarity),
                              Named("Log P value") = wrap(p_value),
                              Named("P-adjusted value") = wrap(p_adjusted),
                              Named("Up-Down score") = wrap(change)
                            );
}

//' HarmonyMarker
//'
//' Find gene marker for a cluster in sparse matrix
//'
//' @param S4_mtx A sparse matrix
//' @param cluster A numeric vector
//' @export
// [[Rcpp::export]]
DataFrame HarmonyMarker(const Rcpp::S4 &S4_mtx, const Rcpp::NumericVector &cluster)
{
    Rcout << "Enter" << std::endl;

    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    const arma::sp_mat &mtx = Rcpp::as<arma::sp_mat>(S4_mtx);
    Rcout << "Done parse" << std::endl;

    std::vector<std::tuple<double, double, double>> res
                    = HarmonyTest(mtx, cluster, total_cnt);
    Rcout << "Done calculate" << std::endl;

    Rcpp::List dim_names = Rcpp::List(S4_mtx.attr("Dimnames"));
    std::vector<std::string> rownames = dim_names[0];
    return PostProcess(res, rownames);
}

//' HarmonyMarkerH5
//'
//' Find gene marker for a cluster in H5 file
//'
//' @param hdf5Path A string path
//' @param cluster A numeric vector
//' @export
// [[Rcpp::export]]
DataFrame HarmonyMarkerH5(
    const std::string &hdf5Path,
    const Rcpp::NumericVector &cluster)
{
    com::bioturing::Hdf5Util oHdf5Util(hdf5Path);
    HighFive::File *file = oHdf5Util.Open(1);

    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    Rcout << "Group1 " << total_cnt[0]
          << "Group2 " << total_cnt[1] << std::endl;

    std::vector<std::tuple<double, double, double>> res
        = HarmonyTest(oHdf5Util, file, cluster, total_cnt);

    Rcout << "Done calculate" << std::endl;
    std::vector<std::string> rownames;
    // Read the barcode slot since this is the transposed matrix
    oHdf5Util.ReadDatasetVector<std::string>(file, GROUP_NAME,
                                            "barcodes", rownames);
    oHdf5Util.Close(file);

    return PostProcess(res, rownames);
}
