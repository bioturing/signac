#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#define HARMONY_EPS 1e-50
#define MINIMAL_SAMPLE 10
#define MINIMAL_BIN 2
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
#include "Hdf5Util.h"
#include "chisq.h"

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace std::placeholders;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

struct GeneResult {
    int gene_id;
    std::string gene_name;

    double d_score; //dissimilarity score
    double b_cnt;

    double log_p_value;
    double perm_p_value;

    double ud_score; //up-down score

    double log_fc;
};

inline double HarmonicMean(double a, double b)
{
    return 2 / (1 / a + 1 / b);
}

inline double Score(int x, int y, int n1, int n2)
{
    if (x == 0 && y == 0)
        return 0;

    double a = (double)x / n1;
    double b = (double)y / n2;

    double result = pow(a-b,2)/(2*(a+b));
    //correction
    double correction = 2 * a * b * (n1 * a * (1 - b) + n2 * b *(1 - a))
                        / (n1 * n2 * pow(a + b, 3));

    return result - correction;
}

inline double LnPvalue(double score, int n1, int n2, int b_cnt)
{
    if (b_cnt <= 0)
        throw std::domain_error("Bin count should be positive");

    int Dof = b_cnt - 1;
    double x = 2 * score * HarmonicMean(n1, n2) + b_cnt - 1;
    return log_chisqr(Dof, x);
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

int GetThreshold(const std::array<int,2> &cnt)
{
    double avg = HarmonicMean(cnt[0],cnt[1]);
    double est_bin = std::pow(avg,1 - GROUPING_RATE);

    int thres = std::max<int>(
                        MINIMAL_SAMPLE,
                        ((cnt[0] + cnt[1]) / est_bin)
                );

    Rcout << cnt[0] << " " << cnt[1] << " " << avg << " " << thres << std::endl;

    if(thres * MINIMAL_BIN > cnt[0] + cnt[1])
        throw std::runtime_error("Not enough bins to compare. "
                                "Please choose larger groups to compare");
    return thres;
}

// Expression matrix needs to be sorted
double ComputeUd_score(
        const std::vector<std::array<int, 2>> &bins,
        const std::array<int, 2> &cnt)
{
    double up   = 0;
    double down = 0;

    const int in  = cnt[0];
    const int out = cnt[1];

    int cl_in  = 0;
    int cl_out = 0;
    int cr_in  = 0;
    int cr_out = 0;

    double prev_exp = 0;
    int n = bins.size();
    for (int i = 0; i < n; ++i) {
        cr_in = cl_in + bins[i][0];
        cr_out = cl_out + bins[i][1];

        double len = (double)bins[i][0] / in;

        if (len > 0) {
            double in_l = (double) cl_in/in;
            double in_r = (double) cr_in/in;

            double out_l = (double) cl_out/out;
            double out_r = (double) cr_out/out;

            up += std::min(len, std::max(0.0, out_l - in_l));
            down += std::min(len, std::max(0.0, in_r - out_r));
        }

        cl_in = cr_in;
        cl_out = cr_out;
    }

    return (up - down) / (up + down);
}


void Resample(
    std::vector<std::array<int, 2>> &bins,
    const std::vector<bool> &group,
    const std::array<int, 2> &cnt)
{
    int j = 0;
    for (int i = 0; i < bins.size(); ++i) {
        int next_j = j + bins[i][0] + bins[i][1];

        int x = 0;

        for (; j < next_j; ++j)
            x += group[j];

        int y = bins[i][0] + bins[i][1] - x;

        bins[i][0] = x;
        bins[i][1] = y;
    }
}

// Bins is the group count for each UMI value after sorted
double ComputeSimilarity(
        const std::vector<std::array<int, 2>> &bins,
        const std::array<int, 2> &cnt)
{
    double d_score = 0;

    int n1 = cnt[0];
    int n2 = cnt[1];

    int n = bins.size();

    int i = 0;

    for (; i < n; ++i) {
        int b1 = bins[i][0];
        int b2 = bins[i][1];

        d_score += Score(b1, b2, n1, n2);
    }

    return d_score;
}


std::vector<std::array<int, 2>> Binning(
        std::vector<std::pair<double, int>> exp,
        const std::array<int, 2> &zero_cnt)
{
    std::vector<std::array<int, 2>> result(exp.size() + 1);    //+1 for zero

    if (exp.size() == 0) {
        result[0] = zero_cnt;
        return result;
    }

    std::sort(exp.begin(), exp.end());

    double p_exp = exp[0].first;

    int i = 0, j = 0;
    for (; i < exp.size(); ++i) {
        double c_exp = exp[i].first;

        if (c_exp >= 0)
            break;

        if (abs(c_exp - p_exp) >= HARMONY_EPS) {
            ++j;
            p_exp = c_exp;
        }
        ++result[j][exp[i].second];
    }

    if (zero_cnt[0] + zero_cnt[1] > 0) {
        if (p_exp < -HARMONY_EPS) {
            ++j;
            p_exp = 0;
        }

        result[j][0] += zero_cnt[0];
        result[j][1] += zero_cnt[1];
    }

    for (; i < exp.size(); ++i) {
        double c_exp = exp[i].first;

        if (abs(c_exp - p_exp) >= HARMONY_EPS) {
            ++j;
            p_exp = c_exp;
        }
        ++result[j][exp[i].second];
    }

    result.resize(j + 1);
    return result;
}


void Grouping(
        std::vector<std::array<int, 2>> &bins,
        int thres)
{
    int b1 = 0;
    int b2 = 0;

    int n = bins.size();

    int i = 0, j = 0;

    for (; i < n; ++i) {
        b1 += bins[i][0];
        b2 += bins[i][1];

        if (b1 + b2 >= thres) {
            bins[j][0] = b1;
            bins[j][1] = b2;

            b1 = b2 = 0;
            ++j;
        }
    }

    if (b1 + b2 < MINIMAL_SAMPLE) {
        --j;
        bins[j][0] += b1;
        bins[j][1] += b2;
    } else {
        bins[j][0] = b1;
        bins[j][1] = b2;
    }

    bins.resize(j + 1);
}

void ProcessGene(
        std::vector<std::pair<double, int>> exp,
        const std::array<int, 2> &cnt,
        const std::array<int, 2> &zero_cnt,
        int thres,
        int perm,
        struct GeneResult &result)
{
    std::vector<std::array<int, 2>> bins = Binning(std::move(exp), zero_cnt);

    result.ud_score = ComputeUd_score(bins, cnt);
    Grouping(bins, thres);

    result.d_score = ComputeSimilarity(bins, cnt);
    result.b_cnt = bins.size();
    result.log_p_value = LnPvalue(result.d_score, cnt[0], cnt[1], bins.size());

    if (bins.size() <= 1) {
        result.perm_p_value = 1;
        return;
    }

    std::vector<bool> group(cnt[0] + cnt[1]);
    std::fill(group.begin(), group.begin() + cnt[0], true);

    int count = 0;
    for(int i = 0; i < perm; ++i) {
        std::random_shuffle(group.begin(), group.end());
        Resample(bins, group, cnt);
        double score = ComputeSimilarity(bins, cnt);
        count += score >= result.d_score;
    }

    result.perm_p_value = (double)count/perm;
}

std::vector<struct GeneResult> HarmonyTest(
        const arma::sp_mat &mtx,
        const Rcpp::NumericVector &cluster,
        const std::array<int, 2> &total_cnt,
        int threshold,
        int perm)
{
    int thres = threshold == 0? GetThreshold(total_cnt) : threshold;
    int n_genes = mtx.n_rows;

    if (cluster.size() != mtx.n_cols)
        throw std::domain_error("Input cluster size is not equal "
                                "to the number of columns in matrix");

    std::vector<std::vector<std::pair<double, int>>> exp(n_genes);

    std::vector<struct GeneResult> res(n_genes);
    arma::sp_mat::const_col_iterator c_it;

    std::vector<std::array<int, 2>> zero_cnt(n_genes, total_cnt);

    std::vector<std::array<double, 2>> total_exp(n_genes);


    for (int i = 0; i < mtx.n_cols; ++i) {
        if (!(int)cluster[i])
            continue;

        int cid = (int)cluster[i] - 1;

        for (c_it = mtx.begin_col(i); c_it != mtx.end_col(i); ++c_it) {
            int r = c_it.row();

            exp[r].push_back({*c_it, cid});
            total_exp[r][cid] += *c_it;
            --zero_cnt[r][cid];
        }
    }


    for (int i = 0; i < n_genes; ++i) {
        res[i].gene_id = i + 1;

        double m1 = total_exp[i][0] / total_cnt[0];
        double m2 = total_exp[i][1] / total_cnt[1];

        res[i].log_fc = log((m1 / m2) + 1);

        ProcessGene(
            std::move(exp[i]),
            total_cnt,
            std::move(zero_cnt[i]),
            thres,
            perm,
            res[i]
        );

        if ((i + 1) % 1000 == 0)
            Rcout << "Processed " << i + 1 << " genes\r";
    }

    return res;
}

std::vector<struct GeneResult> HarmonyTest(
        com::bioturing::Hdf5Util &oHdf5Util,
        HighFive::File *file,
        const Rcpp::NumericVector &cluster,
        const std::array<int, 2> &total_cnt,
        int threshold)
{
    int thres = threshold == 0? GetThreshold(total_cnt) : threshold;

    if (thres < MINIMAL_SAMPLE)
        throw std::runtime_error("Threshold is too small."
            "Maybe the number of cells in one cluster is too small");

    std::vector<int> shape;
    oHdf5Util.ReadDatasetVector<int>(file, GROUP_NAME, "shape", shape);

    int n_genes = shape[1];

    if (cluster.size() != shape[0])
        throw std::domain_error("Input cluster size is not equal to "
                                "the number of columns in matrix");

    std::vector<struct GeneResult> res(n_genes);

    for (int i = 0; i < n_genes; ++i) {
        std::vector<int> col_idx;
        std::vector<double> g_exp;
        oHdf5Util.ReadGeneExpH5(file, GROUP_NAME, i, col_idx,  g_exp);

        std::vector<std::pair<double, int>> exp(col_idx.size());
        std::array<int, 2> zero_cnt = {total_cnt[0], total_cnt[1]};

        for (int k = 0; k < col_idx.size(); ++k) {
            int idx = (int)cluster[col_idx[k]];
            if (idx) {
                exp[k] = {g_exp[k], idx - 1};
                --zero_cnt[idx - 1];
            }
        }

        ProcessGene(
            std::move(exp),
            total_cnt,
            zero_cnt,
            thres,
            0,
            res[i]
        );
    }

    return res;
}

DataFrame PostProcess(
        std::vector<struct GeneResult> &res,
        std::vector<std::string> &rownames)
{
    int n_gene = res.size();
    std::vector<std::pair<double,int>> order(n_gene);

    for (int i = 0; i < n_gene; ++i)
        order[i] = std::make_pair(res[i].log_p_value, i);

    std::sort(order.begin(), order.end());

    std::vector<std::string> g_names(n_gene);
    std::vector<int> g_id(n_gene);

    std::vector<double> d_score(n_gene), ud_score(n_gene), log2_fc(n_gene);
    std::vector<double> log10_pv(n_gene), perm_pv(n_gene), log10_adj_pv(n_gene);
    std::vector<double> b_cnt(n_gene);

    //Adjust p value
    double prev = -std::numeric_limits<double>::infinity();
    for(int i = 0; i < n_gene; ++i) {
        double log_p = order[i].first  + log(n_gene) - log(i + 1);

        if (log_p > 0)
            log_p = 0;

        if (log_p > prev)
            prev = log_p;

        log10_adj_pv[i] = prev * M_LOG10E;
    }

    for(int i = 0; i < n_gene; ++i) {
        int k = order[i].second;

        g_names[i]  = rownames[k];
        g_id[i]     = res[k].gene_id;
        d_score[i]  = res[k].d_score;
        b_cnt[i]    = res[k].b_cnt;
        log10_pv[i] = res[k].log_p_value * M_LOG10E;
        perm_pv[i]  = res[k].perm_p_value;
        ud_score[i] = res[k].ud_score;
        log2_fc[i]  = res[k].log_fc * M_LOG2E;
    }

    Rcout << "Done all" << std::endl;
    return DataFrame::create(
                    Named("Gene ID")                = wrap(g_id),
                    Named("Gene Name")              = wrap(g_names),
                    Named("Dissimilarity")          = wrap(d_score),
                    Named("Bin count")              = wrap(b_cnt),
                    Named("Log10 p value")          = wrap(log10_pv),
                    Named("Perm p value")           = wrap(perm_pv),
                    Named("Log10 adjusted p value") = wrap(log10_adj_pv),
                    Named("Up-Down score")          = wrap(ud_score),
                    Named("Log2 fold change")       = wrap(log2_fc)
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
DataFrame HarmonyMarker(
        const Rcpp::S4 &S4_mtx,
        const Rcpp::NumericVector &cluster,
        int threshold = 0,
        int perm = 0)
{
    Rcout << "Enter" << std::endl;

    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    const arma::sp_mat &mtx = Rcpp::as<arma::sp_mat>(S4_mtx);
    Rcout << "Done parse" << std::endl;

    std::vector<struct GeneResult> res
            = HarmonyTest(mtx, cluster, total_cnt, threshold, perm);
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
    const Rcpp::NumericVector &cluster, int threshold = 0)
{
    com::bioturing::Hdf5Util oHdf5Util(hdf5Path);
    HighFive::File *file = oHdf5Util.Open(1);

    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    Rcout << "Group1 " << total_cnt[0]
          << "Group2 " << total_cnt[1] << std::endl;

    std::vector<struct GeneResult> res
        = HarmonyTest(oHdf5Util, file, cluster, total_cnt, threshold);

    Rcout << "Done calculate" << std::endl;
    std::vector<std::string> rownames;
    // Read the barcode slot since this is the transposed matrix
    oHdf5Util.ReadDatasetVector<std::string>(file, GROUP_NAME,
                                            "barcodes", rownames);
    oHdf5Util.Close(file);

    return PostProcess(res, rownames);
}
