#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#define VENICE_EPS 1e-50
#define MINIMAL_SAMPLE 10
#define MINIMAL_BIN 2
#define GROUPING_RATE 0.6
#define GROUP_NAME "bioturing"

#define C_INSIDE 1
#define C_OUTSIDE 0

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <Rmath.h>

#include <progress.hpp>
#include <progress_bar.hpp>

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


using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace std::placeholders;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]
// [[Rcpp::depends(RcppProgress)]]


struct GeneResult {
    int gene_id;
    std::string gene_name;

    double d_score; //dissimilarity score
    double d_bias; //correction for d_score bias
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

inline std::pair<double, double> Score(int x, int y, int n1, int n2)
{
    if (x == 0 && y == 0)
        return std::make_pair(0.0, 0.0);

    double a = (double)x / n1;
    double b = (double)y / n2;

    double result = pow(a - b, 2) / (2 * (a + b));
    //bias correction
    double correction = 2 * a * b * (x * (1 - b) + y * (1 - a))
                        / (n1 * n2 * pow(a + b, 3));

    return std::make_pair(result, correction);
}

inline double LogPvalue(double score, int n1, int n2, int b_cnt)
{
    if (b_cnt <= 0)
        throw std::domain_error("Bin count should be positive");

    int Dof = b_cnt - 1;
    double x = 2 * score * HarmonicMean(n1, n2);

    return R::pchisq(x, Dof, false, true);
}

void GetTotalCount(
    const Rcpp::NumericVector &cluster,
    std::array<int, 2> &total_cnt)
{
    total_cnt[0] = total_cnt[1] = 0;

    for (int i = 0; i < cluster.size(); ++i) {
        int c = (int)cluster[i];
        if (c == C_INSIDE)
            ++total_cnt[0];
        else if (c == C_OUTSIDE)
            ++total_cnt[1];
    }
}

int GetThreshold(const std::array<int,2> &cnt)
{
    double avg = HarmonicMean(cnt[0],cnt[1]);
    double est_bin = std::pow(avg,1 - GROUPING_RATE);

    int thres = std::max<int>(
                        MINIMAL_SAMPLE,
                        ((cnt[0] + cnt[1]) / est_bin)
                );
#ifdef DEBUG
    Rcout << cnt[0] << " " << cnt[1] << " " << avg << " " << thres << std::endl;
#endif

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
    const std::array<int, 2> &cnt)
{
    int j = 0;

    int n1 = cnt[0];
    int n2 = cnt[1];

    for (int i = 0; i < bins.size(); ++i) {

        int k = bins[i][0] + bins[i][1];

        int r = R::rhyper(n1, n2, k);

        bins[i][0] = r;
        bins[i][1] = k - r;

        n1 -= r;
        n2 -= k - r;
    }
}

// Bins is the group count for each UMI value after sorted
std::pair<double, double> ComputeSimilarity(
        const std::vector<std::array<int, 2>> &bins,
        const std::array<int, 2> &cnt)
{
    double d_score = 0;
    double d_bias = 0;

    int n1 = cnt[0];
    int n2 = cnt[1];

    int n = bins.size();

    int i = 0;

    for (; i < n; ++i) {
        int b1 = bins[i][0];
        int b2 = bins[i][1];

        std::pair<double, double> s = Score(b1, b2, n1, n2);

        d_score += s.first;
        d_bias += s.second;
    }

    return std::make_pair(d_score, d_bias);
}

double ComputeLogFC(
            const std::vector<std::pair<double, int>> & exp,
            const std::array<int, 2> &cnt)
{
    double m[2];

    for (int i = 0; i < exp.size(); ++i)
        m[exp[i].second] += exp[i].first;
    
    m[0] /= cnt[0];
    m[1] /= cnt[1];

    return log(m[0] + 1) - log(m[1] + 1);
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

        if (abs(c_exp - p_exp) >= VENICE_EPS) {
            ++j;
            p_exp = c_exp;
        }
        ++result[j][exp[i].second];
    }

    if (zero_cnt[0] + zero_cnt[1] > 0) {
        if (p_exp < -VENICE_EPS) {
            ++j;
            p_exp = 0;
        }

        result[j][0] += zero_cnt[0];
        result[j][1] += zero_cnt[1];
    }

    for (; i < exp.size(); ++i) {
        double c_exp = exp[i].first;

        if (abs(c_exp - p_exp) >= VENICE_EPS) {
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

double PermPvalue(
            std::vector<std::array<int, 2>> bins,
            const std::array<int, 2> &cnt,
            double d_score,
            int perm)
{
    if (bins.size() <= 1)
        return 1;

    int count = 0;
    for(int i = 0; i < perm; ++i) {
        Resample(bins, cnt);
        double score = ComputeSimilarity(bins, cnt).first;
        count += score >= d_score;
    }

    return (double)count/perm;
}

void ProcessGene(
        std::vector<std::pair<double, int>> exp,
        const std::array<int, 2> &cnt,
        const std::array<int, 2> &zero_cnt,
        int thres,
        int perm,
        struct GeneResult &res)
{
    res.log_fc = ComputeLogFC(exp, cnt);

    std::vector<std::array<int, 2>> bins = Binning(std::move(exp), zero_cnt);

    res.ud_score = ComputeUd_score(bins, cnt);

    Grouping(bins, thres);

    std::pair<double, double> score = ComputeSimilarity(bins, cnt);
    res.d_score = score.first;
    res.d_bias = score.second;

    res.b_cnt = bins.size();
    res.log_p_value = LogPvalue(res.d_score, cnt[0], cnt[1], bins.size());
    res.perm_p_value = PermPvalue(std::move(bins), cnt, res.d_score, perm);
}

std::vector<struct GeneResult> VeniceTest(
        const arma::sp_mat &mtx,
        const Rcpp::NumericVector &cluster,
        int threshold,
        int perm,
        bool display_progress = true)
{
    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    int thres = threshold == 0? GetThreshold(total_cnt) : threshold;
    int n_genes = mtx.n_rows;
    Progress p(n_genes, display_progress);

    if (cluster.size() != mtx.n_cols)
        throw std::domain_error("The length of cluster vector is not equal "
                                "to the number of columns in the matrix");

    std::vector<std::vector<std::pair<double, int>>> exp(n_genes);

    std::vector<struct GeneResult> res(n_genes);
    arma::sp_mat::const_col_iterator c_it;

    std::vector<std::array<int, 2>> zero_cnt(n_genes, total_cnt);

    for (int i = 0; i < mtx.n_cols; ++i) {
        if (Progress::check_abort())
            return {};

        int c = (int)cluster[i];
        int cid;

        if (c == C_INSIDE)
            cid = 0;
        else if (c == C_OUTSIDE)
            cid = 1;
        else
            continue;

        for (c_it = mtx.begin_col(i); c_it != mtx.end_col(i); ++c_it) {
            int r = c_it.row();

            exp[r].push_back({*c_it, cid});
            --zero_cnt[r][cid];
        }
    }


    for (int i = 0; i < n_genes; ++i) {
        res[i].gene_id = i + 1;

        if (Progress::check_abort())
            return {};

        ProcessGene(
            std::move(exp[i]),
            total_cnt,
            std::move(zero_cnt[i]),
            thres,
            perm,
            res[i]
        );

        p.increment();
    }

    return res;
}

std::vector<struct GeneResult> VeniceTest(
        com::bioturing::Hdf5Util &oHdf5Util,
        HighFive::File *file,
        const Rcpp::NumericVector &cluster,
        int threshold,
        int perm)
{
    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

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
            int c = (int)cluster[col_idx[k]];
            int cid;

            if (c == C_INSIDE)
                cid = 0;
            else if (c == C_OUTSIDE)
                cid = 1;
            else
                continue;
             
            exp[k] = {g_exp[k], cid};
            --zero_cnt[cid];
        }

        ProcessGene(
            std::move(exp),
            total_cnt,
            zero_cnt,
            thres,
            perm,
            res[i]
        );
    }

    return res;
}

DataFrame PostProcess(
        std::vector<struct GeneResult> &res,
        std::vector<std::string> &rownames,
        bool correct = true)
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

        if (correct)
            d_score[i] = res[k].d_score - res[k].d_bias;
        else
            d_score[i] = res[k].d_score;
        
        b_cnt[i]    = res[k].b_cnt;
        log10_pv[i] = res[k].log_p_value * M_LOG10E;
        perm_pv[i]  = res[k].perm_p_value;
        ud_score[i] = res[k].ud_score;
        log2_fc[i]  = res[k].log_fc * M_LOG2E;
    }
#ifdef DEBUG
    Rcout << "Done!" << std::endl;
#endif
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

//' VeniceMarker
//'
//' Find gene marker for a cluster in sparse matrix
//'
//' @param S4_mtx A sparse matrix
//' @param cluster A numeric vector
//' @export
// [[Rcpp::export]]
DataFrame VeniceMarker(
        const Rcpp::S4 &S4_mtx,
        const Rcpp::NumericVector &cluster,
        int threshold = 0,
        int perm = 0,
        bool correct = true,
        bool verbose = false)
{
#ifdef DEBUG
    Rcout << "Enter" << std::endl;
#endif

    const arma::sp_mat &mtx = Rcpp::as<arma::sp_mat>(S4_mtx);
#ifdef DEBUG
    Rcout << "Done parse" << std::endl;
#endif
    std::vector<struct GeneResult> res = VeniceTest(mtx, cluster, threshold, perm, verbose);
#ifdef DEBUG
    Rcout << "Done calculate" << std::endl;
#endif
    Rcpp::List dim_names = Rcpp::List(S4_mtx.attr("Dimnames"));
    std::vector<std::string> rownames = dim_names[0];
    return PostProcess(res, rownames, correct);
}

//' VeniceMarkerH5
//'
//' Find gene marker for a cluster in H5 file
//'
//' @param hdf5Path A string path
//' @param cluster A numeric vector
//' @export
// [[Rcpp::export]]
DataFrame VeniceMarkerH5(
    const std::string &hdf5Path,
    const Rcpp::NumericVector &cluster,
    int threshold = 0,
    int perm = 0,
    bool correct = true)
{
    com::bioturing::Hdf5Util oHdf5Util(hdf5Path);
    HighFive::File *file = oHdf5Util.Open(1);

    std::vector<struct GeneResult> res
        = VeniceTest(oHdf5Util, file, cluster, threshold, perm);
#ifdef DEBUG
    Rcout << "Done calculate" << std::endl;
#endif
    std::vector<std::string> rownames;
    // Read the barcode slot since this is the transposed matrix
    oHdf5Util.ReadDatasetVector<std::string>(file, GROUP_NAME,
                                            "barcodes", rownames);
    oHdf5Util.Close(file);

    return PostProcess(res, rownames, correct);
}
