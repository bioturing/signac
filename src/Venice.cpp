#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#define VENICE_EPS 1e-50
#define MINIMAL_SAMPLE 10
#define MINIMAL_BIN 2
#define GROUPING_RATE 0.6

#define C_INSIDE 1
#define C_OUTSIDE 0

#include <boost/sort/spreadsort/spreadsort.hpp>
#include <Rmath.h>

#include <progress.hpp>
#include <progress_bar.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <H5Cpp.h>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

// [[Rcpp::depends(Rhdf5lib)]]
// [[Rcpp::depends(RcppProgress)]]

struct GeneResult {
    int gene_id;

    double d_score; //dissimilarity score
    double d_bias; //correction for d_score bias
    int b_cnt;

    double log_p_value;
    double perm_p_value;

    double ud_score; //up-down score

    double log_fc;

    double pct1; // percent 
    double pct2;
};

class MatrixReader {
private:
    HighFive::File file;
    std::string group_name;
    std::array<int, 2> shape;
    std::vector<uint64_t> ptr;
    HighFive::DataSet indices;
    HighFive::DataSet data;
    
    uint64_t cache_size;

    uint64_t start = 0;
    uint64_t end = 0;
    std::vector<int> indices_v;
    std::vector<double> data_v;

    void read(uint64_t start) {
        assert(start < ptr.size());
        uint64_t l = this->ptr[start];
        uint64_t end = start + 1;
        while (end < this->ptr.size() - 1 && this->ptr[end + 1] - l <= this->cache_size)
            ++end;

        uint64_t r = this->ptr[end];
        uint64_t len = r - l;

        if (len > this->indices_v.size()) {
            this->indices_v.resize(len);
            this->data_v.resize(len);
        }

        this->start = start;
        this->end = end;
        this->indices.select({l}, {len}).read(this->indices_v);
        this->data.select({l}, {len}).read(this->data_v);
    }
public:
    MatrixReader(const std::string &file_name, const std::string &group_name, uint64_t cache_size = 1000000): 
            file(file_name, HighFive::File::ReadOnly),
            group_name(group_name),
            indices {this->file.getDataSet(group_name + "/indices")},
            data {this->file.getDataSet(group_name + "/data")},
            cache_size(cache_size) {
        if(!this->file.exist(group_name))
            Rcpp::stop("Group " + group_name + "doesn't exist in " + this->file.getName());

        // Read index pointer
        if (!this->file.exist(group_name + "/indptr") ||
            !this->file.exist(group_name + "/barcodes") ||
            !this->file.exist(group_name + "/shape"))  {
                
            Rcpp::stop("Group " + group_name + " in file " + this->file.getName() + " doesn't have a correct matrix format");
        }

        this->file.getDataSet(group_name + "/indptr").read(this->ptr);
        this->file.getDataSet(group_name + "/shape").read(this->shape);

        if (cache_size < shape[1])
            Rcpp::warning("Cache size maybe too small for this data, venice will use more memory when needed");

        this->indices_v.resize(cache_size);
        this->data_v.resize(cache_size);
    }


    void get_expression(int gene_id, const Rcpp::IntegerVector &cluster, 
                        std::vector<std::pair<double, bool>> &result) {
        
        assert(gene_id >= start);
        if (gene_id >= end)
            read(gene_id);

        uint64_t l = this->ptr[gene_id] - this->ptr[start];
        uint64_t r = this->ptr[gene_id + 1] - this->ptr[start];
        uint64_t len = r - l;

        result.clear();
        result.reserve(len);
        for (uint64_t i = l; i < r; ++i) {
            int c_id = -1;
            switch (cluster[this->indices_v[i]]) {
                case C_INSIDE:
                    c_id = 0;
                    break;
                case C_OUTSIDE:
                    c_id = 1;
                    break;
                default:
                    continue;
            }
            result.push_back({this->data_v[i], c_id});
        }
    }

    std::vector<std::string> get_names() {
        std::vector<std::string> gene;
        HighFive::DataSet datasetVec = this->file.getDataSet(group_name + "/barcodes");
        HighFive::DataSpace dataSpace = datasetVec.getSpace();
        std::vector<size_t> dims = dataSpace.getDimensions();
        HighFive::DataType dataType = datasetVec.getDataType();
        size_t str_size = H5Tget_size(dataType.getId());
        H5T_str_t str_pad = H5Tget_strpad(dataType.getId());
        H5T_cset_t str_cset = H5Tget_cset(dataType.getId());

        if((str_pad == H5T_STR_NULLTERM) && (str_cset == H5T_CSET_UTF8))
            datasetVec.read(gene);
        else
            datasetVec.read(gene, str_size, str_pad, str_cset, dims[0]);

        return gene;
    }

    std::array<int, 2> get_shape() {
        return this->shape;
    }
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
    double x = 2 * score * (HarmonicMean(n1, n2) - 1);

    return R::pchisq(x, Dof, false, true);
}

void GetTotalCount(
    const Rcpp::IntegerVector &cluster,
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
            const std::vector<std::pair<double, bool>> & exp,
            const std::array<int, 2> &cnt)
{
    double m[2];

    for (int i = 0; i < exp.size(); ++i)
        m[exp[i].second] += exp[i].first;
    
    m[0] /= cnt[0];
    m[1] /= cnt[1];

    return log(m[0] + 1) - log(m[1] + 1);
}

struct rightshift {
  inline long long operator()(const std::pair<double, int> &x, const unsigned offset) const {
    return boost::sort::spreadsort::float_mem_cast<double, long long>(x.first) >> offset;
  }
};

std::array<int, 2> count_missing(const std::vector<std::pair<double, bool>> &exp, const std::array<int, 2> &total_cnt) {
    std::array<int, 2> zero_cnt = total_cnt;
    for (const std::pair<double, bool> &e : exp)
        --zero_cnt[e.second];
    return zero_cnt;
}

std::vector<std::array<int, 2>> Binning(
        std::vector<std::pair<double, bool>> exp,
        const std::array<int, 2> &total_cnt)
{
    std::vector<std::array<int, 2>> result(exp.size() + 1);    //+1 for zero

    if (exp.size() == 0) {
        result[0] = total_cnt;
        return result;
    }

    std::array<int, 2> zero_cnt = count_missing(exp, total_cnt);
    boost::sort::spreadsort::float_sort(exp.begin(), exp.end(), rightshift());

    double p_exp = -std::numeric_limits<double>::infinity();

    int i = 0, j = -1;
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

        if (std::abs(c_exp - p_exp) >= VENICE_EPS) {
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

int CountNumExpressClusterZero(const std::vector<std::pair<double, bool>> &exp) {
  int count = 0;
    for (int i = 0; i < exp.size(); ++i)
        if (exp[i].second)
            ++count;
  return count;
}

bool IsValidPercent(const GeneResult &gene_result, double threshold_pct) {
  if (gene_result.pct1 < threshold_pct && gene_result.pct2 < threshold_pct) 
      return false;
  return true;
}

struct GeneResult ProcessGene(
        std::vector<std::pair<double, bool>> exp,
        const std::array<int, 2> &cnt,
        int thres,
        int perm,
        double threshold_pct)
{
    struct GeneResult res;
    
    int count_cluster0 = CountNumExpressClusterZero(exp);
    res.pct1 = 100.0 * (exp.size() - count_cluster0) / cnt[0];
    res.pct2 = 100.0 * count_cluster0 / cnt[1];
    if (!IsValidPercent(res, threshold_pct))
      return res;
    
    res.log_fc = ComputeLogFC(exp, cnt);

    std::vector<std::array<int, 2>> bins = Binning(std::move(exp), cnt);

    res.ud_score = ComputeUd_score(bins, cnt);

    Grouping(bins, thres);

    std::pair<double, double> score = ComputeSimilarity(bins, cnt);
    res.d_score = score.first;
    res.d_bias = score.second;

    res.b_cnt = bins.size();
    res.log_p_value = LogPvalue(res.d_score, cnt[0], cnt[1], bins.size());
    res.perm_p_value = PermPvalue(std::move(bins), cnt, res.d_score, perm);

    return res;
}

void BuildExpDgT(
        const Rcpp::S4 &mtx, 
        const Rcpp::IntegerVector &cluster,
        const std::array<int, 2> &total_cnt,
        std::vector<std::vector<std::pair<double, bool>>> &exp)
{
    if (!mtx.is("dgTMatrix"))
        throw std::runtime_error("wrong format (should be dgTMatrix)");

    const Rcpp::IntegerVector &dim = mtx.attr("Dim");

    const int n_genes = dim[0];
    const int n_cells = dim[1];

    if (cluster.size() != n_cells)
        throw std::domain_error("The length of cluster vector (" + std::to_string(cluster.size()) + ") is not equal "
                                "to the number of cells (" + std::to_string(n_cells) + ") in the matrix");

    exp.resize(n_genes);    

    const Rcpp::IntegerVector &row_id = mtx.attr("i");
    const Rcpp::IntegerVector &col_id = mtx.attr("j");
    const Rcpp::NumericVector &raw_exp = mtx.attr("x");

    uint64_t n_exp = raw_exp.size();

    for (uint64_t i = 0; i < n_exp; ++i) {
        int col = col_id[i];

        bool cid;

        switch (cluster[col]) {
            case C_INSIDE:
                cid = 0;
                break;
            case C_OUTSIDE:
                cid = 1;
                break;
            default:
                continue;
        }

        int row = row_id[i];
        double x = raw_exp[i];

        exp[row].push_back({x, cid});
    }
}

void BuildExpDgC(
        const Rcpp::S4 &mtx, 
        const Rcpp::IntegerVector &cluster,
        const std::array<int, 2> &total_cnt,
        std::vector<std::vector<std::pair<double, bool>>> &exp)
{
    if (!mtx.is("dgCMatrix"))
        throw std::runtime_error("wrong format (should be dgCMatrix)");

    const Rcpp::IntegerVector &dim = mtx.attr("Dim");

    const int n_genes = dim[0];
    const int n_cells = dim[1];

    if (cluster.size() != n_cells)
        throw std::domain_error("The length of cluster vector (" + std::to_string(cluster.size()) + ") is not equal "
                                "to the number of cells (" + std::to_string(n_cells) + ") in the matrix");

    exp.resize(n_genes);

    const Rcpp::IntegerVector &row_id = mtx.attr("i");
    const Rcpp::IntegerVector &pt = mtx.attr("p");
    const Rcpp::NumericVector &raw_exp = mtx.attr("x");


    for (int col = 0; col < n_cells; ++col) {
        bool cid;

        switch (cluster[col]) {
            case C_INSIDE:
                cid = 0;
                break;
            case C_OUTSIDE:
                cid = 1;
                break;
            default:
                continue;
        }

        for (uint64_t i = pt[col]; i < pt[col + 1]; ++i) {
            int row = row_id[i];
            double x = raw_exp[i];
            exp[row].push_back({x, cid});
        }
    }
}


std::vector<struct GeneResult> VeniceTest(
        const Rcpp::S4 &mtx,
        const Rcpp::IntegerVector &cluster,
        int threshold,
        double threshold_pct,
        int perm,
        bool display_progress = true)
{
    if (!mtx.is("sparseMatrix"))
        throw std::runtime_error("The input expression matrix "
                                 "must be a sparseMatrix");

    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    int thres = threshold == 0? GetThreshold(total_cnt) : threshold;

    std::vector<std::vector<std::pair<double, bool>>> exp;

    const int n_genes = ((Rcpp::IntegerVector)mtx.attr("Dim"))[0];

    Progress p(n_genes, display_progress);

    if(mtx.is("dgTMatrix"))
        BuildExpDgT(mtx, cluster, total_cnt, exp);
    else if (mtx.is("dgCMatrix"))
        BuildExpDgC(mtx, cluster, total_cnt, exp);
    else
        throw std::runtime_error("this matrix format is not support. Please convert to dgTMatrix/dgCMatrix");

    std::vector<struct GeneResult> res;
    res.reserve(n_genes);

    for (int i = 0; i < n_genes; ++i) {

        if (Progress::check_abort()) {
            res.resize(i);
            return res;
        }

        GeneResult gene_result = ProcessGene(
                                  std::move(exp[i]),
                                  total_cnt,
                                  thres,
                                  perm,
                                  threshold_pct);
        
        if (IsValidPercent(gene_result, threshold_pct)) {
            res.push_back(gene_result);
            res.back().gene_id = i + 1;
        }
        
        p.increment();   
    }
    return res;
}

// Venice test for transposed matrix
std::vector<struct GeneResult> VeniceTest(
        MatrixReader &matrix,
        const Rcpp::IntegerVector &cluster,
        int threshold,
        double threshold_pct,
        int perm,
        bool correct,
        bool display_progress = true)
{
    std::array<int, 2> total_cnt;
    GetTotalCount(cluster, total_cnt);

    int thres = threshold == 0? GetThreshold(total_cnt) : threshold;

    if (thres < MINIMAL_SAMPLE)
        throw std::runtime_error("Threshold is too small."
            "Maybe the number of cells in one cluster is too small");

    std::array<int, 2> shape = matrix.get_shape();

    int n_genes = shape[1];
    if (cluster.size() != shape[0])
        throw std::domain_error("Input cluster size is not equal to "
                                "the number of columns in matrix");

    Progress p(n_genes, display_progress);

    std::vector<struct GeneResult> res;
    res.reserve(n_genes);

    for (int i = 0; i < n_genes; ++i) {
        if (Progress::check_abort()) {
            res.resize(i);
            return res;
        }

        std::vector<std::pair<double,bool>> exp;
        matrix.get_expression(i, cluster, exp);

        GeneResult gene_result = ProcessGene(
                                  std::move(exp),
                                  total_cnt,
                                  thres,
                                  perm,
                                  threshold_pct);
        
        if (IsValidPercent(gene_result, threshold_pct)) {
            res.push_back(gene_result);
            res.back().gene_id = i + 1;    
        }
        
        p.increment();
    }

    return res;
}

Rcpp::DataFrame PostProcess(
        const std::vector<struct GeneResult> &res,
        const std::vector<std::string> &gene_names,
        bool correct = true)
{
    int n_gene = res.size();
    std::vector<std::pair<double,int>> order(n_gene);

    for (int i = 0; i < n_gene; ++i)
        order[i] = std::make_pair(res[i].log_p_value, i);

    std::sort(order.begin(), order.end());

    std::vector<std::string> g_names(n_gene);
    std::vector<int> g_id(n_gene);

    std::vector<double> d_score(n_gene), ud_score(n_gene), log2_fc(n_gene), pct1(n_gene), pct2(n_gene);
    std::vector<double> log10_pv(n_gene), perm_pv(n_gene), log10_adj_pv(n_gene);
    std::vector<int> b_cnt(n_gene);

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

        g_names[i]  = gene_names[res[k].gene_id - 1];
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
        pct1[i]     = res[k].pct1;
        pct2[i]     = res[k].pct2;
    }
#ifdef DEBUG
    Rcout << "Done!" << std::endl;
#endif
    return Rcpp::DataFrame::create(
                    Rcpp::Named("Gene ID")                = Rcpp::wrap(g_id),
                    Rcpp::Named("Gene Name")              = Rcpp::wrap(g_names),
                    Rcpp::Named("Dissimilarity")          = Rcpp::wrap(d_score),
                    Rcpp::Named("Bin count")              = Rcpp::wrap(b_cnt),
                    Rcpp::Named("Log10 p value")          = Rcpp::wrap(log10_pv),
                    Rcpp::Named("Perm p value")           = Rcpp::wrap(perm_pv),
                    Rcpp::Named("Log10 adjusted p value") = Rcpp::wrap(log10_adj_pv),
                    Rcpp::Named("Up-Down score")          = Rcpp::wrap(ud_score),
                    Rcpp::Named("Log2 fold change")       = Rcpp::wrap(log2_fc),
                    Rcpp::Named("pct1")                   = Rcpp::wrap(pct1),
                    Rcpp::Named("pct2")                   = Rcpp::wrap(pct2)
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
Rcpp::DataFrame VeniceMarker(
        const Rcpp::S4 &S4_mtx,
        const Rcpp::IntegerVector &cluster,
        int threshold = 0,
        double threshold_pct = 0,
        int perm = 0,
        bool correct = true,
        bool verbose = false)
{

    std::vector<struct GeneResult> res = VeniceTest(S4_mtx, cluster, threshold, threshold_pct, perm, verbose);
#ifdef DEBUG
    Rcout << "Done calculate" << std::endl;
#endif
    Rcpp::List dim_names = Rcpp::List(S4_mtx.attr("Dimnames"));
    std::vector<std::string> rownames = dim_names[0];
    return PostProcess(res, rownames, correct);
}

//' VeniceMarkerTransposedH5
//'
//' Find gene marker for a cluster in tranposed H5 file
//'
//' @param hdf5Path A string path
//' @param cluster A numeric vector
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame VeniceMarkerTransposedH5(
    const std::string &hdf5Path,
    const std::string &group_name,
    const Rcpp::IntegerVector &cluster,
    int threshold = 0,
    double threshold_pct = 0,
    int perm = 0,
    bool correct = true,
    bool verbose = false)
{
    MatrixReader matrix(hdf5Path, group_name);

    std::vector<struct GeneResult> res
        = VeniceTest(matrix, cluster, threshold, threshold_pct, perm, correct, verbose);
#ifdef DEBUG
    Rcout << "Done calculate" << std::endl;
#endif
    return PostProcess(res, matrix.get_names(), correct);
}
