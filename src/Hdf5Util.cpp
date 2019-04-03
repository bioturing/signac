#include "Hdf5Util.h"
using namespace HighFive;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

//' WriteSpMtV2
//'
//' This function is used to write a sparse matrix (as arma)
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
bool WriteSpMtV2(const std::string &filePath, const std::string &groupName, const arma::sp_mat &mat) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.WriteSpMtFromArma(mat, groupName);
}

//' WriteSpMtV1
//'
//' This function is used to write a sparse matrix (as S4)
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
bool WriteSpMtV1(const std::string &filePath, const std::string &groupName, const Rcpp::S4 &mat) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.WriteSpMtFromS4(mat, groupName);
}

//' ReadSpMtV2
//'
//' This function is used to read a sparse matrix from HDF5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
arma::sp_mat ReadSpMtV2(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.ReadSpMtAsArma(groupName);
}

//' ReadSpMtV1
//'
//' This function is used to read a sparse matrix from HDF5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::S4 ReadSpMtV1(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.ReadSpMtAsS4(groupName);
}

//' ReadRowSumSpMt
//'
//' Read rows sums
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ReadRowSumSpMt(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    std::vector<double> sumVec;

    int status = oHdf5Util.ReadVector(sumVec, groupName, "rowsums");
    if(status == -1) {
        std::stringstream ostr;
        ostr << "Can not open file :" << filePath;
        ::Rf_error(ostr.str().c_str());
        throw;
    } else if(status == -2) {
        std::stringstream ostr;
        ostr << "Can not exist dataset :" << groupName;
        ::Rf_error(ostr.str().c_str());
        throw;
    } else if(status == 1) {
        throw;
    }

    if(sumVec.size() == 0) {
        arma::sp_mat mat = oHdf5Util.ReadSpMtAsArma(groupName);
        if(mat.size() == 0) {
            Rcpp::S4 s = oHdf5Util.ReadSpMtAsS4(groupName);
            mat = com::bioturing::Hdf5Util::FastConvertS4ToSparseMT(s);
        }

        sumVec.resize(mat.n_rows);
        for (int i= 0; i< mat.n_cols; i++)
        {
            for (arma::sp_mat::const_col_iterator cij = mat.begin_col(i); cij != mat.end_col(i); ++cij) {
                sumVec[cij.row()] += (*cij);
            }
        }
        oHdf5Util.WriteVector(sumVec, groupName, "rowsums");
    }

    Rcpp::NumericVector result(sumVec.begin(), sumVec.end());
    return result;
}

//' ReadColSumSpMt
//'
//' Read cols sums
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ReadColSumSpMt(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    std::vector<double> sumVec;

    int status = oHdf5Util.ReadVector(sumVec, groupName, "colsums");
    if(status == -1) {
        std::stringstream ostr;
        ostr << "Can not open file :" << filePath;
        ::Rf_error(ostr.str().c_str());
        throw;
    } else if(status == -2) {
        std::stringstream ostr;
        ostr << "Can not exist dataset :" << groupName;
        ::Rf_error(ostr.str().c_str());
        throw;
    } else if(status == 1) {
        throw;
    }

    if(sumVec.size() == 0) {
        arma::sp_mat mat = oHdf5Util.ReadSpMtAsArma(groupName);
        if(mat.size() == 0) {
            Rcpp::S4 s = oHdf5Util.ReadSpMtAsS4(groupName);
            mat = com::bioturing::Hdf5Util::FastConvertS4ToSparseMT(s);
        }

        sumVec.resize(mat.n_cols);
        com::bioturing::SumColumWorker<std::vector<double>> sumColWorker(&mat, sumVec);
        RcppParallel::parallelFor(0, mat.n_cols, sumColWorker);
        oHdf5Util.WriteVector(sumVec, groupName, "colsums");
    }

    Rcpp::NumericVector result(sumVec.begin(), sumVec.end());
    return result;
}

