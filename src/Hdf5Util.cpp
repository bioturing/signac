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
