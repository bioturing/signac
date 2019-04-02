#include "Hdf5Util.h"
using namespace HighFive;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

// [[Rcpp::export]]
bool WriteSpMtFromArma(const std::string &filePath, const std::string &groupName, const arma::sp_mat &mat) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.WriteSpMtFromArma(mat, groupName);
}

// [[Rcpp::export]]
bool WriteSpMtFromS4(const std::string &filePath, const std::string &groupName, const Rcpp::S4 &mat) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.WriteSpMtFromS4(mat, groupName);
}

// [[Rcpp::export]]
arma::sp_mat ReadSpMt(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.ReadSpMtAsArma(groupName);
}

// [[Rcpp::export]]
arma::sp_mat ReadSpMtAsS4(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.ReadSpMtAsS4(groupName);
}
