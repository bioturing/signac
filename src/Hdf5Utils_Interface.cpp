#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG
#define ARMA_USE_HDF5

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]
// [[Rcpp::depends(BH)]]
#include "Hdf5Util.h"
#include "Hdf5Utils_Interface.h"
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

//' ReadH5Vector
//'
//' This function is used to read a vector from hdf5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
std::vector<int> ReadH5Vector(const std::string &filePath,
                                 const std::string &groupName,
                                 const std::string &datasetName) {
    std::vector<int> vec;
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    oHdf5Util.ReadDatasetVector(groupName, datasetName, vec);
    return vec;
}

//' ReadH5VectorRange
//'
//' This function is used to read a vector[i:j] from hdf5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
void ReadGeneExpH5(const std::string &filePath,
                   const std::string &groupName,
                   const int g_idx,
                   std::vector<int> &col_idx,
                   std::vector<double> &g_exp) {
    std::vector<int> vec = ReadH5Vector(filePath, groupName, "indptr");
    ReadH5VectorRange(filePath, groupName, "indices", (unsigned int)vec[g_idx],
                      (unsigned int)vec[g_idx + 1], col_idx);
    ReadH5VectorRange(filePath, groupName, "data", (unsigned int)vec[g_idx],
                      (unsigned int)vec[g_idx + 1], g_exp);
}


