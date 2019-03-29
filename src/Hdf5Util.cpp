#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <hdf5.h>
#include "CommonUtil.h"
#include "highfive/H5File.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

// [[Rcpp::export]]
void FastCreateH5File(const std::string &file_path) {
    hid_t file_id;
    file_id = H5Fcreate(file_path.data(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose(file_id);
}
