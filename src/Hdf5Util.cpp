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

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

/*
h5::fd_t FastOpenH5File(const std::string &file_path) {
    h5::fd_t file;
    try {
        file = h5::create(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    catch (std::exception& err) {
#ifdef DEBUG
        std::cerr << "Open file with error=" << err.what() << std::endl;
#endif
    }
    return file;
}
*/

