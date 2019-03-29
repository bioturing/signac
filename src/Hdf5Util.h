#ifndef HDF5_UTIL
#define HDF5_UTIL

#include <RcppArmadillo.h>
#include <unordered_map>
#include <fstream>
#include <string>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

using namespace Rcpp;

void FastCreateH5File(const std::string &file_path);

#endif //HDF5_UTIL
