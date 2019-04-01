#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include "CommonUtil.h"
#include <h5cpp/hdf5.hpp>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]



