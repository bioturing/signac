#ifndef HDF5_UTIL_INTERFACE
#define HDF5_UTIL_INTERFACE

#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]
// [[Rcpp::depends(BH)]]
#include "Hdf5Util.h"

/*
 * @abstract: Simple function to read a numeric vector from HDF5
 * @param filePath: path to the hdf5 file
 * @param groupName: group name e.g "matrix" or "bioturing"
 * @param datasetName: dataset name e.g "shape" or "indices" or "data"
 * @example:
 * std::vector<int> shape = ReadH5Vector("inst/v3input.h5", "matrix", "shape");
 */
Rcpp::NumericVector ReadH5Vector(const std::string &filePath,
                                 const std::string &groupName,
                                 const std::string &datasetName);

/*
 * @abstract: Read a numeric vector[i:j) from HDF5. This function will not be
 * exported;
 * @param filePath: path to the hdf5 file
 * @param groupName: group name e.g "matrix" or "bioturing"
 * @param datasetName: dataset name e.g "shape" or "indices" or "data"
 * @param start: the first index of the vector that will be read
 * @param end: the first index of the vector that will not be read
 * @example
 * std::vector<int> col_indices;
 * ReadH5VectorRange("inst/v3input.h5", "matrix", 0, 100, col_indices);
 * std::vector<double> gene_exp;
 * ReadH5VectorRange("inst/v3input.h5", "matrix", 0, 100, gene_exp);
 */
template<typename T>
void ReadH5VectorRange(const std::string &filePath,
                                 const std::string &groupName,
                                 const std::string &datasetName,
                                 const int start,
                                 const int end,
                                 std::vector<T> &vec);

/*
 * @abstract: Read the expression vector of one gene
 * @param filePath: path to the hdf5 file
 * @param groupName: group name e.g "matrix" or "bioturing"
 * @param g_idx: index of the gene (0-based)
 * @param col_idx: indices of columns that have the expression of the gene > 0
 * @param g_exp: vector of expression value of the gene. (|g_exp| == |col_idx|)
 * @example
 * std::vector<int> col_idx;
 * std::vector<double> g_exp;
 * ReadGeneExpH5("inst/v3input.h5", "matrix", 1, col_idx, g_exp);
 */
void ReadGeneExpH5(const std::string &filePath,
                   const std::string &groupName,
                   const int g_idx,
                   std::vector<int> &col_idx,
                   std::vector<double> &g_exp);
#endif //HDF5_UTIL_INTERFACE
