#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG
#define ARMA_USE_HDF5

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]
// [[Rcpp::depends(BH)]]
#include "Hdf5Util.h"

//' FastMatMult
//'
//' This function is used to add two matrix
//'
//' @param mat1 A matrix
//' @param mat2 A matrix
//' @export
// [[Rcpp::export]]
arma::mat FastMatMult(const arma::mat &mat1, const arma::mat &mat2) {
    return mat1 * mat2;
}

//' FastGetRowsOfMat
//'
//' This function is used to get some rows of matrix
//'
//' @param mat A matrix
//' @param vec A row vector
//' @export
// [[Rcpp::export]]
arma::mat FastGetRowsOfMat(const arma::mat &mat, arma::uvec vec) {
    arma::uvec new_vec(vec.size());
    PerformRVector(vec, (int)mat.n_rows, new_vec);
    return mat.rows(new_vec);
}

// [[Rcpp::export]]
arma::mat FastGetColsOfMat(const arma::mat &mat, arma::uvec vec) {
    arma::uvec new_vec(vec.size());
    PerformRVector(vec, (int)mat.n_cols, new_vec);
    return mat.cols(new_vec);
}

//' FastGetSubMat
//'
//' This function is used to get subview matrix
//'
//' @param mat A matrix
//' @param vec A row vector
//' @param vec A col vector
//' @export
// [[Rcpp::export]]
arma::mat FastGetSubMat(const arma::mat &mat, arma::uvec rvec, arma::uvec cvec) {
    arma::uvec new_rvec(rvec.size());
    PerformRVector(rvec, (int)mat.n_rows, new_rvec);
    arma::uvec new_cvec(cvec.size());
    PerformRVector(cvec, (int)mat.n_cols, new_cvec);
    return mat.submat(new_rvec, new_cvec);
}
