#include "Hdf5Util.h"

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat FastMatMult(const arma::mat &mat1, const arma::mat &mat2) {
    return mat1 * mat2;
}

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

// [[Rcpp::export]]
arma::mat FastGetSubMat(const arma::mat &mat, arma::uvec rvec, arma::uvec cvec) {
    arma::uvec new_rvec(rvec.size());
    PerformRVector(rvec, (int)mat.n_rows, new_rvec);
    arma::uvec new_cvec(cvec.size());
    PerformRVector(cvec, (int)mat.n_cols, new_cvec);
    return mat.submat(new_rvec, new_cvec);
}
