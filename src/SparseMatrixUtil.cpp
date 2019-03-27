#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <hdf5.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

// [[Rcpp::export]]
arma::sp_mat FastCreateSparseMat(int nrow, int ncol) {
    return arma::speye(nrow, ncol);
}

// [[Rcpp::export]]
Rcpp::List FastStatsOfSparseMat(const arma::sp_mat &mat) {
    return Rcpp::List::create(mat.n_rows, mat.n_cols, mat.n_elem, mat.n_nonzero);
}

// [[Rcpp::export]]
arma::sp_mat FastCreateFromTriplet(arma::urowvec vec1, arma::urowvec vec2, arma::colvec vec_val) {
    arma::umat loc = arma::join_vert(vec1, vec2);
    arma::sp_mat sp(loc, vec_val);
    return sp;
}

// [[Rcpp::export]]
arma::sp_mat FastConvertToSparseMat(const SEXP &s) {
    return Rcpp::as<arma::sp_mat>(s);
}

// [[Rcpp::export]]
Rcpp::List FastConvertToTripletMat(const SEXP &s) {
    return Rcpp::simple_triplet_matrix(Rcpp::as<arma::sp_mat>(s));
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatSqrt(const arma::sp_mat &mat) {
    return arma::sqrt(mat);
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatMult(const arma::sp_mat &mat1, const arma::sp_mat &mat2) {
    return mat1 * mat2;
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatAddition(const arma::sp_mat &mat1, const arma::sp_mat &mat2) {
    return mat1 + mat2;
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatMultWithNum(const arma::sp_mat &mat, const int &k) {
    return k * mat;
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatSymmatl(const arma::sp_mat &mat) {
    return arma::symmatl(mat);
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatTranspose(const arma::sp_mat &mat) {
    return arma::trans(mat);
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatTrimatu(const arma::sp_mat &mat) {
    return arma::trimatu(mat);
}

// [[Rcpp::export]]
int FastSparseMatTrace(const arma::sp_mat &mat) {
    return arma::trace(mat);
}

// [[Rcpp::export]]
arma::sp_mat FastConvertToDiagonalSparseMat(arma::sp_mat &mat) {
    mat.diag().ones();
    return mat;
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatSquare(const arma::sp_mat &mat) {
    return arma::square(mat);
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatRepmat(const arma::sp_mat &mat, const int &i, const int &j) {
    return arma::repmat(mat, i, j);
}

// [[Rcpp::export]]
arma::sp_mat FastSparseMatSign(const arma::sp_mat &mat) {
    return arma::sign(mat);
}
