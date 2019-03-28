#ifndef SPARSE_MATRIX_UTIL
#define SPARSE_MATRIX_UTIL

#include <RcppArmadillo.h>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;

arma::sp_mat FastConvertToSparseMat(const SEXP &s);
Rcpp::List FastConvertToTripletMat(const SEXP &s);
arma::sp_mat FastCreateSparseMat(int nrow, int ncol);
arma::sp_mat FastCreateFromTriplet(const arma::urowvec &vec1, const arma::urowvec &vec2, const arma::colvec &vec_val);
Rcpp::List FastStatsOfSparseMat(const arma::sp_mat &mat);
arma::sp_mat FastSparseMatTranspose(const arma::sp_mat &mat);
arma::sp_mat FastSparseMatSqrt(const arma::sp_mat &mat);
arma::sp_mat FastSparseMatMult(const arma::sp_mat &mat1, const arma::sp_mat &mat2);
arma::sp_mat FastSparseMatSymmatl(const arma::sp_mat &mat);
arma::sp_mat FastSparseMatAddition(const arma::sp_mat &mat1, const arma::sp_mat &mat2);
arma::sp_mat FastSparseMatMultWithNum(const arma::sp_mat &mat, const int &k);
arma::sp_mat FastSparseMatTrimatu(const arma::sp_mat &mat);
int FastSparseMatTrace(const arma::sp_mat &mat);
arma::sp_mat FastConvertToDiagonalSparseMat(arma::sp_mat &mat);
arma::sp_mat FastSparseMatSquare(const arma::sp_mat &mat);
arma::sp_mat FastSparseMatRepmat(const arma::sp_mat &mat, const int &i, const int &j);
arma::sp_mat FastSparseMatSign(const arma::sp_mat &mat);
arma::sp_mat FastSparseMatMultSD(const arma::sp_mat &mat1, const arma::mat &mat2);
arma::sp_mat FastSparseMatMultDS(const arma::mat &mat1, const arma::sp_mat &mat2);
arma::sp_mat FastSparseMatMultDD(const arma::mat &mat1, const arma::mat &mat2);
arma::sp_mat FastGetRowOfSparseMat(const arma::sp_mat &mat, const int &i);
arma::sp_mat FastGetColOfSparseMat(const arma::sp_mat &mat, const int &j);
arma::sp_mat FastGetRowsOfSparseMat(const arma::sp_mat &mat, const int &start, const int &end);
arma::sp_mat FastGetColsOfSparseMat(const arma::sp_mat &mat, const int &start, const int &end);
arma::sp_mat FastGetSubSparseMat(const arma::sp_mat &mat, const arma::urowvec &rvec, const arma::ucolvec &cvec);
arma::sp_mat FastGetSubSparseMatByRows(const arma::sp_mat &mat, const arma::urowvec &rvec);
arma::sp_mat FastGetSubSparseMatByCols(const arma::sp_mat &mat, const arma::ucolvec &cvec);

#endif //SPARSE_MATRIX_UTIL
