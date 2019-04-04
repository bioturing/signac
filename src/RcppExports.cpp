// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "Signac_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// FastGetCurrentDate
Rcpp::Date FastGetCurrentDate();
RcppExport SEXP _Signac_FastGetCurrentDate() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(FastGetCurrentDate());
    return rcpp_result_gen;
END_RCPP
}
// FastDiffVector
arma::uvec FastDiffVector(const arma::uvec& a, const arma::uvec& b);
RcppExport SEXP _Signac_FastDiffVector(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(FastDiffVector(a, b));
    return rcpp_result_gen;
END_RCPP
}
// FastRandVector
arma::uvec FastRandVector(int num);
RcppExport SEXP _Signac_FastRandVector(SEXP numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type num(numSEXP);
    rcpp_result_gen = Rcpp::wrap(FastRandVector(num));
    return rcpp_result_gen;
END_RCPP
}
// HarmonyMarker
DataFrame HarmonyMarker(Rcpp::S4& S4_mtx, const Rcpp::NumericVector& in_idx, const Rcpp::Nullable<Rcpp::NumericVector>& out_idx);
RcppExport SEXP _Signac_HarmonyMarker(SEXP S4_mtxSEXP, SEXP in_idxSEXP, SEXP out_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type S4_mtx(S4_mtxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type in_idx(in_idxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericVector>& >::type out_idx(out_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(HarmonyMarker(S4_mtx, in_idx, out_idx));
    return rcpp_result_gen;
END_RCPP
}
// WriteSpMtV2
bool WriteSpMtV2(const std::string& filePath, const std::string& groupName, const arma::sp_mat& mat);
RcppExport SEXP _Signac_WriteSpMtV2(SEXP filePathSEXP, SEXP groupNameSEXP, SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filePath(filePathSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type groupName(groupNameSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteSpMtV2(filePath, groupName, mat));
    return rcpp_result_gen;
END_RCPP
}
// WriteSpMtV1
bool WriteSpMtV1(const std::string& filePath, const std::string& groupName, const Rcpp::S4& mat);
RcppExport SEXP _Signac_WriteSpMtV1(SEXP filePathSEXP, SEXP groupNameSEXP, SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filePath(filePathSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type groupName(groupNameSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteSpMtV1(filePath, groupName, mat));
    return rcpp_result_gen;
END_RCPP
}
// ReadSpMtV2
arma::sp_mat ReadSpMtV2(const std::string& filePath, const std::string& groupName);
RcppExport SEXP _Signac_ReadSpMtV2(SEXP filePathSEXP, SEXP groupNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filePath(filePathSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type groupName(groupNameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadSpMtV2(filePath, groupName));
    return rcpp_result_gen;
END_RCPP
}
// ReadSpMtV1
Rcpp::S4 ReadSpMtV1(const std::string& filePath, const std::string& groupName);
RcppExport SEXP _Signac_ReadSpMtV1(SEXP filePathSEXP, SEXP groupNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filePath(filePathSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type groupName(groupNameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadSpMtV1(filePath, groupName));
    return rcpp_result_gen;
END_RCPP
}
// ReadRowSumSpMt
Rcpp::NumericVector ReadRowSumSpMt(const std::string& filePath, const std::string& groupName);
RcppExport SEXP _Signac_ReadRowSumSpMt(SEXP filePathSEXP, SEXP groupNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filePath(filePathSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type groupName(groupNameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadRowSumSpMt(filePath, groupName));
    return rcpp_result_gen;
END_RCPP
}
// ReadColSumSpMt
Rcpp::NumericVector ReadColSumSpMt(const std::string& filePath, const std::string& groupName);
RcppExport SEXP _Signac_ReadColSumSpMt(SEXP filePathSEXP, SEXP groupNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filePath(filePathSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type groupName(groupNameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadColSumSpMt(filePath, groupName));
    return rcpp_result_gen;
END_RCPP
}
// FastMatMult
arma::mat FastMatMult(const arma::mat& mat1, const arma::mat& mat2);
RcppExport SEXP _Signac_FastMatMult(SEXP mat1SEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(FastMatMult(mat1, mat2));
    return rcpp_result_gen;
END_RCPP
}
// FastGetRowsOfMat
arma::mat FastGetRowsOfMat(const arma::mat& mat, arma::uvec vec);
RcppExport SEXP _Signac_FastGetRowsOfMat(SEXP matSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetRowsOfMat(mat, vec));
    return rcpp_result_gen;
END_RCPP
}
// FastGetColsOfMat
arma::mat FastGetColsOfMat(const arma::mat& mat, arma::uvec vec);
RcppExport SEXP _Signac_FastGetColsOfMat(SEXP matSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetColsOfMat(mat, vec));
    return rcpp_result_gen;
END_RCPP
}
// FastGetSubMat
arma::mat FastGetSubMat(const arma::mat& mat, arma::uvec rvec, arma::uvec cvec);
RcppExport SEXP _Signac_FastGetSubMat(SEXP matSEXP, SEXP rvecSEXP, SEXP cvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type rvec(rvecSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cvec(cvecSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetSubMat(mat, rvec, cvec));
    return rcpp_result_gen;
END_RCPP
}
// FastCreateSparseMat
arma::sp_mat FastCreateSparseMat(int nrow, int ncol);
RcppExport SEXP _Signac_FastCreateSparseMat(SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(FastCreateSparseMat(nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// FastStatsOfSparseMat
Rcpp::List FastStatsOfSparseMat(const arma::sp_mat& mat);
RcppExport SEXP _Signac_FastStatsOfSparseMat(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastStatsOfSparseMat(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastCreateFromTriplet
arma::sp_mat FastCreateFromTriplet(const arma::urowvec& vec1, const arma::urowvec& vec2, const arma::colvec& vec_val);
RcppExport SEXP _Signac_FastCreateFromTriplet(SEXP vec1SEXP, SEXP vec2SEXP, SEXP vec_valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::urowvec& >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< const arma::urowvec& >::type vec2(vec2SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type vec_val(vec_valSEXP);
    rcpp_result_gen = Rcpp::wrap(FastCreateFromTriplet(vec1, vec2, vec_val));
    return rcpp_result_gen;
END_RCPP
}
// FastConvertToSparseMat
arma::sp_mat FastConvertToSparseMat(const SEXP& s);
RcppExport SEXP _Signac_FastConvertToSparseMat(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(FastConvertToSparseMat(s));
    return rcpp_result_gen;
END_RCPP
}
// FastConvertToTripletMat
Rcpp::List FastConvertToTripletMat(const SEXP& s);
RcppExport SEXP _Signac_FastConvertToTripletMat(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(FastConvertToTripletMat(s));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatSqrt
arma::sp_mat FastSparseMatSqrt(const arma::sp_mat& mat);
RcppExport SEXP _Signac_FastSparseMatSqrt(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatSqrt(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatMult
arma::sp_mat FastSparseMatMult(const arma::sp_mat& mat1, const arma::sp_mat& mat2);
RcppExport SEXP _Signac_FastSparseMatMult(SEXP mat1SEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatMult(mat1, mat2));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatAddition
arma::sp_mat FastSparseMatAddition(const arma::sp_mat& mat1, const arma::sp_mat& mat2);
RcppExport SEXP _Signac_FastSparseMatAddition(SEXP mat1SEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatAddition(mat1, mat2));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatMultWithNum
arma::sp_mat FastSparseMatMultWithNum(const arma::sp_mat& mat, const int& k);
RcppExport SEXP _Signac_FastSparseMatMultWithNum(SEXP matSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatMultWithNum(mat, k));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatSymmatl
arma::sp_mat FastSparseMatSymmatl(const arma::sp_mat& mat);
RcppExport SEXP _Signac_FastSparseMatSymmatl(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatSymmatl(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatTranspose
arma::sp_mat FastSparseMatTranspose(const arma::sp_mat& mat);
RcppExport SEXP _Signac_FastSparseMatTranspose(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatTranspose(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatTrimatu
arma::sp_mat FastSparseMatTrimatu(const arma::sp_mat& mat);
RcppExport SEXP _Signac_FastSparseMatTrimatu(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatTrimatu(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatTrace
int FastSparseMatTrace(const arma::sp_mat& mat);
RcppExport SEXP _Signac_FastSparseMatTrace(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatTrace(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastConvertToDiagonalSparseMat
arma::sp_mat FastConvertToDiagonalSparseMat(arma::sp_mat& mat);
RcppExport SEXP _Signac_FastConvertToDiagonalSparseMat(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastConvertToDiagonalSparseMat(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatSquare
arma::sp_mat FastSparseMatSquare(const arma::sp_mat& mat);
RcppExport SEXP _Signac_FastSparseMatSquare(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatSquare(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatRepmat
arma::sp_mat FastSparseMatRepmat(const arma::sp_mat& mat, const int& i, const int& j);
RcppExport SEXP _Signac_FastSparseMatRepmat(SEXP matSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const int& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const int& >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatRepmat(mat, i, j));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatSign
arma::sp_mat FastSparseMatSign(const arma::sp_mat& mat);
RcppExport SEXP _Signac_FastSparseMatSign(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatSign(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatMultSD
arma::sp_mat FastSparseMatMultSD(const arma::sp_mat& mat1, const arma::mat& mat2);
RcppExport SEXP _Signac_FastSparseMatMultSD(SEXP mat1SEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatMultSD(mat1, mat2));
    return rcpp_result_gen;
END_RCPP
}
// FastSparseMatMultDD
arma::sp_mat FastSparseMatMultDD(const arma::mat& mat1, const arma::mat& mat2);
RcppExport SEXP _Signac_FastSparseMatMultDD(SEXP mat1SEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(FastSparseMatMultDD(mat1, mat2));
    return rcpp_result_gen;
END_RCPP
}
// FastGetRowOfSparseMat
arma::sp_mat FastGetRowOfSparseMat(const arma::sp_mat& mat, const int& i);
RcppExport SEXP _Signac_FastGetRowOfSparseMat(SEXP matSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const int& >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetRowOfSparseMat(mat, i));
    return rcpp_result_gen;
END_RCPP
}
// FastGetColOfSparseMat
arma::sp_mat FastGetColOfSparseMat(const arma::sp_mat& mat, const int& j);
RcppExport SEXP _Signac_FastGetColOfSparseMat(SEXP matSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const int& >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetColOfSparseMat(mat, j));
    return rcpp_result_gen;
END_RCPP
}
// FastGetRowsOfSparseMat
arma::sp_mat FastGetRowsOfSparseMat(const arma::sp_mat& mat, const int& start, const int& end);
RcppExport SEXP _Signac_FastGetRowsOfSparseMat(SEXP matSEXP, SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const int& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const int& >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetRowsOfSparseMat(mat, start, end));
    return rcpp_result_gen;
END_RCPP
}
// FastGetColsOfSparseMat
arma::sp_mat FastGetColsOfSparseMat(const arma::sp_mat& mat, const int& start, const int& end);
RcppExport SEXP _Signac_FastGetColsOfSparseMat(SEXP matSEXP, SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const int& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const int& >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetColsOfSparseMat(mat, start, end));
    return rcpp_result_gen;
END_RCPP
}
// FastGetSubSparseMat
arma::sp_mat FastGetSubSparseMat(const arma::sp_mat& mat, const arma::urowvec& rrvec, const arma::ucolvec& ccvec, const bool& need_perform_row, const bool& need_perform_col);
RcppExport SEXP _Signac_FastGetSubSparseMat(SEXP matSEXP, SEXP rrvecSEXP, SEXP ccvecSEXP, SEXP need_perform_rowSEXP, SEXP need_perform_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::urowvec& >::type rrvec(rrvecSEXP);
    Rcpp::traits::input_parameter< const arma::ucolvec& >::type ccvec(ccvecSEXP);
    Rcpp::traits::input_parameter< const bool& >::type need_perform_row(need_perform_rowSEXP);
    Rcpp::traits::input_parameter< const bool& >::type need_perform_col(need_perform_colSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetSubSparseMat(mat, rrvec, ccvec, need_perform_row, need_perform_col));
    return rcpp_result_gen;
END_RCPP
}
// FastGetSubSparseMatByRows
arma::sp_mat FastGetSubSparseMatByRows(const arma::sp_mat& mat, const arma::urowvec& rvec);
RcppExport SEXP _Signac_FastGetSubSparseMatByRows(SEXP matSEXP, SEXP rvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::urowvec& >::type rvec(rvecSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetSubSparseMatByRows(mat, rvec));
    return rcpp_result_gen;
END_RCPP
}
// FastGetSubSparseMatByCols
arma::sp_mat FastGetSubSparseMatByCols(const arma::sp_mat& mat, const arma::ucolvec& cvec);
RcppExport SEXP _Signac_FastGetSubSparseMatByCols(SEXP matSEXP, SEXP cvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::ucolvec& >::type cvec(cvecSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetSubSparseMatByCols(mat, cvec));
    return rcpp_result_gen;
END_RCPP
}
// FastGetSumSparseMatByRows
Rcpp::NumericVector FastGetSumSparseMatByRows(const arma::sp_mat& mat, const arma::urowvec& rvec);
RcppExport SEXP _Signac_FastGetSumSparseMatByRows(SEXP matSEXP, SEXP rvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::urowvec& >::type rvec(rvecSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetSumSparseMatByRows(mat, rvec));
    return rcpp_result_gen;
END_RCPP
}
// FastGetSumSparseMatByCols
Rcpp::NumericVector FastGetSumSparseMatByCols(const arma::sp_mat& mat, const arma::ucolvec& cvec);
RcppExport SEXP _Signac_FastGetSumSparseMatByCols(SEXP matSEXP, SEXP cvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::ucolvec& >::type cvec(cvecSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetSumSparseMatByCols(mat, cvec));
    return rcpp_result_gen;
END_RCPP
}
// FastGetSumSparseMatByAllRows
Rcpp::NumericVector FastGetSumSparseMatByAllRows(arma::sp_mat& mat);
RcppExport SEXP _Signac_FastGetSumSparseMatByAllRows(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetSumSparseMatByAllRows(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastGetSumSparseMatByAllCols
Rcpp::NumericVector FastGetSumSparseMatByAllCols(arma::sp_mat& mat);
RcppExport SEXP _Signac_FastGetSumSparseMatByAllCols(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetSumSparseMatByAllCols(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastGetMedianSparseMatByAllRows
Rcpp::NumericVector FastGetMedianSparseMatByAllRows(arma::sp_mat& mat);
RcppExport SEXP _Signac_FastGetMedianSparseMatByAllRows(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetMedianSparseMatByAllRows(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastGetMedianSparseMatByAllCols
Rcpp::NumericVector FastGetMedianSparseMatByAllCols(arma::sp_mat& mat);
RcppExport SEXP _Signac_FastGetMedianSparseMatByAllCols(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastGetMedianSparseMatByAllCols(mat));
    return rcpp_result_gen;
END_RCPP
}
// FastConvertS4ToSparseMT
arma::sp_mat FastConvertS4ToSparseMT(Rcpp::S4& mat);
RcppExport SEXP _Signac_FastConvertS4ToSparseMT(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(FastConvertS4ToSparseMT(mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Signac_FastGetCurrentDate", (DL_FUNC) &_Signac_FastGetCurrentDate, 0},
    {"_Signac_FastDiffVector", (DL_FUNC) &_Signac_FastDiffVector, 2},
    {"_Signac_FastRandVector", (DL_FUNC) &_Signac_FastRandVector, 1},
    {"_Signac_Harmony", (DL_FUNC) &_Signac_Harmony, 3},
    {"_Signac_WriteSpMtV2", (DL_FUNC) &_Signac_WriteSpMtV2, 3},
    {"_Signac_WriteSpMtV1", (DL_FUNC) &_Signac_WriteSpMtV1, 3},
    {"_Signac_ReadSpMtV2", (DL_FUNC) &_Signac_ReadSpMtV2, 2},
    {"_Signac_ReadSpMtV1", (DL_FUNC) &_Signac_ReadSpMtV1, 2},
    {"_Signac_ReadRowSumSpMt", (DL_FUNC) &_Signac_ReadRowSumSpMt, 2},
    {"_Signac_ReadColSumSpMt", (DL_FUNC) &_Signac_ReadColSumSpMt, 2},
    {"_Signac_FastMatMult", (DL_FUNC) &_Signac_FastMatMult, 2},
    {"_Signac_FastGetRowsOfMat", (DL_FUNC) &_Signac_FastGetRowsOfMat, 2},
    {"_Signac_FastGetColsOfMat", (DL_FUNC) &_Signac_FastGetColsOfMat, 2},
    {"_Signac_FastGetSubMat", (DL_FUNC) &_Signac_FastGetSubMat, 3},
    {"_Signac_FastCreateSparseMat", (DL_FUNC) &_Signac_FastCreateSparseMat, 2},
    {"_Signac_FastStatsOfSparseMat", (DL_FUNC) &_Signac_FastStatsOfSparseMat, 1},
    {"_Signac_FastCreateFromTriplet", (DL_FUNC) &_Signac_FastCreateFromTriplet, 3},
    {"_Signac_FastConvertToSparseMat", (DL_FUNC) &_Signac_FastConvertToSparseMat, 1},
    {"_Signac_FastConvertToTripletMat", (DL_FUNC) &_Signac_FastConvertToTripletMat, 1},
    {"_Signac_FastSparseMatSqrt", (DL_FUNC) &_Signac_FastSparseMatSqrt, 1},
    {"_Signac_FastSparseMatMult", (DL_FUNC) &_Signac_FastSparseMatMult, 2},
    {"_Signac_FastSparseMatAddition", (DL_FUNC) &_Signac_FastSparseMatAddition, 2},
    {"_Signac_FastSparseMatMultWithNum", (DL_FUNC) &_Signac_FastSparseMatMultWithNum, 2},
    {"_Signac_FastSparseMatSymmatl", (DL_FUNC) &_Signac_FastSparseMatSymmatl, 1},
    {"_Signac_FastSparseMatTranspose", (DL_FUNC) &_Signac_FastSparseMatTranspose, 1},
    {"_Signac_FastSparseMatTrimatu", (DL_FUNC) &_Signac_FastSparseMatTrimatu, 1},
    {"_Signac_FastSparseMatTrace", (DL_FUNC) &_Signac_FastSparseMatTrace, 1},
    {"_Signac_FastConvertToDiagonalSparseMat", (DL_FUNC) &_Signac_FastConvertToDiagonalSparseMat, 1},
    {"_Signac_FastSparseMatSquare", (DL_FUNC) &_Signac_FastSparseMatSquare, 1},
    {"_Signac_FastSparseMatRepmat", (DL_FUNC) &_Signac_FastSparseMatRepmat, 3},
    {"_Signac_FastSparseMatSign", (DL_FUNC) &_Signac_FastSparseMatSign, 1},
    {"_Signac_FastSparseMatMultSD", (DL_FUNC) &_Signac_FastSparseMatMultSD, 2},
    {"_Signac_FastSparseMatMultDD", (DL_FUNC) &_Signac_FastSparseMatMultDD, 2},
    {"_Signac_FastGetRowOfSparseMat", (DL_FUNC) &_Signac_FastGetRowOfSparseMat, 2},
    {"_Signac_FastGetColOfSparseMat", (DL_FUNC) &_Signac_FastGetColOfSparseMat, 2},
    {"_Signac_FastGetRowsOfSparseMat", (DL_FUNC) &_Signac_FastGetRowsOfSparseMat, 3},
    {"_Signac_FastGetColsOfSparseMat", (DL_FUNC) &_Signac_FastGetColsOfSparseMat, 3},
    {"_Signac_FastGetSubSparseMat", (DL_FUNC) &_Signac_FastGetSubSparseMat, 5},
    {"_Signac_FastGetSubSparseMatByRows", (DL_FUNC) &_Signac_FastGetSubSparseMatByRows, 2},
    {"_Signac_FastGetSubSparseMatByCols", (DL_FUNC) &_Signac_FastGetSubSparseMatByCols, 2},
    {"_Signac_FastGetSumSparseMatByRows", (DL_FUNC) &_Signac_FastGetSumSparseMatByRows, 2},
    {"_Signac_FastGetSumSparseMatByCols", (DL_FUNC) &_Signac_FastGetSumSparseMatByCols, 2},
    {"_Signac_FastGetSumSparseMatByAllRows", (DL_FUNC) &_Signac_FastGetSumSparseMatByAllRows, 1},
    {"_Signac_FastGetSumSparseMatByAllCols", (DL_FUNC) &_Signac_FastGetSumSparseMatByAllCols, 1},
    {"_Signac_FastGetMedianSparseMatByAllRows", (DL_FUNC) &_Signac_FastGetMedianSparseMatByAllRows, 1},
    {"_Signac_FastGetMedianSparseMatByAllCols", (DL_FUNC) &_Signac_FastGetMedianSparseMatByAllCols, 1},
    {"_Signac_FastConvertS4ToSparseMT", (DL_FUNC) &_Signac_FastConvertS4ToSparseMT, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_Signac(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
