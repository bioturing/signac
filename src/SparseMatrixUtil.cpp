#include "Hdf5Util.h"

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]

//' FastCreateSparseMat
//'
//' Generate a sparse matrix with the elements along the main diagonal set to one and off-diagonal elements set to zero
//'
//' @param nrow An integer number
//' @param ncol An integer number
//' @export
// [[Rcpp::export]]
arma::sp_mat FastCreateSparseMat(int nrow, int ncol) {
    return arma::speye(nrow, ncol);
}

//' FastStatsOfSparseMat
//'
//' Get statistic of sparse matrix
//'
//' @param mat An sparse matrix
//' @export
// [[Rcpp::export]]
Rcpp::List FastStatsOfSparseMat(const arma::sp_mat &mat) {
    return Rcpp::List::create(mat.n_rows, mat.n_cols, mat.n_elem, mat.n_nonzero);
}

//' FastCreateFromTriplet
//'
//' Create sparse matrix from triplet
//'
//' @param vec1 An integer vector
//' @param vec2 An integer vector
//' @param vec_val An double vector
//' @export
// [[Rcpp::export]]
arma::sp_mat FastCreateFromTriplet(const arma::urowvec &vec1, const arma::urowvec &vec2, const arma::colvec &vec_val) {
    arma::umat loc = arma::join_vert(vec1, vec2);
    arma::sp_mat sp(loc, vec_val);
    return sp;
}

//' FastConvertToSparseMat
//'
//' Convert SEXP to sparse matrix
//'
//' @param s A SEXP type
//' @export
// [[Rcpp::export]]
arma::sp_mat FastConvertToSparseMat(const SEXP &s) {
    return Rcpp::as<arma::sp_mat>(s);
}

//' FastConvertToTripletMat
//'
//' Get list triplet of SEXP
//'
//' @param s A SEXP type
//' @export
// [[Rcpp::export]]
Rcpp::List FastConvertToTripletMat(const SEXP &s) {
    return Rcpp::simple_triplet_matrix(Rcpp::as<arma::sp_mat>(s));
}

//' FastSparseMatSqrt
//'
//' SQRT element i,j in a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatSqrt(const arma::sp_mat &mat) {
    return arma::sqrt(mat);
}

//' FastSparseMatMult
//'
//' Multiply two sparse matrix
//'
//' @param mat1 A sparse matrix
//' @param mat2 A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatMult(const arma::sp_mat &mat1, const arma::sp_mat &mat2) {
    return mat1 * mat2;
}

//' FastSparseMatAddition
//'
//' Add two sparse matrix
//'
//' @param mat1 A sparse matrix
//' @param mat2 A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatAddition(const arma::sp_mat &mat1, const arma::sp_mat &mat2) {
    return mat1 + mat2;
}

//' FastSparseMatMultWithNum
//'
//' Multiply sparse matrix with number
//'
//' @param mat A sparse matrix
//' @param k A integer number
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatMultWithNum(const arma::sp_mat &mat, const int &k) {
    return k * mat;
}

//' FastSparseMatSymmatl
//'
//' Symmatl sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatSymmatl(const arma::sp_mat &mat) {
    return arma::symmatl(mat);
}

//' FastSparseMatTranspose
//'
//' Transpose sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatTranspose(const arma::sp_mat &mat) {
    return arma::trans(mat);
}

//' FastSparseMatTrimatu
//'
//' Trimatu sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatTrimatu(const arma::sp_mat &mat) {
    return arma::trimatu(mat);
}

//' FastSparseMatTrace
//'
//' Trace sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
int FastSparseMatTrace(const arma::sp_mat &mat) {
    return arma::trace(mat);
}

//' FastConvertToDiagonalSparseMat
//'
//' Get diag sparse matrix from a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastConvertToDiagonalSparseMat(arma::sp_mat &mat) {
    mat.diag().ones();
    return mat;
}

//' FastSparseMatSquare
//'
//' Square a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatSquare(const arma::sp_mat &mat) {
    return arma::square(mat);
}

//' FastSparseMatRepmat
//'
//' Repmat a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatRepmat(const arma::sp_mat &mat, const int &i, const int &j) {
    return arma::repmat(mat, i, j);
}

//' FastSparseMatSign
//'
//' Sign a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatSign(const arma::sp_mat &mat) {
    return arma::sign(mat);
}

//' FastSparseMatMultSD
//'
//' Multiply a sparse matrix with a dense matrix
//'
//' @param mat1 A sparse matrix
//' @param mat2 A dense matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatMultSD(const arma::sp_mat &mat1, const arma::mat &mat2) {
    arma::sp_mat temp2(mat2);
    arma::sp_mat result(mat1 * temp2);
    return result;
}

//' FastSparseMatMultDD
//'
//' Multiply two dense matrix
//'
//' @param mat1 A dense matrix
//' @param mat2 A dense matrix
//' @export
// [[Rcpp::export]]
arma::sp_mat FastSparseMatMultDD(const arma::mat &mat1, const arma::mat &mat2) {
    arma::sp_mat temp1(mat1);
    arma::sp_mat temp2(mat2);
    arma::sp_mat result(temp1 * temp2);
    return result;
}

//' FastGetRowOfSparseMat
//'
//' Get row of sparse matrix
//'
//' @param mat A dense matrix
//' @param i A integer number
//' @export
// [[Rcpp::export]]
arma::sp_mat FastGetRowOfSparseMat(const arma::sp_mat &mat, const int &i) {
    int irow = i;

    try {
        PerformRIndex(i, (int)mat.n_rows, irow);
        return mat.row(irow - 1);
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("Signac exception (unknown reason)");
    }

    arma::sp_mat emat;
    return emat;
}

//' FastGetColOfSparseMat
//'
//' Get col of sparse matrix
//'
//' @param mat A dense matrix
//' @param i A integer number
//' @export
// [[Rcpp::export]]
arma::sp_mat FastGetColOfSparseMat(const arma::sp_mat &mat, const int &j) {
    int icol = j;

    try {
        PerformRIndex(j, (int)mat.n_cols, icol);
        return mat.col(icol);
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("Signac exception (unknown reason)");
    }

    arma::sp_mat emat;
    return emat;
}

//' FastGetRowsOfSparseMat
//'
//' Get rows of sparse matrix
//'
//' @param mat A sparse matrix
//' @param start A integer number
//' @param end A integer number
//' @export
// [[Rcpp::export]]
arma::sp_mat FastGetRowsOfSparseMat(const arma::sp_mat &mat, const int &start, const int &end) {
    int irow_start = start;
    int irow_end = end;

    try {
        PerformRMultiIndex(start, (int)mat.n_rows, irow_start, end, (int)mat.n_rows, irow_end);
        return mat.rows(irow_start, irow_end);
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("Signac exception (unknown reason)");
    }

    arma::sp_mat emat;
    return emat;
}

//' FastGetColsOfSparseMat
//'
//' Get cols of sparse matrix
//'
//' @param mat A sparse matrix
//' @param start A integer number
//' @param end A integer number
//' @export
// [[Rcpp::export]]
arma::sp_mat FastGetColsOfSparseMat(const arma::sp_mat &mat, const int &start, const int &end) {
    int icol_start = start;
    int icol_end = end;

    try {
        PerformRMultiIndex(start, (int)mat.n_cols, icol_start, end, (int)mat.n_cols, icol_end);
        return mat.cols(icol_start, icol_end);
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("Signac exception (unknown reason)");
    }

    arma::sp_mat emat;
    return emat;
}

//' FastGetSubSparseMat
//'
//' Get subview of a sparse matrix
//'
//' @param mat A sparse matrix
//' @param rrvec A row vector
//' @param ccvec A col vector
//' @export
// [[Rcpp::export]]
arma::sp_mat FastGetSubSparseMat(const arma::sp_mat &mat, const arma::urowvec &rrvec, const arma::ucolvec &ccvec) {
    arma::urowvec rvec(rrvec.size());
    PerformRVector(rrvec, (int)mat.n_rows, rvec);
    arma::ucolvec cvec(ccvec.size());
    std::size_t total_rows = rvec.size();
    std::size_t total_cols = cvec.size();
    try {
        PerformRVector(ccvec, (int)mat.n_cols, cvec);

        bool found = false;
        std::size_t n = 0;
        std::size_t p = 0;
        std::size_t found_idx = 0;

        arma::vec new_val(mat.n_nonzero);
        arma::uvec new_rvec(mat.n_nonzero);
        arma::uvec new_cvec(total_cols + 1);
        new_cvec(p) = 0;

        for (auto const& j: cvec) {
            for (std::size_t k = mat.col_ptrs[j]; k < mat.col_ptrs[j + 1]; k++) {
                found = false;
                found_idx = 0;
                while (!found && found_idx < total_rows) {
                    if (mat.row_indices[k] == rvec.at(found_idx)) {
                        found = true;
                    }
                    found_idx++;
                }

                if (found) {
                    new_val(n) = mat.values[k];
                    new_rvec(n) = found_idx - 1;
                    n++;
                }
            }

            p++;
            new_cvec(p) = n;
        }
        new_cvec(p) = n ;

        new_val.reshape(n, 1);
        new_rvec.reshape(n, 1);

        return arma::sp_mat(new_rvec, new_cvec, new_val, total_rows, total_cols);
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("Signac exception (unknown reason)");
    }

    arma::sp_mat emat;
    return emat;
}

//' FastGetSubSparseMatByRows
//'
//' Get subview of a sparse matrix by rows
//'
//' @param mat A sparse matrix
//' @param rvec A row vector
//' @export
// [[Rcpp::export]]
arma::sp_mat FastGetSubSparseMatByRows(const arma::sp_mat &mat, const arma::urowvec &rvec) {
    arma::ucolvec cvec(mat.n_cols);
    for(int i = 0; i< mat.n_cols; i++) {
        cvec(i) = i;
    }
    return FastGetSubSparseMat(mat, rvec, cvec);
}

//' FastGetSubSparseMatByCols
//'
//' Get subview of a sparse matrix by cols
//'
//' @param mat A sparse matrix
//' @param cvec A col vector
//' @export
// [[Rcpp::export]]
arma::sp_mat FastGetSubSparseMatByCols(const arma::sp_mat &mat, const arma::ucolvec &cvec) {
    arma::urowvec rvec(mat.n_rows);
    for(int i = 0; i< mat.n_rows; i++) {
        rvec(i) = i;
    }
    return FastGetSubSparseMat(mat, rvec, cvec);
}

//' FastGetSumSparseMatByRows
//'
//' Sum of rows in a sparse matrix
//'
//' @param mat A sparse matrix
//' @param rvec A col vector
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector FastGetSumSparseMatByRows(const arma::sp_mat &mat, const arma::urowvec &rvec) {
    Rcpp::NumericVector result(rvec.size());
    arma::urowvec rrvec(rvec.size());
    PerformRVector(rvec, (int)mat.n_rows, rrvec);

    for (int i = 0; i< rvec.size(); i++)
    {
        for (arma::sp_mat::const_row_iterator rij = mat.begin_row(i); rij != mat.end_row(i); ++rij) {
            result[rij.row()] += (*rij);
        }
    }

    return result;
}

//' FastGetSumSparseMatByCols
//'
//' Sum of cols in a sparse matrix
//'
//' @param mat A sparse matrix
//' @param cvec A col vector
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector FastGetSumSparseMatByCols(const arma::sp_mat &mat, const arma::ucolvec &cvec) {
    Rcpp::NumericVector result(cvec.size());
    arma::ucolvec ccvec(cvec.size());

    try {
        PerformRVector(cvec, (int)mat.n_cols, ccvec);
        for (int i = 0; i< cvec.size(); i++)
        {
            for (arma::sp_mat::const_col_iterator cij = mat.begin_col(i); cij != mat.end_col(i); ++cij) {
                result[cij.col()] += (*cij);
            }
        }

        return result;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("Signac exception (unknown reason)");
    }

    return result;
}

//' FastGetSumSparseMatByAllRows
//'
//' Sum all rows in a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector FastGetSumSparseMatByAllRows(arma::sp_mat &mat) {
    Rcpp::NumericVector result(mat.n_rows);
    for (int i= 0; i< mat.n_cols; i++)
    {
        for (arma::sp_mat::const_col_iterator cij = mat.begin_col(i); cij != mat.end_col(i); ++cij) {
            result[cij.row()] += (*cij);
        }
    }
    return result;
}

struct SumColumWorker : public RcppParallel::Worker
{
    const arma::sp_mat *input;
    Rcpp::NumericVector &output;

    SumColumWorker(const arma::sp_mat *input, Rcpp::NumericVector &output)
        : input(input), output(output) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (int i= begin; i< end; ++i)
        {
            for (arma::sp_mat::const_col_iterator cij = input->begin_col(i); cij != input->end_col(i); ++cij) {
                output[cij.col()] += (*cij);
            }
        }
    }
};

//' FastGetSumSparseMatByAllCols
//'
//' Sum all cols in a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector FastGetSumSparseMatByAllCols(arma::sp_mat &mat) {
    Rcpp::NumericVector result(mat.n_cols);

    SumColumWorker sumColWorker(&mat, result);
    RcppParallel::parallelFor(0, mat.n_cols, sumColWorker);

    return result;
}

//' FastGetMedianSparseMatByAllRows
//'
//' Median all rows in a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector FastGetMedianSparseMatByAllRows(arma::sp_mat &mat) {
    Rcpp::NumericVector result(mat.n_cols);
    arma::sp_mat tmat = mat.t();
    for (int i= 0; i< tmat.n_cols; i++)
    {
        arma::vec ccvec(tmat.col(i).n_nonzero);
        unsigned int index = 0;
        for (arma::sp_mat::const_col_iterator cij = tmat.begin_col(i); cij != tmat.end_col(i); ++cij) {
            ccvec[index++] = (*cij);
        }
        result[i] += arma::median(ccvec);
    }
    return result;
}

//' FastGetMedianSparseMatByAllCols
//'
//' Median all cols in a sparse matrix
//'
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector FastGetMedianSparseMatByAllCols(arma::sp_mat &mat) {
    Rcpp::NumericVector result(mat.n_cols);
    for (int i= 0; i< mat.n_cols; i++)
    {
        arma::vec ccvec(mat.col(i).n_nonzero);
        unsigned int index = 0;
        for (arma::sp_mat::const_col_iterator cij = mat.begin_col(i); cij != mat.end_col(i); ++cij) {
            ccvec[index++] = (*cij);
        }
        result[i] += arma::median(ccvec);
    }
    return result;
}
