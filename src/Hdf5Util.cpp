#define ARMA_USE_CXX11
#define ARMA_NO_DEBUG
#define ARMA_USE_HDF5

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rhdf5lib)]]
// [[Rcpp::depends(BH)]]
#include "Hdf5Util.h"

//' WriteSpMtAsSpMat
//'
//' This function is used to write a sparse ARMA matrix
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
bool WriteSpMtAsSpMat(const std::string &filePath, const std::string &groupName, const arma::sp_mat &mat) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.WriteSpMtFromArma(mat, groupName);
}

//' WriteSpMtAsS4
//'
//' This function is used to write a sparse S4 matrix
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
void WriteSpMtAsS4(const std::string &filePath, const std::string &groupName, const Rcpp::S4 &mat) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    oHdf5Util.WriteSpMtFromS4(mat, groupName);
}

//' ReadSpMtAsSPMat
//'
//' This function is used to read a sparse matrix from HDF5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
arma::sp_mat ReadSpMtAsSPMat(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.ReadSpMtAsArma(groupName);
}

//' ReadSpMtAsS4
//'
//' This function is used to read a sparse matrix from HDF5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::S4 ReadSpMtAsS4(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.ReadSpMtAsS4(groupName);
}

//' ReadRowSumSpMt
//'
//' Read rows sums
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ReadRowSumSpMt(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    std::vector<double> sumVec;

    std::string datasetName = oHdf5Util.getRowsumDatasetName();
    oHdf5Util.ReadDatasetVector<double>(groupName, datasetName, sumVec);
    if(sumVec.size() == 0) {
        arma::sp_mat mat = oHdf5Util.ReadSpMtAsArma(groupName);
        if(mat.size() == 0) {
            Rcpp::S4 s = oHdf5Util.ReadSpMtAsS4(groupName);
            mat = Rcpp::as<arma::sp_mat>(s);
        }

        sumVec.resize(mat.n_rows);
        for (int i= 0; i< mat.n_cols; i++)
        {
            for (arma::sp_mat::const_col_iterator cij = mat.begin_col(i); cij != mat.end_col(i); ++cij) {
                sumVec[cij.row()] += (*cij);
            }
        }
        oHdf5Util.WriteDatasetVector(groupName, datasetName, sumVec);
    }

    //Rcpp::NumericVector result(sumVec.begin(), sumVec.end());
    return Rcpp::wrap(sumVec);
}

//' ReadColSumSpMt
//'
//' Read cols sums
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ReadColSumSpMt(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    std::vector<double> sumVec;

    std::string datasetName = oHdf5Util.getColsumDatasetName();
    oHdf5Util.ReadDatasetVector(groupName, datasetName, sumVec);
    if(sumVec.size() == 0) {
        arma::sp_mat mat = oHdf5Util.ReadSpMtAsArma(groupName);
        if(mat.size() == 0) {
            Rcpp::S4 s = oHdf5Util.ReadSpMtAsS4(groupName);
            mat = Rcpp::as<arma::sp_mat>(s);
        }

        sumVec.resize(mat.n_cols);
        com::bioturing::SumColumWorker<std::vector<double>> sumColWorker(&mat, sumVec);
        RcppParallel::parallelFor(0, mat.n_cols, sumColWorker);
        oHdf5Util.WriteDatasetVector(groupName, datasetName, sumVec);
    }

    //Rcpp::NumericVector result(sumVec.begin(), sumVec.end());
    //return result;
    return Rcpp::wrap(sumVec);
}

//' GetListAttributes
//'
//' Get list attribute of a dataset
//'
//' @param filePath A HDF5 path
//' @param groupName A string (HDF5 dataset)
//' @param datasetName A dataset name
//' @export
// [[Rcpp::export]]
Rcpp::StringVector GetListAttributes(const std::string &filePath, const std::string &groupName, const std::string &datasetName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    std::vector<std::string> dataVec;
    oHdf5Util.GetListAttributes(groupName, datasetName, dataVec);
    Rcpp::StringVector result(dataVec.begin(), dataVec.end());
    return result;
}

//' GetListObjectNames
//'
//' Get list object of a group
//'
//' @param filePath A HDF5 path
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::StringVector GetListObjectNames(const std::string &filePath, const std::string &groupName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    std::vector<std::string> dataVec;
    oHdf5Util.GetListObjectNames(groupName, dataVec);
    Rcpp::StringVector result(dataVec.begin(), dataVec.end());
    return result;
}

//' GetListRootObjectNames
//'
//' Get list groups
//'
//' @param filePath An HDF5 path
//' @export
// [[Rcpp::export]]
Rcpp::StringVector GetListRootObjectNames(const std::string &filePath) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    std::vector<std::string> dataVec;
    oHdf5Util.GetListRootObjectNames(dataVec);
    Rcpp::StringVector result(dataVec.begin(), dataVec.end());
    return result;
}

//' Read10XH5
//'
//' Get list triplet of SEXP
//'
//' @param s A SEXP type
//' @export
// [[Rcpp::export]]
Rcpp::List Read10XH5Content(const std::string &filePath, const bool &use_names, const bool &unique_features) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    return oHdf5Util.Read10XH5(filePath, use_names, unique_features, false);
}
