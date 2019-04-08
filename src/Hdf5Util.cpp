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
void WriteSpMtAsSpMat(const std::string &filePath, const std::string &groupName, const arma::sp_mat &mat) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    oHdf5Util.WriteSpMtFromArma(mat, groupName);
}

//' WriteSpMtAsSpMatFromS4
//'
//' This function is used to write a sparse ARMA matrix from S4
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @param mat A sparse matrix
//' @export
// [[Rcpp::export]]
void WriteSpMtAsSpMatFromS4(const std::string &filePath, const std::string &groupName, const Rcpp::S4 &mat) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(-1);
    oHdf5Util.WriteSpMtFromArma(file, mat, groupName);
    oHdf5Util.Close(file);
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
    HighFive::File *file = oHdf5Util.Open(-1);
    oHdf5Util.WriteSpMtFromS4(file, mat, groupName);
    oHdf5Util.Close(file);
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
    HighFive::File *file = oHdf5Util.Open(1);
    Rcpp::S4 s = oHdf5Util.ReadSpMtAsS4(file, groupName);
    oHdf5Util.Close(file);
    return s;
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
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<double> sumVec;

    std::string datasetName = oHdf5Util.getRowsumDatasetName();
    oHdf5Util.ReadDatasetVector<double>(file, groupName, datasetName, sumVec);
    if(sumVec.size() == 0) {
        arma::sp_mat mat = oHdf5Util.ReadSpMtAsArma(groupName);
        if(mat.size() == 0) {
            Rcpp::S4 s = oHdf5Util.ReadSpMtAsS4(file, groupName);
            mat = Rcpp::as<arma::sp_mat>(s);
        }

        sumVec.resize(mat.n_rows);
        for (int i= 0; i< mat.n_cols; i++)
        {
            for (arma::sp_mat::const_col_iterator cij = mat.begin_col(i); cij != mat.end_col(i); ++cij) {
                sumVec[cij.row()] += (*cij);
            }
        }
        oHdf5Util.WriteDatasetVector(file, groupName, datasetName, sumVec);
    }

    oHdf5Util.Close(file);
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
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<double> sumVec;

    std::string datasetName = oHdf5Util.getColsumDatasetName();
    oHdf5Util.ReadDatasetVector(file, groupName, datasetName, sumVec);
    if(sumVec.size() == 0) {
        arma::sp_mat mat = oHdf5Util.ReadSpMtAsArma(groupName);
        if(mat.size() == 0) {
            Rcpp::S4 s = oHdf5Util.ReadSpMtAsS4(file, groupName);
            mat = Rcpp::as<arma::sp_mat>(s);
        }

        sumVec.resize(mat.n_cols);
        com::bioturing::SumColumWorker<std::vector<double>> sumColWorker(&mat, sumVec);
        RcppParallel::parallelFor(0, mat.n_cols, sumColWorker);
        oHdf5Util.WriteDatasetVector(file, groupName, datasetName, sumVec);
    }

    oHdf5Util.Close(file);
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
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<std::string> dataVec;
    oHdf5Util.GetListAttributes(file, groupName, datasetName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
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
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<std::string> dataVec;
    oHdf5Util.GetListObjectNames(file, groupName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
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
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<std::string> dataVec;
    oHdf5Util.GetListRootObjectNames(file, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
}

//' Read10XH5
//'
//' Get list triplet of SEXP
//'
//' @param filePath A fiel path
//' @param use_names Use names flag
//' @param unique_features Unique features flag
//' @export
// [[Rcpp::export]]
Rcpp::List Read10XH5Content(const std::string &filePath, const bool &use_names, const bool &unique_features) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(1);
    Rcpp::List arrData = oHdf5Util.Read10XH5(file, filePath, use_names, unique_features);
    oHdf5Util.Close(file);
    return arrData;
}

//' WriteRootDataset
//'
//' Write a string vector to root group
//'
//' @param filePath A fiel path
//' @param datasetName A dataset name
//' @param datasetVal A string vector
//' @export
// [[Rcpp::export]]
void WriteRootDataset(const std::string &filePath, const std::string &datasetName, const std::vector<std::string> &datasetVal) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(-1);
    oHdf5Util.WriteRootDataset(file, datasetName, datasetVal);
    oHdf5Util.Close(file);
}

//' ReadRootDataset
//'
//' Read a string vector from root group
//'
//' @param filePath A fiel path
//' @param datasetName A dataset name
//' @export
// [[Rcpp::export]]
Rcpp::StringVector ReadRootDataset(const std::string &filePath, const std::string &datasetName) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<std::string> dataVec;
    oHdf5Util.ReadRootDataset(file, datasetName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
}

<<<<<<< 32cd1c3f8b12ba0cb6431f2769b3f82bcf1ec79a
<<<<<<< 12516e0bde4aa096932420c0ddd8c6c401ba4fea
//' ReadIntegerVector
//'
//' This function is used to read a integer vector from hdf5 file
=======
//' ReadH5VectorRange
//'
//' This function is used to read a vector[i:j] from hdf5 file
>>>>>>> Update HDF5 function
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
<<<<<<< 12516e0bde4aa096932420c0ddd8c6c401ba4fea
Rcpp::NumericVector ReadIntegerVector(const std::string &filePath,
                                      const std::string &groupName,
                                      const std::string &datasetName) {
    std::vector<int> dataVec;
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(1);
    oHdf5Util.ReadDatasetVector<int>(file, groupName, datasetName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
}

//' ReadDoubleVector
//'
//' This function is used to read a double vector from hdf5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ReadDoubleVector(const std::string &filePath,
                                     const std::string &groupName,
                                     const std::string &datasetName) {
    std::vector<double> dataVec;
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(1);
    oHdf5Util.ReadDatasetVector<double>(file, groupName, datasetName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
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

    std::string vecGroupName = groupName + "_" + "rowsum";
    oHdf5Util.ReadVector(sumVec, vecGroupName);
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
        oHdf5Util.WriteVector(sumVec, vecGroupName);
    }

    Rcpp::NumericVector result(sumVec.begin(), sumVec.end());
    return result;
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

    std::string vecGroupName = groupName + "_" + "colsum";
    oHdf5Util.ReadVector(sumVec, vecGroupName);
    if(sumVec.size() == 0) {
        arma::sp_mat mat = oHdf5Util.ReadSpMtAsArma(groupName);
        if(mat.size() == 0) {
            Rcpp::S4 s = oHdf5Util.ReadSpMtAsS4(groupName);
            mat = Rcpp::as<arma::sp_mat>(s);
        }

        sumVec.resize(mat.n_cols);
        com::bioturing::SumColumWorker<std::vector<double>> sumColWorker(&mat, sumVec);
        RcppParallel::parallelFor(0, mat.n_cols, sumColWorker);
        oHdf5Util.WriteVector(sumVec, vecGroupName);
    }

    Rcpp::NumericVector result(sumVec.begin(), sumVec.end());
    return result;
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
    oHdf5Util.GetListAttributes(file, groupName, datasetName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
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
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<std::string> dataVec;
    oHdf5Util.GetListObjectNames(file, groupName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
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
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<std::string> dataVec;
    oHdf5Util.GetListRootObjectNames(file, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
}

//' Read10XH5
//'
//' Get list triplet of SEXP
//'
//' @param filePath A fiel path
//' @param use_names Use names flag
//' @param unique_features Unique features flag
//' @export
// [[Rcpp::export]]
Rcpp::List Read10XH5Content(const std::string &filePath, const bool &use_names, const bool &unique_features) {
    return Rcpp::List::create();
=======
void ReadGeneExpH5(const std::string &filePath,
                   const std::string &groupName,
                   const int g_idx,
                   std::vector<int> &col_idx,
                   std::vector<double> &g_exp) {
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(1);
    std::vector<int> vec;
    oHdf5Util.ReadDatasetVector<int>(file, groupName, "indptr", vec);
    oHdf5Util.ReadDatasetRangeVector<int>(file, groupName, "indices", (unsigned int)vec[g_idx],
                      (unsigned int)vec[g_idx + 1], col_idx);
    oHdf5Util.ReadDatasetRangeVector<double>(file, groupName, "data", (unsigned int)vec[g_idx],
                      (unsigned int)vec[g_idx + 1], g_exp);
    oHdf5Util.Close(file);
>>>>>>> Update HDF5 function
}

=======
>>>>>>> Update hdf5 FOR CLEANING
//' ReadIntegerVector
//'
//' This function is used to read a integer vector from hdf5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ReadIntegerVector(const std::string &filePath,
                                      const std::string &groupName,
                                      const std::string &datasetName) {
    std::vector<int> dataVec;
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(1);
    oHdf5Util.ReadDatasetVector<int>(file, groupName, datasetName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
}

//' ReadDoubleVector
//'
//' This function is used to read a double vector from hdf5 file
//'
//' @param filePath A string (HDF5 path)
//' @param groupName A string (HDF5 dataset)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ReadDoubleVector(const std::string &filePath,
                                     const std::string &groupName,
                                     const std::string &datasetName) {
    std::vector<double> dataVec;
    com::bioturing::Hdf5Util oHdf5Util(filePath);
    HighFive::File *file = oHdf5Util.Open(1);
    oHdf5Util.ReadDatasetVector<double>(file, groupName, datasetName, dataVec);
    oHdf5Util.Close(file);
    return Rcpp::wrap(dataVec);
}
