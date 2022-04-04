#ifndef HDF5_UTIL
#define HDF5_UTIL

#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif
#define GENOME_SEPARATOR "_"

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include "CommonUtil.h"
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <H5Cpp.h>

using namespace H5;
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

namespace com {
namespace bioturing {

template <typename T>
struct SumColumWorker : public RcppParallel::Worker
{
    const arma::sp_mat *input;
    T &output;

    SumColumWorker(const arma::sp_mat *input, T &output)
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

class Hdf5Util {
public:
    Hdf5Util(const std::string &file_name_) {
        file_name = std::string(file_name_.data(), file_name_.size());
    }

    ~Hdf5Util() {}

    std::string getRowsumDatasetName() {
        return "rowsums";
    }

    std::string getColsumDatasetName() {
        return "colsums";
    }

    template <typename T>
    void WriteDatasetVector(HighFive::File *file, const std::string &groupName, const std::string &datasetName, const std::vector<T> &datasetVec) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not write dataset, open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == true) {
                std::stringstream ostr;
                ostr << "Existing group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == true) {
                std::stringstream ostr;
                ostr << "Existing dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            HighFive::DataSet datasetGroup = file->createDataSet<T>(groupName, HighFive::DataSpace::From(datasetVec));
            datasetGroup.write(datasetVec);
            file->flush();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteVector HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    template <typename T>
    void ReadDatasetVector(HighFive::File *file, const std::string &groupName, const std::string &datasetName, std::vector<T> &vvec) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read dataset, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::DataSet datasetVec = file->getDataSet(groupName + "/" + datasetName);

#ifdef DEBUG
            HighFive::DataSpace dataSpace = datasetVec.getSpace();
            auto dims = dataSpace.getDimensions();
            HighFive::DataType dataType = datasetVec.getDataType();
            size_t str_size = H5Tget_size(dataType.getId());
            std::stringstream ostr;
            ostr << "Read dataset (T) : [" << groupName << "/" << datasetName << "] with str_size=" << str_size << ", dims[0]=" << dims[0];
            ::Rf_warning(ostr.str().c_str());
#endif

            datasetVec.read(vvec);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadDatatypeVector (T) HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    template<size_t MAX_LEN>
    void ReadDatasetVector(HighFive::File *file, const std::string &groupName, const std::string &datasetName, std::vector<std::string> &vvec) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read dataset, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            HighFive::DataSet dset = file->getDataSet(groupName + "/" + datasetName);
            if (dset.getDataType().isFixedLenStr()) {
                if (dset.getDataType().getSize() > MAX_LEN) {
                    ::Rf_error("Dataset length exceeds the specified maxiumum length");
                    Close(file);
                    throw;
                }
                HighFive::FixedLenStringArray<MAX_LEN> tmp;
                dset.read(tmp);
                vvec.resize(tmp.size());
                for (int i = 0; i < tmp.size(); ++i)
                    vvec[i] = tmp[i];
            } else if (dset.getDataType().isVariableStr()) {
                dset.read(vvec);
            } else {
                ::Rf_error("Dataset is neither a fixed-length strings array nor a variable-length strings array");
                Close(file);
                throw;
            }
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadDatatypeVector (STRING) HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    void ReadDatasetInfo(HighFive::File *file, const std::string &groupName, const std::string &datasetName, std::vector<unsigned int> &vvec) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read info, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            HighFive::DataSet datasetVec = file->getDataSet(groupName + "/" + datasetName);
            HighFive::DataSpace dataSpace = datasetVec.getSpace();
            auto ndim = dataSpace.getNumberDimensions();
            auto dims = dataSpace.getDimensions();
            HighFive::DataType dataType = datasetVec.getDataType();
            size_t str_size = H5Tget_size(dataType.getId());
            H5T_str_t str_pad = H5Tget_strpad(dataType.getId());
            H5T_cset_t str_cset = H5Tget_cset(dataType.getId());

            vvec.push_back(datasetVec.getStorageSize());
            vvec.push_back(ndim);
            vvec.push_back(dims[0]);
            vvec.push_back(dims[1]);
            vvec.push_back(str_size);
            vvec.push_back(str_pad);
            vvec.push_back(str_cset);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadDatatypeInfo HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    template <typename T>
    void WriteRootDataset(HighFive::File *file, const std::string &datasetName, const std::vector<T> &datasetVal) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not write dataset, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(datasetName) == true) {
                std::stringstream ostr;
                ostr << "Existing group :" << datasetName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            HighFive::DataSet datasetGroup = file->createDataSet<T>(datasetName, HighFive::DataSpace::From(datasetVal));
            datasetGroup.write(datasetVal);
            file->flush();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteRootDataset HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    template <typename T>
    void ReadRootDataset(HighFive::File *file, const std::string &datasetName, std::vector<T> &datasetVal) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read dataset, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            HighFive::DataSet datasetVec = file->getDataSet(datasetName);

#ifdef DEBUG
            HighFive::DataSpace dataSpace = datasetVec.getSpace();
            auto dims = dataSpace.getDimensions();
            HighFive::DataType dataType = datasetVec.getDataType();
            size_t str_size = H5Tget_size(dataType.getId());
            std::stringstream ostr;
            ostr << "Read dataset (T) : [" << datasetName << "] with str_size=" << str_size << ", dims[0]=" << dims[0];
            ::Rf_warning(ostr.str().c_str());
#endif

            datasetVec.read(datasetVal);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadRootDataset HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    void ReadRootDataset(HighFive::File *file, const std::string &datasetName, std::vector<std::string> &datasetVal) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read dataset, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }
            //TODO: Change this?
            ReadDatasetVector<256>(file, "", datasetName, datasetVal);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadRootDataset HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    void GetListAttributes(HighFive::File *file, const std::string &groupName, const std::string &datasetName, std::vector<std::string> &arrList) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read attributes, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            HighFive::DataSet datasetData = file->getDataSet(groupName + "/" + datasetName);
            arrList = datasetData.listAttributeNames();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "GetListAttributes HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    void GetListRootObjectNames(HighFive::File *file, std::vector<std::string> &arrList) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read object names, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            arrList = file->listObjectNames();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "GetListObjectNames HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    bool CheckRootAttribute(HighFive::File *file, const std::string &attName) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not check attribute, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        bool isExist = false;
        try {
            isExist = file->hasAttribute(attName);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "CheckRootAttribute HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
        return isExist;
    }

    void GetListRootAttributes(HighFive::File *file, std::vector<std::string> &arrList) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read root attributes, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            arrList = file->listAttributeNames();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "GetListRootAttributes HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    void GetListObjectNames(HighFive::File *file, const std::string &groupName, std::vector<std::string> &arrList) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read objects, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            HighFive::Group datasetGroup = file->getGroup(groupName);
            arrList = datasetGroup.listObjectNames();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "GetListObjectNames HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    void WriteSpMtFromArma(HighFive::File *file, const Rcpp::S4 &s, const std::string &groupName) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not write dataset, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        std::string filePath = GetH5FilePathOfGroupName(groupName);
        try {
            if(file->exist(groupName) == true) {
                std::stringstream ostr;
                ostr << "Existing group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            arma::sp_mat mat = Rcpp::as<arma::sp_mat>(s);
            if(mat.save(filePath) == true) {
                if((file->exist(groupName + "/features") == true) && (file->exist(groupName + "/barcodes") == true)) {
                    return;
                }

                std::string feature_slot;
                if(file->exist(groupName + "/features") == true) {
                    ::Rf_warning("FORMAT_VERSION >= 3");
                    feature_slot = "features/id";
                    if(file->exist(groupName + "/" + feature_slot) == false) {
                        feature_slot = "features/name";
                    }
                } else {
                    ::Rf_warning("FORMAT_VERSION < 3");
                    feature_slot = "genes";
                    if(file->exist(groupName + "/" + feature_slot) == false) {
                        feature_slot = "gene_names";
                    }
                }

                if((file->exist(groupName + "/" + feature_slot) == true) && (file->exist(groupName + "/barcodes") == true)) {
                    return;
                }

                //Get Dimnames data
                Rcpp::List dim_names = s.slot("Dimnames");

                //Write rownames data
                Rcpp::CharacterVector rownames = dim_names[0];
                std::vector<std::string> arrRowNames(rownames.begin(), rownames.end());
                HighFive::DataSet datasetRowNames = file->createDataSet<std::string>(groupName + "/features", HighFive::DataSpace::From(arrRowNames));
                datasetRowNames.write(arrRowNames);

                //Write colnames data
                Rcpp::CharacterVector colnames = dim_names[1];
                std::vector<std::string> arrColNames(colnames.begin(), colnames.end());
                HighFive::DataSet datasetColNames = file->createDataSet<std::string>(groupName + "/barcodes", HighFive::DataSpace::From(arrColNames));
                datasetColNames.write(arrColNames);

                file->flush();
            }
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteSpMtFromArma (S4) HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    void WriteSpMtFromArma(const arma::sp_mat &mat, const std::string &groupName) {
        std::string filePath = GetH5FilePathOfGroupName(groupName);
        try {
            mat.save(filePath);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteSpMtFromArma HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    void WriteSpMtFromS4(HighFive::File *file, const Rcpp::S4 &mat, const std::string &groupName) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not write dataset, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            Rcpp::IntegerVector dims = mat.slot("Dim");
            arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
            arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
            arma::vec x = Rcpp::as<arma::vec>(mat.slot("x"));
            Rcpp::List dim_names = mat.slot("Dimnames");

            // Write group name data
            file->createGroup(groupName);

            // Write DIM data
            std::vector<unsigned int> arrDims(dims.begin(), dims.end());
            HighFive::DataSet datasetDim = file->createDataSet<unsigned int>(groupName + "/shape", HighFive::DataSpace::From(arrDims));
            datasetDim.write(arrDims);

            //Write i data
            std::vector<unsigned int> arrI(i.begin(), i.end());
            HighFive::DataSet datasetI = file->createDataSet<unsigned int>(groupName + "/indices", HighFive::DataSpace::From(arrI));
            datasetI.write(arrI);

            //Write p data
            std::vector<unsigned int> arrP(p.begin(), p.end());
            HighFive::DataSet datasetP = file->createDataSet<unsigned int>(groupName + "/indptr", HighFive::DataSpace::From(arrP));
            datasetP.write(arrP);

            //Write x data
            std::vector<double> arrX(x.begin(), x.end());
            HighFive::DataSet datasetX = file->createDataSet<double>(groupName + "/data", HighFive::DataSpace::From(arrX));
            datasetX.write(arrX);

            //Write rownames data
            Rcpp::CharacterVector rownames = dim_names[0];
            std::vector<std::string> arrRowNames(rownames.begin(), rownames.end());
            HighFive::DataSet datasetRowNames = file->createDataSet<std::string>(groupName + "/features", HighFive::DataSpace::From(arrRowNames));
            datasetRowNames.write(arrRowNames);

            //Write colnames data
            Rcpp::CharacterVector colnames = dim_names[1];
            std::vector<std::string> arrColNames(colnames.begin(), colnames.end());
            HighFive::DataSet datasetColNames = file->createDataSet<std::string>(groupName + "/barcodes", HighFive::DataSpace::From(arrColNames));
            datasetColNames.write(arrColNames);

            file->flush();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteSpMtFromS4 HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }
    }

    bool CheckFileExist(const std::string &file_path)
    {
        std::ifstream infile(file_path);
        return infile.good();
    }

    arma::sp_mat ReadSpMtAsArma(const std::string &groupName) {
        arma::sp_mat mat;
        std::string filePath = GetH5FilePathOfGroupName(groupName);
        if(CheckFileExist(filePath) == true)
        {
            mat.load(filePath);
            return mat;
        }
        return mat;
    }

    Rcpp::S4 ReadSpMtAsS4(HighFive::File *file, const std::string &groupName) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read sparse matrix, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        std::string klass = "dgCMatrix";
        Rcpp::S4 s(klass);

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                Close(file);
                throw;
            }

            std::string feature_slot;
            if(file->exist(groupName + "/features") == true) {
                if((file->exist(groupName + "/features/id") == true) || (file->exist(groupName + "/features/name") == true)) {
                    feature_slot = "features/id";
                    if(file->exist(groupName + "/" + feature_slot) == false) {
                        feature_slot = "features/name";
                    }
                } else {
                    feature_slot = "features";
                }
            } else {
                feature_slot = "genes";
                if(file->exist(groupName + "/" + feature_slot) == false) {
                    feature_slot = "gene_names";
                }
            }

            std::vector<std::string> arrDatasetName = {"data", "indices", "indptr", "shape", feature_slot};
            for(const std::string &datasetName : arrDatasetName) {
                if(file->exist(groupName + "/" + datasetName) == false) {
                    std::stringstream ostr;
                    ostr << "Can not exist dataset : " << datasetName << " in " << groupName;
                    ::Rf_error(ostr.str().c_str());
                    Close(file);
                    throw;
                }
            }

            std::vector<int> arrDims;
            ReadDatasetVector<int>(file, groupName, "shape", arrDims);
            std::vector<int> arrIndices;
            ReadDatasetVector<int>(file, groupName, "indices", arrIndices);
            std::vector<int> arrIndptr;
            ReadDatasetVector<int>(file, groupName, "indptr", arrIndptr);
            std::vector<double> arrData;
            ReadDatasetVector<double>(file, groupName, "data", arrData);
            std::vector<std::string> arrFeature;
            ReadDatasetVector(file, groupName, feature_slot, arrFeature);
            std::vector<std::string> arrBarcode;

            s.slot("p") = std::move(arrIndptr);
            s.slot("i") = std::move(arrIndices);
            s.slot("x") = std::move(arrData);
            s.slot("Dim") = std::move(arrDims);

            if (file->exist(groupName + "/" + "barcodes")) {
              ReadDatasetVector(file, groupName, "barcodes", arrBarcode);
              s.slot("Dimnames") = Rcpp::List::create(arrFeature, arrBarcode);
            } else{
              arrBarcode.resize(arrDims[1], "col");
              s.slot("Dimnames") = Rcpp::List::create(arrFeature, arrBarcode);
            }
            return s;
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadSpMtAsS4 in HDF5 format, error=" << err.what();
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }

        return s;
    }

    Rcpp::List Read10XH5(HighFive::File *file, const std::string &filePath, const bool &use_names, const bool &unique_features) {
        if(file == nullptr) {
            std::stringstream ostr;
            ostr << "Can not read dataset, please open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        Environment pkg = Environment::namespace_env("Matrix");
        Function fMatrix = pkg["sparseMatrix"];

        Rcpp::List arrList = Rcpp::List::create();
        try {
            std::vector<std::string> genomes;
            GetListRootObjectNames(file, genomes);

            for(std::string &groupName : genomes) {
                if(file->exist(groupName) == false) {
                    std::stringstream ostr;
                    ostr << "Can not exist group :" << groupName;
                    ::Rf_error(ostr.str().c_str());
                    Close(file);
                    throw;
                }

                std::string feature_slot;
                if(file->exist(groupName + "/features") == true) {
                    ::Rf_warning("FORMAT_VERSION >= 3");
                    feature_slot = "features/id";
                    if(file->exist(groupName + "/" + feature_slot) == false) {
                        feature_slot = "features/name";
                    }
                } else {
                    ::Rf_warning("FORMAT_VERSION < 3");
                    feature_slot = "genes";
                    if(file->exist(groupName + "/" + feature_slot) == false) {
                        feature_slot = "gene_names";
                    }
                }

                std::vector<std::string> arrDatasetName = {"data", "indices", "indptr", "shape", feature_slot, "barcodes"};
                for(const std::string &datasetName : arrDatasetName) {
                    if(file->exist(groupName + "/" + datasetName) == false) {
                        std::stringstream ostr;
                        ostr << "Can not exist dataset : " << datasetName << " in " << groupName;
                        ::Rf_error(ostr.str().c_str());
                        Close(file);
                        throw;
                    }
                }

                std::vector<int> arrDims;
                ReadDatasetVector<int>(file, groupName, "shape", arrDims);
                std::vector<int> arrIndices;
                ReadDatasetVector<int>(file, groupName, "indices", arrIndices);
                std::vector<int> arrIndptr;
                ReadDatasetVector<int>(file, groupName, "indptr", arrIndptr);
                std::vector<double> arrData;
                ReadDatasetVector<double>(file, groupName, "data", arrData);
                std::vector<std::string> arrFeature;
                ReadDatasetVector(file, groupName, feature_slot, arrFeature);
                std::vector<std::string> arrBarcode;
                ReadDatasetVector(file, groupName, "barcodes", arrBarcode);

                if (unique_features == true) {
                    std::vector<std::string>::iterator it;
                    it = std::unique (arrFeature.begin(), arrFeature.end());
                    arrFeature.resize(std::distance(arrFeature.begin(),it));
                }

                std::vector<std::string> arrFeatureType;
                if(file->exist(groupName + "/features/feature_type") == true) {
                    ReadDatasetVector(file, groupName, "features/feature_type", arrFeatureType);
                }

                std::vector<std::string> arrFeatureGenome;
                if(file->exist(groupName + "/features/genome") == true) {
                    ReadDatasetVector(file, groupName, "features/genome", arrFeatureGenome);
                }

                Rcpp::List arrFeatureGenomeList = Rcpp::List::create();
                if(arrFeatureGenome.size() > 0) {
                    std::vector<std::string>::iterator it;
                    it = std::unique (arrFeatureGenome.begin(), arrFeatureGenome.end());
                    arrFeatureGenome.resize(std::distance(arrFeatureGenome.begin(),it));

                    std::map<std::string, std::vector<unsigned int>> arrFeatureGenomeMap;
                    for(int i = 0; i < arrFeatureGenome.size(); i++) {
                        std::vector<unsigned int> vecFeature;
                        arrFeatureGenomeMap[arrFeatureGenome[i]] = vecFeature;
                    }

                    unsigned int jLoop = 0;
                    for (const std::string & featureName : arrFeature)
                    {
                        std::vector<std::string> arrItem;
                        boost::split(arrItem, featureName, boost::is_any_of(GENOME_SEPARATOR));
                        arrFeatureGenomeMap[arrItem[0]].push_back(jLoop++);
                    }

                    for(int i = 0; i < arrFeatureGenome.size(); i++) {
                        arrFeatureGenomeList[arrFeatureGenome[i]] = std::move(arrFeatureGenomeMap[arrFeatureGenome[i]]);
                    }
                }

                std::transform(std::begin(arrIndices),std::end(arrIndices),std::begin(arrIndices),[](int x){return x+1;});
                S4 mat = fMatrix(
                  Named("i", arrIndices),
                  Named("p", arrIndptr),
                  Named("x", arrData),
                  Named("dims", arrDims),
                  Named("giveCsparse", false),
                  Named("dimnames",
                  Rcpp::List::create(arrFeature, arrBarcode))
                );

                arrList[groupName] = Rcpp::List::create(
                  Named("mat") = mat,
                  Named("feature_type") = arrFeatureType,
                  Named("feature_genome") = arrFeatureGenomeList);
            }
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "Read10XH5 HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            Close(file);
            throw;
        }

        return arrList;
    }

    HighFive::File *Open(const int &mode) {
        HighFive::File *file = nullptr;
        try {
            int new_mode = HighFive::File::OpenOrCreate;
            switch(mode) {
            case -1:
                new_mode = HighFive::File::OpenOrCreate;
                break;
            default:
                new_mode = HighFive::File::ReadOnly;
            break;
            }
            file = new HighFive::File(file_name, new_mode);
        } catch (HighFive::Exception& err) {
            std::cerr << "Can not open HDF5 file, error=" << err.what() << std::endl;
        }
        return file;
    }

    void Close(HighFive::File *file) {
        if(file != nullptr) {
            delete file;
        }
    }

private:
    std::string file_name;

    std::string GetH5FilePathOfGroupName(const std::string &groupName) {
        std::vector<std::string> arrPath;
        boost::split(arrPath, file_name, boost::is_any_of(PATH_SEPARATOR));

        if(arrPath.size() == 0) {
            std::stringstream ostr;
            ostr << "Invalid file path: " << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        std::string fileName = arrPath[(arrPath.size()-1)];
        std::string fileGroupName = groupName + "_" + arrPath[(arrPath.size()-1)];
        std::string fileGroupPath = boost::replace_all_copy(file_name, fileName, fileGroupName);
        return fileGroupPath;
    }
};

} // namespace bioturing
} // namespace com
#endif //HDF5_UTIL
