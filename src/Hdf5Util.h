#ifndef HDF5_UTIL
#define HDF5_UTIL

#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include "CommonUtil.h"
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

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
    void WriteDatasetVector(const std::string &groupName, const std::string &datasetName, const std::vector<T> &datasetVec) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == true) {
                std::stringstream ostr;
                ostr << "Existing group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == true) {
                std::stringstream ostr;
                ostr << "Existing dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::DataSet datasetGroup = file->createDataSet<T>(groupName, HighFive::DataSpace::From(datasetVec));
            datasetGroup.write(datasetVec);
            file->flush();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteVector HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    template <typename T>
    void ReadDatasetVector(const std::string &groupName, const std::string &datasetName, std::vector<T> &vvec) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::DataSet datasetVec = file->getDataSet(groupName + "/" + datasetName);
            datasetVec.read(vvec);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadDatasetVector HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    template <typename T>
    void WriteDatasetSingle(const std::string &groupName, const std::string &datasetName, const T &datasetVal) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == true) {
                std::stringstream ostr;
                ostr << "Existing group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == true) {
                std::stringstream ostr;
                ostr << "Exist dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::DataSet datasetGroup = file->createDataSet<T>(groupName, HighFive::DataSpace::From(datasetVal));
            datasetGroup.write(datasetVal);
            file->flush();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteDatasetSingle HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    template <typename T>
    void ReadDatasetSingle(const std::string &groupName, const std::string &datasetName, T &datasetVal) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::DataSet datasetData = file->getDataSet(groupName + "/" + datasetName);
            datasetData.read(datasetVal);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadDatasetSingle HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    template <typename T>
    void WriteRootDataset(const std::string &datasetName, const T &datasetVal) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(datasetName) == true) {
                std::stringstream ostr;
                ostr << "Existing group :" << datasetName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::DataSet datasetGroup = file->createDataSet<T>(datasetName, HighFive::DataSpace::From(datasetVal));
            datasetGroup.write(datasetVal);
            file->flush();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteRootDataset HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    template <typename T>
    void ReadRootDataset(const std::string &datasetName, T &datasetVal) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::DataSet datasetData = file->getDataSet(datasetName);
            datasetData.read(datasetVal);
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadRootDataset HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    void GetListAttributes(const std::string &groupName, const std::string &datasetName, std::vector<std::string> &arrList) {
        boost::shared_ptr<HighFive::File> file = Open(1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            if(file->exist(groupName + "/" + datasetName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist dataset :" << datasetName << "in " << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::DataSet datasetData = file->getDataSet(groupName + "/" + datasetName);
            arrList = datasetData.listAttributeNames();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "GetListAttributes HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    void GetListRootObjectNames(std::vector<std::string> &arrList) {
        boost::shared_ptr<HighFive::File> file = Open(1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            arrList = file->listObjectNames();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "GetListObjectNames HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    bool CheckRootAttribute(const std::string &attName) {
        boost::shared_ptr<HighFive::File> file = Open(1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
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
            throw;
        }
        return isExist;
    }

    void GetListRootAttributes(std::vector<std::string> &arrList) {
        boost::shared_ptr<HighFive::File> file = Open(1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            arrList = file->listAttributeNames();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "GetListRootAttributes HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    void GetListObjectNames(const std::string &groupName, std::vector<std::string> &arrList) {
        boost::shared_ptr<HighFive::File> file = Open(1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == false) {
                std::stringstream ostr;
                ostr << "Can not exist group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            HighFive::Group datasetGroup = file->getGroup(groupName);
            arrList = datasetGroup.listObjectNames();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "GetListObjectNames HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }
    }

    bool WriteSpMtFromArma(const arma::sp_mat &mat, const std::string &groupName) {
        std::string filePath = GetH5FilePathOfGroupName(groupName);
        ::Rf_error(filePath.c_str());
        return mat.save(filePath);
    }

    void WriteSpMtFromS4(const Rcpp::S4 &mat, const std::string &groupName) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        try {
            if(file->exist(groupName) == true) {
                std::stringstream ostr;
                ostr << "Existing group :" << groupName;
                ::Rf_error(ostr.str().c_str());
                throw;
            }

            Rcpp::IntegerVector dims = mat.slot("Dim");
            arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
            arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
            arma::vec x = Rcpp::as<arma::vec>(mat.slot("x"));

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
            HighFive::DataSet datasetX = file->createDataSet<double>(groupName + "/data", HighFive::DataSpace::From(arrP));
            datasetX.write(arrX);

            file->flush();
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "WriteSpMtFromS4 HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
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

    Rcpp::S4 ReadSpMtAsS4(const std::string &groupName){
        boost::shared_ptr<HighFive::File> file = Open(1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
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
                throw;
            }

            std::vector<std::string> arrDatasetName = {"data", "indices", "indptr", "shape"};
            for(const std::string &datasetName : arrDatasetName) {
                if(file->exist(groupName + "/" + datasetName) == false) {
                    std::stringstream ostr;
                    ostr << "Can not exist dataset : " << datasetName << " in " << groupName;
                    ::Rf_error(ostr.str().c_str());
                    throw;
                }
            }

            HighFive::DataSet datasetShape = file->getDataSet(groupName + "/shape");
            HighFive::DataSet datasetIndices = file->getDataSet(groupName + "/indices");
            HighFive::DataSet datasetIndptr = file->getDataSet(groupName + "/indptr");
            HighFive::DataSet datasetData = file->getDataSet(groupName + "/data");

            std::vector<int> arrDims;
            datasetShape.read(arrDims);
            std::vector<int> arrIndices;
            datasetIndices.read(arrIndices);
            std::vector<int> arrIndptr;
            datasetIndptr.read(arrIndptr);
            std::vector<double> arrData;
            datasetData.read(arrData);

            s.slot("p") = std::move(arrIndptr);
            s.slot("i") = std::move(arrIndices);
            s.slot("x") = std::move(arrData);
            s.slot("Dim") = std::move(arrDims);
            return s;
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "ReadSpMtAsS4 in HDF5 format, error=" << err.what();
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        return s;
    }

    Rcpp::List Read10XH5(const std::string &filePath, const bool &use_names, const bool &unique_features, const bool &give_csparse) {
        boost::shared_ptr<HighFive::File> file = Open(1);

        if(file.get() == nullptr) {
            std::stringstream ostr;
            ostr << "Can not open file :" << file_name;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        Rcpp::List arrList = Rcpp::List::create();
        try {
            std::vector<std::string> genomes;
            GetListRootObjectNames(genomes);

            for(std::string &groupName : genomes) {
                if(file->exist(groupName) == false) {
                    std::stringstream ostr;
                    ostr << "Can not exist group :" << groupName;
                    ::Rf_error(ostr.str().c_str());
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
                std::cout << "arrDatasetName size=" << arrDatasetName.size() << std::endl;
                for(const std::string &datasetName : arrDatasetName) {
                    if(file->exist(groupName + "/" + datasetName) == false) {
                        std::stringstream ostr;
                        ostr << "Can not exist dataset : " << datasetName << " in " << groupName;
                        ::Rf_error(ostr.str().c_str());
                        throw;
                    }
                }

                std::cout << "getDataSet shape" << std::endl;
                HighFive::DataSet datasetShape = file->getDataSet(groupName + "/shape");
                HighFive::DataSet datasetIndices = file->getDataSet(groupName + "/indices");
                HighFive::DataSet datasetIndptr = file->getDataSet(groupName + "/indptr");
                HighFive::DataSet datasetData = file->getDataSet(groupName + "/data");
                HighFive::DataSet datasetFeature = file->getDataSet(groupName + "/" + feature_slot);
                HighFive::DataSet datasetBarcode = file->getDataSet(groupName + "/barcodes");

                std::vector<int> arrDims;
                datasetShape.read(arrDims);
                std::vector<int> arrIndices;
                datasetIndices.read(arrIndices);
                std::vector<int> arrIndptr;
                datasetIndptr.read(arrIndptr);
                std::vector<double> arrData;
                datasetData.read(arrData);

                std::vector<std::string> arrFeature;
                //datasetFeature.read(arrFeature);
                std::vector<std::string> arrBarcode;
                //datasetBarcode.read(arrBarcode);

                std::cout << "arrIndptr size=" << arrIndptr.size() << std::endl;
                std::string klass = "dgTMatrix";
                Rcpp::S4 s(klass);
                s.slot("p") = std::move(arrIndptr);
                s.slot("i") = std::move(arrIndices);
                s.slot("x") = std::move(arrData);
                s.slot("Dim") = std::move(arrDims);

                if (unique_features == true) {
                    std::vector<std::string>::iterator it;
                    it = std::unique (arrFeature.begin(), arrFeature.end());
                    arrFeature.resize(std::distance(arrFeature.begin(),it));
                }
                s.slot("Dimnames") = Rcpp::List::create(arrFeature, arrBarcode);

                if(file->exist(groupName + "/features/feature_type") == false) {
                    HighFive::DataSet datasetFeatureType = file->getDataSet(groupName + "/features/feature_type");
                    std::vector<std::string> arrFeatureType;
                    datasetFeatureType.read(arrBarcode);
                    s.slot("factors") = Rcpp::List::create(arrFeatureType);
                }

                arrList[groupName] = s;
            }
        } catch (HighFive::Exception& err) {
            std::stringstream ostr;
            ostr << "Read10XH5 HDF5 format, error=" << err.what() ;
            ::Rf_error(ostr.str().c_str());
            throw;
        }

        return arrList;
    }

private:
    std::string file_name;

    boost::shared_ptr<HighFive::File> Open(const int &mode) {
        boost::shared_ptr<HighFive::File> file;
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
            file.reset(new HighFive::File(file_name, new_mode));
        } catch (HighFive::Exception& err) {
            std::cerr << "Can not open HDF5 file, error=" << err.what() << std::endl;
        }
        return file;
    }

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
