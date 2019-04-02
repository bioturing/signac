#ifndef HDF5_UTIL
#define HDF5_UTIL

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include "CommonUtil.h"
#include <highfive/H5File.hpp>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

namespace com {
namespace bioturing {

class Hdf5Util {
public:
    Hdf5Util(const std::string &file_name_) {
        file_name = std::string(file_name_.data(), file_name_.size());
    }

    ~Hdf5Util() {}

    template <typename T>
    bool WriteVector(const std::vector<T> &arrVec, const std::string &groupName) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            return false;
        }

        try {
            HighFive::DataSet datasetGroup = file->createDataSet<double>(groupName, HighFive::DataSpace::From(arrVec));
            datasetGroup.write(arrVec);
            file->flush();
            return true;
        } catch (HighFive::Exception& err) {
            std::cerr << "Write Matrix to HDF5 format, error=" << err.what() << std::endl;
        }

        return false;
    }

    bool WriteString(const std::string &stringName, const std::string &groupName) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            return false;
        }

        try {
            HighFive::DataSet datasetGroup = file->createDataSet<std::string>(groupName, HighFive::DataSpace::From(stringName));
            datasetGroup.write(stringName);
            file->flush();
            return true;
        } catch (HighFive::Exception& err) {
            std::cerr << "Write String to HDF5 format, error=" << err.what() << std::endl;
        }

        return false;
    }

    void GetListAttributes(const std::string &groupName, std::vector<std::string> &arrList) {
        boost::shared_ptr<HighFive::File> file = Open(HighFive::File::ReadOnly);

        if(file.get() == nullptr) {
            return;
        }

        try {
            HighFive::DataSet datasetGroup = file->getDataSet(groupName);
            arrList = datasetGroup.listAttributeNames();
        } catch (HighFive::Exception& err) {
            std::cerr << "Get list attributes of group in HDF5 format, error=" << err.what() << std::endl;
        }
    }

    bool WriteSpMtFromArma(const arma::sp_mat &mat, const std::string &groupName) {
        boost::filesystem::path groupPath;
        GetH5FilePathOfGroupName(groupName, groupPath);
        return mat.save(groupPath.c_str());
    }

    bool WriteSpMtFromS4(const Rcpp::S4 &mat, const std::string &groupName) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            return false;
        }

        try {
            Rcpp::IntegerVector dims = mat.slot("Dim");
            arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
            arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
            arma::vec x = Rcpp::as<arma::vec>(mat.slot("x"));

            // Write name data
            HighFive::DataSet datasetGroup = file->createDataSet<std::string>(groupName, HighFive::DataSpace::From(groupName));
            datasetGroup.write(groupName);

            // Write DIM data
            std::vector<unsigned int> arrDims(dims.begin(), dims.end());
            HighFive::Attribute datasetDim = datasetGroup.createAttribute<unsigned int>(
                "shape", HighFive::DataSpace::From(arrDims));
            datasetDim.write(arrDims);

            //Write i data
            std::vector<unsigned int> arrI(i.begin(), i.end());
            HighFive::Attribute datasetI = datasetGroup.createAttribute<unsigned int>(
                "indices", HighFive::DataSpace::From(arrI));
            datasetI.write(arrI);

            //Write p data
            std::vector<unsigned int> arrP(p.begin(), p.end());
            HighFive::Attribute datasetP = datasetGroup.createAttribute<double>(
                "indptr", HighFive::DataSpace::From(arrP));
            datasetP.write(arrP);

            //Write x data
            std::vector<double> arrX(x.begin(), x.end());
            HighFive::Attribute datasetX = datasetGroup.createAttribute<double>(
                "data", HighFive::DataSpace::From(arrX));
            datasetX.write(arrX);

            file->flush();
            return true;
        } catch (HighFive::Exception& err) {
            std::cerr << "WriteSpMtFromS4 in HDF5 format, error=" << err.what() << std::endl;
        }
        return false;
    }

    arma::sp_mat ReadSpMtAsArma(const std::string &groupName) {
        arma::sp_mat mat;
        boost::filesystem::path groupPath;
        GetH5FilePathOfGroupName(groupName, groupPath);
        mat.load(groupPath.c_str());
        return mat;
    }

    arma::sp_mat ReadSpMtAsS4(const std::string &groupName){
        boost::shared_ptr<HighFive::File> file = Open(HighFive::File::ReadOnly);
        arma::sp_mat mat;

        if(file.get() == nullptr) {
            return mat;
        }

        try {
            HighFive::DataSet datasetGroup = file->getDataSet(groupName);

            if(datasetGroup.hasAttribute("shape") == false) {
                std::stringstream ostr;
                ostr << "Can not see shape entry in [" << groupName << "]";
                throw std::range_error(ostr.str());
            }
            std::vector<unsigned int> arrDims;
            datasetGroup.getAttribute("shape").read(arrDims);

            if(datasetGroup.hasAttribute("indices") == false) {
                std::stringstream ostr;
                ostr << "Can not see indices entry in [" << groupName << "]";
                throw std::range_error(ostr.str());
            }
            arma::urowvec idx;
            datasetGroup.getAttribute("indices").read(idx);

            if(datasetGroup.hasAttribute("indptr") == false) {
                std::stringstream ostr;
                ostr << "Can not see indptr entry in [" << groupName << "]";
                throw std::range_error(ostr.str());
            }
            arma::ucolvec ptr;
            datasetGroup.getAttribute("indptr").read(ptr);

            if(datasetGroup.hasAttribute("data") == false) {
                std::stringstream ostr;
                ostr << "Can not see data entry in [" << groupName << "]";
                throw std::range_error(ostr.str());
            }
            arma::colvec x;
            datasetGroup.getAttribute("data").read(x);

            arma::sp_mat res(idx, ptr, x, arrDims[0], arrDims[1]);
            return res;
        } catch (HighFive::Exception& err) {
            std::cerr << "ReadSpMtAsS4 in HDF5 format, error=" << err.what() << std::endl;
        }

        return mat;
    }

private:
    std::string file_name;

    boost::shared_ptr<HighFive::File> Open(const int &mode) {
        boost::shared_ptr<HighFive::File> file;
        try {
            int new_mode = mode;
            if(new_mode == -1) {
                new_mode = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Overwrite;
            }
            file.reset(new HighFive::File(file_name, new_mode));
        } catch (HighFive::Exception& err) {
            std::cerr << "Can not open HDF5 file, error=" << err.what() << std::endl;
        }
        return file;
    }

    void GetH5FilePathOfGroupName(const std::string &groupName, boost::filesystem::path &groupPath) {
        boost::filesystem::path filePath(file_name);
        boost::filesystem::path dirPath(filePath.parent_path());
        std::string fileName = boost::replace_all_copy(groupName, "/", "_");
        fileName = boost::replace_all_copy(fileName, "\\/", "_");
        std::stringstream ss;
        ss << filePath.filename() << "_" << fileName << ".h5";
        boost::filesystem::path fileH5Name(ss.str());
        groupPath = dirPath / fileH5Name;
    }
};

} // namespace bioturing
} // namespace com
#endif //HDF5_UTIL
