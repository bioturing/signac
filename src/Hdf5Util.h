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
        boost::shared_ptr<HighFive::File> file = Open(1);

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
            return true;
        } catch (HighFive::Exception& err) {
            std::cerr << "WriteSpMtFromS4 in HDF5 format, error=" << err.what() << std::endl;
        }
        return false;
    }

    int WriteVector(const std::vector<double> &vvec, const std::string &groupName, const std::string &datasetName) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            return -1;
        }

        try {
            if(file->exist(groupName) == false) {
                return -2;
            }

            HighFive::DataSet datasetVec = file->createDataSet<double>(groupName + "/" + datasetName, HighFive::DataSpace::From(vvec));
            datasetVec.write(vvec);
            file->flush();

            return 0;
        } catch (HighFive::Exception& err) {
            std::cerr << "WriteVector in HDF5 format, error=" << err.what() << std::endl;
        }
        return 1;
    }

    int ReadVector(std::vector<double> &vvec, const std::string &groupName, const std::string &datasetName) {
        boost::shared_ptr<HighFive::File> file = Open(-1);

        if(file.get() == nullptr) {
            return -1;
        }

        try {
            if(file->exist(groupName) == false) {
                return -2;
            }

            if(file->exist(groupName + "/" + datasetName) == true) {
                HighFive::DataSet datasetVec = file->getDataSet(groupName + "/" + datasetName);
                datasetVec.read(vvec);
            }

            return 0;
        } catch (HighFive::Exception& err) {
            std::cerr << "ReadVector in HDF5 format, error=" << err.what() << std::endl;
        }
        return 1;
    }

    arma::sp_mat ReadSpMtAsArma(const std::string &groupName) {
        arma::sp_mat mat;
        boost::filesystem::path groupPath;
        GetH5FilePathOfGroupName(groupName, groupPath);
        if(boost::filesystem::exists(groupPath.c_str()) == true)
        {
            mat.load(groupPath.c_str());
            return mat;
        }
        return mat;
    }

    Rcpp::S4 ReadSpMtAsS4(const std::string &groupName){
        boost::shared_ptr<HighFive::File> file = Open(1);
        std::string klass = "dgCMatrix";
        Rcpp::S4 s(klass);

        if(file.get() == nullptr) {
            return s;
        }

        try {
            HighFive::DataSet datasetShape = file->getDataSet(groupName + "/shape");
            HighFive::DataSet datasetIndices = file->getDataSet(groupName + "/indices");
            HighFive::DataSet datasetIndptr = file->getDataSet(groupName + "/indptr");
            HighFive::DataSet datasetData = file->getDataSet(groupName + "/data");

            std::vector<int> arrDims;
            datasetShape.read(arrDims);
            std::vector<int> arrD1;
            datasetIndices.read(arrD1);
            std::vector<int> arrD2;
            datasetIndptr.read(arrD2);
            std::vector<double> arrD3;
            datasetData.read(arrD3);

            s.slot("i") = std::move(arrD1);
            s.slot("p") = std::move(arrD2);
            s.slot("x") = std::move(arrD3);
            s.slot("Dim") = std::move(arrDims);
            return s;
        } catch (HighFive::Exception& err) {
            std::cerr << "ReadSpMtAsS4 in HDF5 format, error=" << err.what() << std::endl;
        }

        return s;
    }

    static arma::sp_mat FastConvertS4ToSparseMT(Rcpp::S4 &mat) {
        IntegerVector dims = mat.slot("Dim");
        arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
        arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
        arma::vec x = Rcpp::as<arma::vec>(mat.slot("x"));

        int nrow = dims[0], ncol = dims[1];
        arma::sp_mat res(nrow, ncol);

        arma::access::rw(res.values) = arma::memory::acquire_chunked<double>(x.size() + 1);
        arma::arrayops::copy(arma::access::rwp(res.values), x.begin(), x.size() + 1);

        arma::access::rw(res.row_indices) = arma::memory::acquire_chunked<arma::uword>(i.size() + 1);
        arma::arrayops::copy(arma::access::rwp(res.row_indices), i.begin(), i.size() + 1);

        arma::access::rw(res.col_ptrs) = arma::memory::acquire<arma::uword>(p.size() + 2);
        arma::arrayops::copy(arma::access::rwp(res.col_ptrs), p.begin(), p.size() + 1);

        arma::access::rwp(res.col_ptrs)[p.size()+1] = std::numeric_limits<arma::uword>::max();

        arma::access::rw(res.n_nonzero) = x.size();
        return res;
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

    void GetH5FilePathOfGroupName(const std::string &groupName, boost::filesystem::path &groupPath) {
        boost::filesystem::path filePath(file_name);
        boost::filesystem::path dirPath(filePath.parent_path());
        std::string fileName = boost::replace_all_copy(groupName, "/", "_");
        fileName = boost::replace_all_copy(fileName, "\\/", "_");
        std::stringstream ss;
        ss << filePath.stem().string() << "_" << fileName << ".h5";
        boost::filesystem::path fileH5Name(ss.str());
        groupPath = dirPath / fileH5Name;
    }
};

} // namespace bioturing
} // namespace com
#endif //HDF5_UTIL
