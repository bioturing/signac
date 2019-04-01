#ifndef HDF5_MODEL
#define HDF5_MODEL

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include "CommonUtil.h"
#include <highfive/H5File.hpp>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

namespace com {
namespace bioturing {

class Hdf5Model {
public:
    Hdf5Model();
    ~Hdf5Model();

    // open and close file
    bool openFile(const std::string &filePath, const file::AccessFlags &iFlag,
                  std::string &fileID);
    std::string createFile(const std::string &fileName, std::string &fileID);
    void closeFile(const std::string &filePath);

private:
    std::mutex mapMutex;
    _HDF5_FILE_IMAP arrFileIMap;
    _HDF5_FILE_MUTEX arrFileMutex;

    bool accessDataset(const file::File &fileInstance,
                       const std::string &groupPath, const std::string &dataName,
                       node::Dataset &dataset);
};

} // namespace bioturing
} // namespace com

#endif //HDF5_MODEL
