#ifndef HDF5_MODEL
#define HDF5_MODEL

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include "CommonUtil.h"
#include <h5cpp/hdf5.hpp>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace ::hdf5;

using _HDF5_FILE_IMAP = std::map<std::string, file::File>;
using _HDF5_FILE_MUTEX = std::map<std::string, std::mutex *>;

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

    // read functions
    bool readIntegerBlock(const std::string &fileID, const std::string &groupPath,
                          const std::string &dataName, int64_t iOffset,
                          int64_t iLimit, std::vector<int64_t> &arrReturnData);
    bool readFloatBlock(const std::string &fileID, const std::string &groupPath,
                        const std::string &dataName, int64_t iOffset,
                        int64_t iLimit, std::vector<double> &arrReturnData);

    // write functions
    bool writeIntegerBlock(const std::string &fileID,
                           const std::string &groupPath,
                           const std::string &dataName,
                           const std::vector<int64_t> &data);
    bool writeFloatBlock(const std::string &fileID, const std::string &groupPath,
                         const std::string &dataName,
                         const std::vector<double> &data);
    bool writeFloatMatrix(const std::string &fileID, const std::string &groupPath,
                          const std::string &dataName,
                          const arma::sp_mat &matrix);

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
