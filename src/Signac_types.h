#include <vector>
#include <map>
#include <string>
#include <RcppArmadillo.h>
#include <highfive/H5File.hpp>

typedef std::map<std::string,SEXP> LazyFrameT;
typedef std::map<std::string,SEXP>::iterator LazyFrameIteratorT;
