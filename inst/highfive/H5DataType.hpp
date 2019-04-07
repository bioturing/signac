/*
 *  Copyright (c), 2017, Adrien Devresse <adrien.devresse@epfl.ch>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 *
 */
#ifndef H5DATATYPE_HPP
#define H5DATATYPE_HPP

#include <H5Tpublic.h>
#include "H5Object.hpp"


namespace HighFive {

struct TypeMapper;

///
/// \brief HDF5 Data Type
///
class DataType : public Object {
  public:
    DataType();

    bool operator==(const DataType& other) const;

    bool operator!=(const DataType& other) const;

  protected:
    friend class Attribute;
    friend class File;
    friend class DataSet;
};

///
/// \brief create an HDF5 DataType from a C++ type
///
///  Support only basic data type
///
template <typename T>
class AtomicType : public DataType {
  public:
    AtomicType();

    AtomicType(size_t str_size, H5T_str_t str_pad, H5T_cset_t str_cset) {
        _hid = H5Tcopy(H5T_C_S1);

        if (H5Tset_size(_hid, str_size) < 0) {
            HDF5ErrMapper::ToException<DataTypeException>(
                "Unable to define datatype size to variable");
        }

        if (H5Tset_strpad(_hid, str_pad) < 0) {
            HDF5ErrMapper::ToException<DataTypeException>(
                "Unable to define datatype pad to variable");
        }

        // define encoding to UTF-8 by default
        H5Tset_cset(_hid, str_cset);
    }

    typedef T basic_type;
};
}

#include "bits/H5DataType_misc.hpp"

#endif // H5DATATYPE_HPP
