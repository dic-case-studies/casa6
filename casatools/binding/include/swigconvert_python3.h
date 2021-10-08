#ifndef _casac_conversions_python3_h__
#define _casac_conversions_python3_h__
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdcasa/record.h>
#include <string>
#include <casaswig_types.h>

namespace casac {

    PyObject *fetch_dict_value( PyObject *dict, const char *key );

}

#endif
