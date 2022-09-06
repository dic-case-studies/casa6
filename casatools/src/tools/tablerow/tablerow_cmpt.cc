//# tablerow_cmpt.cc
//# Copyright (C) 2022
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#

#include <map>
#include <cstring>
#include <functional>
#include <Python.h>
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <swigconvert_python.h>
#include <table_cmpt.h>

#include <stdint.h>
#include <iostream>
#include <tables/Tables/TableProxy.h>
#include <tables/Tables/TableRowProxy.h>
#include <stdcasa/StdCasa/CasacSupport.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Containers/Record.h>
#include <tablerow_cmpt.h>

using std::function;
using casa::toRecord;
using casa::fromRecord;
using casacore::LogIO;
using casacore::LogOrigin;
using casacore::AipsError;
using casacore::TableRowProxy;
using casacore::Array;
using casacore::Record;
using casacore::Complex;
using casacore::DComplex;
using casacore::String;

static bool _tablerow_initialize_numpy( ) {
    static bool initialized = false;
    if ( initialized == false ) {
        import_array();
        initialized = true;
    }
    return initialized;
}
static bool numpy_initialized = _tablerow_initialize_numpy( );


namespace casac {

    // constructor used by from python to construct a tablerow object
    tablerow::tablerow( const casac::table *_table, const std::vector<std::string> &_columnnames, bool _exclude ) :
        itsLog(new casacore::LogIO)
    {
        if ( ! _table ||
             ! _table->itsTable ||
             ! _table->itsTable->isReadable( ) ) throw AipsError( "invalid table passed for parameter one" );
        itsTable = (table*) _table;
        itsProxy = _table->itsTable;
        itsRow.reset( new TableRowProxy( *itsProxy, static_cast<casacore::Vector<casacore::String>>(_columnnames), _exclude ) );
    }

    // constructor used by table class (in table_cmpt.cc) to return a
    // tablerow for fetching one or more rows
    tablerow::tablerow( table *tb, std::shared_ptr<casacore::TableProxy> myTable,
                        const std::vector<std::string> &columnnames, bool exclude ) :
        itsLog(new casacore::LogIO), itsProxy(myTable), itsTable(tb)
    {
        if ( ! tb ||
             ! tb->itsTable ||
             ! tb->itsTable->isReadable( ) ) throw AipsError( "invalid table passed for parameter one" );
        itsRow.reset( new TableRowProxy( *itsProxy, static_cast<casacore::Vector<casacore::String>>(columnnames), exclude ) );
    }

    // check to see if tablerow can be modified
    bool tablerow::iswritable( ) {
        *itsLog << LogOrigin(__func__,"");
        try {
            if ( itsRow ) return itsRow->isWritable( ) && itsProxy->isWritable( );
        } catch (AipsError x) {
            *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg( ) << LogIO::POST;
            RETHROW(x);
        }
        reset( );
        throw AipsError( "use of uninitialized table row" );
    }
    // magic function which casacore also supplies
    bool tablerow::_iswritable( ) { return iswritable( ); }

    // fetch one row
    record *tablerow::get( long rownr ) {
        *itsLog << LogOrigin(__func__,"");
        try {
            if ( itsRow ) return fromRecord( itsRow->get( rownr ) );
        } catch (AipsError x) {
            *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg( ) << LogIO::POST;
            RETHROW(x);
        }
        reset( );
        throw AipsError( "use of uninitialized table row" );
    }
    // magic function which casacore also supplies
    record *tablerow::_get( long rownr ) { return get(rownr); }

    // replace one row
    bool tablerow::put( long rownr, const record &value, bool matchingfields) {
        *itsLog << LogOrigin(__func__,"");
        try {
            if ( itsRow ) {
                std::unique_ptr<Record> tab(toRecord( value ));
                itsRow->put( rownr, *tab, matchingfields );
                return true;
            }
        } catch (AipsError x) {
            *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg( ) << LogIO::POST;
            RETHROW(x);
        }
        reset( );
        throw AipsError( "use of uninitialized table row" );
    }
    // magic function casacore also supplies
    bool tablerow::_put( long rownr, const record &value, bool matchingfields) { return put( rownr, value, matchingfields ); }

    // magic function for checking the length (it is not completely clear in which
    // contexts this function is used)
    long tablerow::__len__( ) {
        *itsLog << LogOrigin(__func__,"");
        try {
            if ( itsProxy ) return itsProxy->nrows( );
        } catch (AipsError x) {
            *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg( ) << LogIO::POST;
            RETHROW(x);
        }
        reset( );
        throw AipsError( "use of uninitialized table row" );
    }

    // convert a boolean value to a PyObject
    static inline PyObject *toPy( bool b ) {
        if ( b ) { Py_INCREF(Py_True); return Py_True; }
        else { Py_INCREF(Py_False); return Py_False; }
    }

    // convert numeric scalars to a PyObject
#define PY_NUM_SCALAR( CASACORE_TYPE, NUMPY_TYPE )                                 \
    static inline PyObject *toPy( CASACORE_TYPE i ) {                              \
        static PyObject *itemlen = PyLong_FromLong(sizeof(i));                     \
        return PyArray_Scalar( &i, PyArray_DescrFromType(NUMPY_TYPE), itemlen );   \
    }

    PY_NUM_SCALAR( int8_t, NPY_INT8 )
    PY_NUM_SCALAR( uint8_t, NPY_UINT8 )
    PY_NUM_SCALAR( int16_t, NPY_INT16 )
    PY_NUM_SCALAR( uint16_t, NPY_UINT16 )
    PY_NUM_SCALAR( int32_t, NPY_INT32 )
    PY_NUM_SCALAR( uint32_t, NPY_UINT32 )
    PY_NUM_SCALAR( int64_t, NPY_INT64 )
    PY_NUM_SCALAR( uint64_t, NPY_UINT64 )
    PY_NUM_SCALAR( float, NPY_FLOAT )
    PY_NUM_SCALAR( double, NPY_DOUBLE )
    PY_NUM_SCALAR( Complex, NPY_COMPLEX64 )
    PY_NUM_SCALAR( DComplex, NPY_COMPLEX128 )

    // convert a string to a PyObject
    static inline PyObject *toPy( const String &s ) { return PyUnicode_FromString(s.c_str( )); }

    // convert an array of strings to a PyObject
    static inline PyObject *toPy( const Array<String> &a ) {
        auto shape = a.shape( );
        size_t stringlen = std::accumulate( a.begin( ), a.end( ), (size_t) 0, []( size_t tally, const String &s ) { return s.size( ) > tally ? s.length( ) : tally; } );
        size_t memlen = a.nelements( ) * stringlen * sizeof(uint32_t);
        void *mem = PyDataMem_NEW(memlen);
        uint32_t *ptr = reinterpret_cast<uint32_t*>(mem);
        for ( const auto &str : a ) {
            for ( size_t i=0; i < stringlen; ++i ) {
                *ptr++ = i < str.size( ) ? (unsigned char) str[i] : 0;
            }
        }
        return PyArray_New( &PyArray_Type, shape.nelements( ), (npy_intp*) shape.storage( ), NPY_UNICODE, nullptr, mem, stringlen*sizeof(uint32_t), NPY_ARRAY_OWNDATA | NPY_ARRAY_FARRAY, nullptr );
    }

    // convert numeric arrays to PyObjects
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Allocating the result object with:
    //
    //   PyObject *ndarray = PyArray_New( &PyArray_Type, shape.nelements( ), (npy_intp*) shape.storage( ), NUMPY_TYPE,
    //                                    nullptr, nullptr, 0, NPY_ARRAY_FARRAY, nullptr );
    //
    // and then filling it after the fact with:
    //
    //   bool free_storage = false;
    //   auto storage = a.getStorage( free_storage );
    //   std::memcpy( PyArray_DATA( reinterpret_cast<PyArrayObject *>(ndarray)),
    //                              storage, a.nelements( ) * sizeof(CASACORE_TYPE) );
    //   PyArray_ENABLEFLAGS( reinterpret_cast<PyArrayObject *>(ndarray), NPY_ARRAY_OWNDATA );
    //   if ( free_storage ) delete storage;
    //
    // worked on RHEL7 + python 3.6 but fails on macos 10.15 + python 3.8
    //
#define PY_NUM_ARRAY( CASACORE_TYPE, NUMPY_TYPE )                                                                         \
    static inline PyObject *toPy( const Array<CASACORE_TYPE> &a ) {                                                       \
        auto shape = a.shape( );                                                                                          \
        size_t memlen = a.nelements( ) * sizeof(CASACORE_TYPE);                                                           \
        auto *mem = PyDataMem_NEW(memlen);                                                                                \
        auto *ptr = reinterpret_cast<CASACORE_TYPE*>(mem);                                                                \
        for ( const auto &ele : a ) {                                                                                     \
            *ptr++ = ele;                                                                                                 \
        }                                                                                                                 \
        auto result = PyArray_New( &PyArray_Type, shape.nelements( ), (npy_intp*) shape.storage( ), NUMPY_TYPE, nullptr,  \
                                   mem, sizeof(CASACORE_TYPE), NPY_ARRAY_OWNDATA | NPY_ARRAY_FARRAY, nullptr );           \
        /*** setting NPY_ARRAY_OWNDATA here is required to avoid memory leak of allocated data (rhel7 + python 3.6) ***/  \
        PyArray_ENABLEFLAGS( reinterpret_cast<PyArrayObject*>(result), NPY_ARRAY_OWNDATA );                               \
        return result;                                                                                                    \
    }

    PY_NUM_ARRAY( bool, NPY_BOOL )
    PY_NUM_ARRAY( int8_t, NPY_INT8 )
    PY_NUM_ARRAY( uint8_t, NPY_UINT8 )
    PY_NUM_ARRAY( int16_t, NPY_INT16 )
    PY_NUM_ARRAY( uint16_t, NPY_UINT16 )
    PY_NUM_ARRAY( int32_t, NPY_INT32 )
    PY_NUM_ARRAY( uint32_t, NPY_UINT32 )
    PY_NUM_ARRAY( int64_t, NPY_INT64 )
    PY_NUM_ARRAY( uint64_t, NPY_UINT64 )
    PY_NUM_ARRAY( float, NPY_FLOAT )
    PY_NUM_ARRAY( double, NPY_DOUBLE )
    PY_NUM_ARRAY( Complex, NPY_COMPLEX64 )
    PY_NUM_ARRAY( DComplex, NPY_COMPLEX128 )

    static PyObject *toPy( const casacore::Record &rec ) {
        using namespace casacore;

        // build map from table cell types to conversion functions
        std::map<int,function<PyObject*(size_t i)>> function_map = { {TpBool,[&](size_t i) ->PyObject* { return toPy(rec.asBool(i)); }},
                                                                     {TpChar,[&](size_t i) ->PyObject* { return toPy(rec.asuChar(i)); }},
                                                                     {TpUChar,[&](size_t i) ->PyObject* { return toPy(rec.asuChar(i)); }},
                                                                     {TpShort,[&](size_t i) ->PyObject* { return toPy(rec.asShort(i)); }},
                                                                     {TpUShort,[&](size_t i) ->PyObject* { return toPy(rec.asShort(i)); }},
                                                                     {TpInt,[&](size_t i) ->PyObject* { return toPy(rec.asInt(i)); }},
                                                                     {TpUInt,[&](size_t i) ->PyObject* { return toPy(rec.asuInt(i)); }},
                                                                     {TpInt,[&](size_t i) ->PyObject* { return toPy(rec.asInt(i)); }},
                                                                     {TpUInt,[&](size_t i) ->PyObject* { return toPy(rec.asuInt(i)); }},
                                                                     {TpInt64,[&](size_t i) ->PyObject* { return toPy((int64_t)rec.asInt64(i)); }},
                                                                     {TpFloat,[&](size_t i) ->PyObject* { return toPy(rec.asFloat(i)); }},
                                                                     {TpDouble,[&](size_t i) ->PyObject* { return toPy(rec.asDouble(i)); }},
                                                                     {TpComplex,[&](size_t i) ->PyObject* { return toPy(rec.asComplex(i)); }},
                                                                     {TpDComplex,[&](size_t i) ->PyObject* { return toPy(rec.asDComplex(i)); }},
                                                                     {TpArrayBool,[&](size_t i) ->PyObject* { return toPy(rec.asArrayBool(i)); }},
                                                                     {TpArrayUChar,[&](size_t i) ->PyObject* { return toPy(rec.asArrayuChar(i)); }},
                                                                     {TpArrayChar,[&](size_t i) ->PyObject* { return toPy(rec.asArrayuChar(i)); }},
                                                                     {TpArrayShort,[&](size_t i) ->PyObject* { return toPy(rec.asArrayShort(i)); }},
                                                                     {TpArrayUShort,[&](size_t i) ->PyObject* { return toPy(rec.asArrayShort(i)); }},
                                                                     {TpArrayInt,[&](size_t i) ->PyObject* { return toPy(rec.asArrayInt(i)); }},
                                                                     {TpArrayUInt,[&](size_t i) ->PyObject* { return toPy(rec.asArrayuInt(i)); }},
                                                                     {TpArrayFloat,[&](size_t i) ->PyObject* { return toPy(rec.asArrayFloat(i)); }},
                                                                     {TpArrayDouble,[&](size_t i) ->PyObject* { return toPy(rec.asArrayDouble(i)); }},
                                                                     {TpArrayComplex,[&](size_t i) ->PyObject* { return toPy(rec.asArrayComplex(i)); }},
                                                                     {TpArrayDComplex,[&](size_t i) ->PyObject* { return toPy(rec.asArrayDComplex(i)); }},
                                                                     {TpString,[&](size_t i) ->PyObject* { return toPy(rec.asString(i)); }},
                                                                     {TpRecord,[&](size_t i) ->PyObject* { return toPy(rec.asRecord(i)); }}
        };

        // create result
        auto result = PyDict_New( );
        if ( result == nullptr ) throw PyExc_MemoryError;
        // loop through record fields
        for ( uInt i=0; i < rec.nfields( ); ++i ) {
            auto func = function_map.find( rec.dataType(i) );
            // lookup conversion function
            if ( func != function_map.end( ) ) {
                auto newobj = func->second(i);
                auto name = PyUnicode_FromString(rec.name(i).c_str( ));
                // set field in result
                if ( PyDict_SetItem( result, name, newobj ) != 0 ) {
                    Py_DECREF(result);
                    Py_DECREF(newobj);
                    Py_DECREF(name);
                    throw PyExc_ValueError;
                }
                Py_DECREF(newobj);
                Py_DECREF(name);
            } else {
                Py_DECREF(result);
                throw PyExc_TypeError;
            }
        }
        return result;
    }

    PyObj* tablerow::__getitem__( PyObj *rownr ) {
        PyObject *obj = (PyObject*) rownr;
        if ( PyNumber_Check(obj) ) {
            // index indicates a single row
            if ( itsProxy && itsRow ) {
                auto pylong = PyNumber_Long(obj);
                auto index = PyLong_AsLong(pylong);
                Py_DECREF(pylong);
                if ( index >= 0 && index < itsProxy->nrows( ) )
                    return toPy( itsRow->get( index ) );
                else
                    throw PyExc_IndexError;
            }
            throw PyExc_IndexError;
        } else if ( PySlice_Check(obj) ) {

            // index indicates a slice
            if ( itsProxy && itsRow ) {
                Py_ssize_t start, stop, step;
                if ( PySlice_Unpack( obj, &start, &stop, &step ) < 0 ) {
                    throw PyExc_IndexError;
                }
                auto slice_length = PySlice_AdjustIndices( itsProxy->nrows( ), &start, &stop, step );
                auto result = PyList_New( slice_length );
                for ( ssize_t i=0, row=start; i < slice_length; ++i, row += step ) {
                    if ( row < 0 || row >= itsProxy->nrows( ) ) throw PyExc_IndexError;
                    PyObject *newobj = 0;
                    PyGILState_STATE state;  // Needed for PyGILState_Ensure() and PyGILState_Release()
                    state = PyGILState_Ensure( );
                    try {
                        newobj = toPy( itsRow->get( row ) );
                        PyGILState_Release(state);
                    } catch (...) {
                        PyGILState_Release(state);
                    }
                    if ( PyList_SetItem( result, i, newobj ) < 0 ) {
                        Py_DECREF(result);
                        Py_DECREF(newobj);
                        throw PyExc_ValueError;
                    }
                }
                return result;
            }
            throw PyExc_IndexError;
        } else {
            throw PyExc_IndexError;
        }
        return 0;
    }

    tablerow::~tablerow( ) {
        if ( itsTable ) itsTable->remove_tablerow(this);
        itsProxy.reset( );
        itsTable = 0;
    }

    void tablerow::done( ) {
        if ( itsTable ) itsTable->remove_tablerow(this);
        itsRow.reset( );
        itsProxy.reset( );
        itsTable = 0;
    }

    void tablerow::reset( ) {
        itsRow.reset( );
        itsProxy.reset( );
        itsTable = 0;
    }

}
