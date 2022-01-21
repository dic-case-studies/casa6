#include <swigconvert_python3.h>

namespace casac {

    PyObject *fetch_dict_value( PyObject *dict, const char *key ) {
        PyObject *result = PyDict_GetItemString( dict, key );
        if ( ! result ) {
            PyObject *pyby_key = PyBytes_FromStringAndSize( key, strlen(key) );
            result = PyDict_GetItem( dict, pyby_key );
            Py_DECREF( pyby_key );
            if ( ! result ) {
                PyObject *pyba_key = PyByteArray_FromStringAndSize( key, strlen(key) );
                result = PyDict_GetItem( dict, pyba_key );
                Py_DECREF( pyba_key );
            }
        }
        return result;
    }

}
