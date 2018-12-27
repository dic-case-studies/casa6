##
## CASAtools source file (src/python/coercetype.py)
##
from __future__ import absolute_import
from .platform import str_encode
import sys
import os
import numpy

class CasaCoerce:

    def __init__(self, *args, **kwargs):
        self.ctsys = None

    def set_ctsys(self, util_obj):
        self.ctsys = util_obj

    def to_list(self,value):
        if isinstance(value,str):
            return value
        try:
            i = iter(value)
            return list(value)
        except TypeError:
            return value

    def to_int(self,value):
        if isinstance(value,numpy.int32):
            return int(value)
        return value

    def to_float(self,value):
        if isinstance(value,int) or isinstance(value,numpy.int32) or isinstance(value,numpy.int64) or \
           isinstance(value,numpy.float64) or isinstance(value,numpy.float32):
            return float(value)
        return value

    def to_intvec(self,value):
        if isinstance(value,int) or isinstance(value,numpy.int32):
            return [int(value)]
        return value

    def to_floatvec(self,value):
        if type(value) in [float,int,numpy.int32,numpy.int64,numpy.float32,numpy.float64]:
            return [float(value)]
        if isinstance(value,list):
            if all( [ isinstance(v,float) or isinstance(v,numpy.float64) or isinstance(v,numpy.float32) or \
                      isinstance(v,int) or isinstance(v,numpy.int32) or isinstance(v,numpy.int64) for v in value]):
                return [float(v) for v in value]
        if isinstance(value,numpy.ndarray) and len(value.shape) <= 1:
            if value.dtype.type in [ numpy.float32, numpy.float64, numpy.int32, numpy.int64, int ]:
                return [float(v) for v in value]
        return value

    def to_strvec(self,value):
        if isinstance(value,str):
            return [value]
        return value

    def to_intarray(self,value):
        if isinstance(value,int) or isinstance(value,numpy.int32):
            return numpy.array([int(value)])
        return value

    def to_floatarray(self,value):
        if type(value) in [float,int,numpy.int32,numpy.int64]:
            return numpy.array([float(value)])
        if isinstance(value,list):
            if all([isinstance(v,int) or isinstance(v,numpy.int32) or isinstance(v,numpy.int64) for v in value]):
                return numpy.array([float(v) for v in value])
        if isinstance(value,numpy.ndarray):
            if value.dtype.type in [ numpy.float32, numpy.float64, numpy.int32, numpy.int64, int ]:
                return value.astype(numpy.float64)
        return value

    def to_strarray(self,value):
        if isinstance(value,str):
            return numpy.array([value])
        if isinstance(value,list):
            if all([isinstance(v,str) for v in value]):
                return numpy.array(value)
        return value

    def to_variant(self,value):
        if isinstance(value,numpy.int32):
            return int(value)
        if isinstance(value,numpy.int64):
            return float(value)
        if isinstance(value,numpy.float32) or isinstance(value,numpy.float64):
            return float(value)
        return value

    def expand_path(self,value):
        if not isinstance(value,str):
            return value
        if len(value) == 0 or value.startswith("/") or value.startswith("./") or value.startswith("../"):
            return value
        if os.path.exists(value):
            return value

        if self.ctsys is None:
            sys.exit("configuration error in CasaCoerce.expand_path( )...")

        ###
        ### cerberus validation is not reentrant...
        ###
        return self.ctsys._swigobj.resolve(str_encode(value))

    def expand_pathvec(self,value):
        if not isinstance(value,list):
            return value
        if any([not isinstance(v,str) for v in value]):
            return value

        if self.ctsys is None:
            sys.exit("configuration error in CasaCoerce.expand_pathvec( )...")

        ###
        ### cerberus validation is not reentrant...
        ###
        new_value = map(lambda v: v if len(v) == 0 or v.startswith("/") or v.startswith("./") or v.startswith("../") else self.ctsys._swigobj.resolve(str_encode(v)), value)
        return list(new_value)

coerce = CasaCoerce( )

