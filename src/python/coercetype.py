##
## CASAtools source file (src/python/coercetype.py)
##
from __future__ import absolute_import
import sys
import os
import numpy

class CasaCoerce:

    def __init__(self, *args, **kwargs):
        self.ctsys = None

    def set_ctsys(self, util_obj):
        self.ctsys = util_obj

    def to_strvec(self,value):
        if isinstance(value,str):
            return [value]
        return value

    def to_int(self,value):
        if isinstance(value,numpy.int32):
            return int(value)
        return value

    def to_intvec(self,value):
        if isinstance(value,int) or isinstance(value,numpy.int32):
            return [int(value)]
        return value

    def to_float(self,value):
        if isinstance(value,int):
            return float(value)
        return value

    def to_floatvec(self,value):
        if isinstance(value,float):
            return [value]
        if isinstance(value,list) or isinstance(value,numpy.ndarray):
            return [float(v) for v in value]
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
        return self.ctsys._swigobj.resolve(value)

coerce = CasaCoerce( )

