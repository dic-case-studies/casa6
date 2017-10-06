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

    def to_int(self,value):
        if isinstance(value,numpy.int32):
            return int(value)
        return value

    def to_intvec(self,value):
        if isinstance(value,list) or isinstance(value,numpy.ndarray):
            ## there is the issue of precision loss from 64bit to 32bit here...
            if all([isinstance(v,numpy.int32) or isinstance(v,numpy.int64) for v in value]):
                return [int(v) for v in value]
        return value

    def to_float(self,value):
        if isinstance(value,int):
            return float(value)
        return value

    def to_floatvec(self,value):
        if isinstance(value,list) or isinstance(value,numpy.ndarray):
            return [float(v) for v in value]
        return value

    def expand_path(self,value):
        if os.path.exists(value):
            return value

        if self.ctsys is None:
            sys.exit("configuration error in CasaCoerce.expand_path( )...")

        for i in self.ctsys.getpath( ):
            if os.path.exists( i + os.sep + value ):
                return i + os.sep + value

        return value

coerce = CasaCoerce( )

