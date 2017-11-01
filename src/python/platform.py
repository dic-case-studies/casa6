##
## CASAtools source file (src/python/platform.py)
##
from __future__ import absolute_import
import sys as _sys

pyversion = float(_sys.version_info[0]) + float(_sys.version_info[1]) / 10

if pyversion < 3:
    str_encode = str
    str_decode = str
else:
    def str_encode(s):
        bytes(s,_sys.getdefaultencoding())
    def str_decode(bs):
        bs.decode(_sys.getdefaultencoding(),"strict")
