##
## CASAtools source file (src/python/platform.py)
##
from __future__ import absolute_import
import sys as _sys
import numpy as _numpy

pyversion = float(_sys.version_info[0]) + float(_sys.version_info[1]) / 10

if pyversion < 3:
    numpy_string_coding = "S"
    str_encode = str
    str_decode = str
    def byte_encode(bs,enc):
        if isinstance(bs,bytearray):
            return bs
        else:
            return bytearray(bs)
    def byte_decode(bs,enc):
        return str(bs)
else:
    numpy_string_coding = "U"
    def str_encode(s):
        if isinstance(s,list):
            return [ val if isinstance(val,bytearray) else bytes(str(val),_sys.getdefaultencoding()) for val in s ]
        else:
            return s if isinstance(s,bytearray) else bytes(str(s),_sys.getdefaultencoding())

    def str_decode(bs):
        if isinstance(bs,list):
            return [ val.decode(_sys.getdefaultencoding( ),"strict") if isinstance(val,bytearray) or isinstance(val,bytes) else val for val in bs ]
        else:
            return bs.decode(_sys.getdefaultencoding( ),"strict") if isinstance(bs,bytearray) or isinstance(bs,bytes) else bs

    def byte_encode(bs,enc):
        if isinstance(bs,bytearray):
            return bs
        else:
            return bytearray(bs,enc)
    def byte_decode(bs,enc):
        if isinstance(bs,bytearray):
            return bs.decode(enc,"strict")
        else:
            return str(bs)

str2bytes = str_encode
bytes2str = str_decode

def dict_encode(d):
    coerce =  {
        "string": lambda x: str_encode(str(x)),
        "int": lambda x: int(x),
        "bool": lambda x: bool(x),
        "float": lambda x: float(x),
        "double": lambda x: float(x),
        "floatArray": lambda x: [float(y) for y in x],
        "doubleArray": lambda x: [float(y) for y in x],
        "intArray": lambda x: [int(y) for y in x],
        "boolArray": lambda x: [bool(y) for y in x],
        "byteArray": lambda x: bytearray(x),
        "bytearray": lambda x: bytearray(x),
        "record": lambda x: dict(x),
        "dict": lambda x: dict(x),
    }
    
    if isinstance(d,dict):
        if len(d) == 2 and 'type' in d and 'value' in d and d['type'] in coerce:
            return coerce[d['type']](d['value'])
        else:
            return dict(list(map(encode,list(d.items( )))))
    else:
        return encode(d)

def dict_decode(d):
    coerce =  {
        "string": lambda x: str_decode(x),
        "int": lambda x: int(x),
        "bool": lambda x: bool(x),
        "float": lambda x: float(x),
        "double": lambda x: float(x),
        "floatArray": lambda x: [float(y) for y in x],
        "doubleArray": lambda x: [float(y) for y in x],
        "intArray": lambda x: [int(y) for y in x],
        "boolArray": lambda x: [bool(y) for y in x],
        "byteArray": lambda x: x if isinstance(x,bytearray) else byte_encode(str(x),_sys.getdefaultencoding( )),
        "bytearray": lambda x: x if isinstance(x,bytearray) else byte_encode(str(x),_sys.getdefaultencoding( )),
        "record": lambda x: dict(x),
        "dict": lambda x: dict(x),
    }
    
    if isinstance(d,dict):
        if len(d) == 2 and 'type' in d and 'value' in d and d['type'] in coerce:
            return coerce[d['type']](d['value'])
        else:
            return dict(list(map(decode,list(d.items( )))))
    else:
        return decode(d)

def isiterable(v):
    try:
        _iter = iter(v)
        return True
    except TypeError:
        return False
    
def encode(v):
    if isinstance(v,str):
        return str_encode(v)
    elif isinstance(v,dict):
        return dict_encode(v)
    elif isiterable(v):
        if isinstance(v,list):
            return [encode(x) for x in v]
        elif isinstance(v,tuple):
            return tuple([encode(x) for x in v])
        elif isinstance(v,_numpy.ndarray) and v.dtype.type is _numpy.str_: 
            return _numpy.array([encode(x) for x in v])
        else:
            return v
    else:
        return v

def decode(v):
    if isinstance(v,bytearray) or isinstance(v,bytes):
        return str_decode(v)
    elif isinstance(v,dict):
        return dict_decode(v)
    elif isiterable(v):
        if isinstance(v,list):
            return [decode(x) for x in v]
        elif isinstance(v,tuple):
            return tuple([decode(x) for x in v])
        elif isinstance(v,_numpy.ndarray) and v.dtype.type is _numpy.string_:
            size = max([15]+[len(x) for x in v.ravel( )])
            return _numpy.array([decode(x) for x in v],dtype='|%s%d' % (numpy_string_coding,size+1))
        else:
            return v
    else:
        return v

