from __future__ import absolute_import 
from CASAtools.__cerberus__ import Validator

class CasaValidator(Validator):
    def _validate_type_cInt(self,value):
        if isinstance(value,int):
            return True
    def _validate_type_cFloat(self,value):
        if isinstance(value,float):
            return True
    def _validate_type_cStr(self,value):
        if isinstance(value,str):
            return True
    def _validate_type_cBool(self,value):
        if isinstance(value,bool):
            return True
    def _validate_type_cIntVec(self,value):
        if isinstance(value,list) and \
           all([self._validate_type_cInt(e) for e in value]):
            return True
    def _validate_type_cFloatVec(self,value):
        if isinstance(value,list) and \
           all([self._validate_type_cFloat(e) for e in value]):
            return True
    def _validate_type_cStrVec(self,value):
        if isinstance(value,list) and \
           all([self._validate_type_cStr(e) for e in value]):
            return True
    def _validate_type_cBoolVec(self,value):
        if isinstance(value,list) and \
           all([self._validate_type_cBool(e) for e in value]):
            return True
    def _validate_type_cVariant(self,value):
        if isinstance(value,str) or \
           isinstance(value,int) or \
           isinstance(value,bool) or \
           isinstance(value,dict) or \
           isinstance(value,list) or \
           isinstance(value,float):
            return True


validator = CasaValidator( )
