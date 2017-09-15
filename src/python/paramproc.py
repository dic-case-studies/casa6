from __future__ import absolute_import 
from CASAtools.__cerberus__ import Validator

class CasaValidator(Validator):
    def _validate_type_cfloat(self,value):                                                                                                                                               
        if isinstance(value,float):
            return True


validator = CasaValidator( )
