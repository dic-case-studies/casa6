##################### generated by xml-casa (v2) from synthesisimstore.xml ##########
##################### ad7f52ee22006d49bc399d2801fb74d2 ##############################
from __future__ import absolute_import 
from .__casac__ import synthesisimstore as _synthesisimstore
from .platform import str_encode as _str_encode
from .platform import str_decode as _str_decode
from .typecheck import validator as _pc
from .coercetype import coerce as _coerce


class synthesisimstore:
    ### self
    def __init__(self, *args, **kwargs):
        """This is used to construct {tt synthesisimstore} tool.
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _synthesisimstore()

    def done(self):
        """
        """
        return self._swigobj.done()

