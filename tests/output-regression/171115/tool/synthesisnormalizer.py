##################### generated by xml-casa (v2) from synthesisnormalizer.xml #######
##################### e7758c1bfca4799c53d8fddadd53d3e4 ##############################
from __future__ import absolute_import 
from .__casac__ import synthesisnormalizer as _synthesisnormalizer
from .platform import str_encode as _str_ec
from .platform import str_decode as _str_dc
from .platform import dict_encode as _dict_ec
from .platform import dict_decode as _dict_dc
from .platform import encode as _any_ec
from .platform import decode as _any_dc
from .typecheck import validator as _pc
from .coercetype import coerce as _coerce
from .synthesisimstore import synthesisimstore as _wrap_synthesisimstore

class synthesisnormalizer:
    ### self
    def __init__(self, *args, **kwargs):
        """This is used to construct {tt synthesisnormalizer} tool.
        """
        self._swigobj = kwargs.get('swig_object',None)
        if self._swigobj is None:
            self._swigobj = _synthesisnormalizer()

    def setupnormalizer(self, normpars={ }):
        """
        """
        schema = {'normpars': {'type': 'cDict'}}
        doc = {'normpars': normpars}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _setupnormalizer_result = self._swigobj.setupnormalizer(_dict_ec(_pc.document['normpars']))
        return _setupnormalizer_result

    def gatherweightdensity(self):
        """
        """
        _gatherweightdensity_result = self._swigobj.gatherweightdensity()
        return _gatherweightdensity_result

    def scatterweightdensity(self):
        """
        """
        _scatterweightdensity_result = self._swigobj.scatterweightdensity()
        return _scatterweightdensity_result

    def gatherpsfweight(self):
        """
        """
        _gatherpsfweight_result = self._swigobj.gatherpsfweight()
        return _gatherpsfweight_result

    def gatherresidual(self):
        """
        """
        _gatherresidual_result = self._swigobj.gatherresidual()
        return _gatherresidual_result

    def dividepsfbyweight(self):
        """
        """
        _dividepsfbyweight_result = self._swigobj.dividepsfbyweight()
        return _dividepsfbyweight_result

    def normalizeprimarybeam(self):
        """
        """
        _normalizeprimarybeam_result = self._swigobj.normalizeprimarybeam()
        return _normalizeprimarybeam_result

    def divideresidualbyweight(self):
        """
        """
        _divideresidualbyweight_result = self._swigobj.divideresidualbyweight()
        return _divideresidualbyweight_result

    def dividemodelbyweight(self):
        """
        """
        _dividemodelbyweight_result = self._swigobj.dividemodelbyweight()
        return _dividemodelbyweight_result

    def multiplymodelbyweight(self):
        """
        """
        _multiplymodelbyweight_result = self._swigobj.multiplymodelbyweight()
        return _multiplymodelbyweight_result

    def scattermodel(self):
        """
        """
        _scattermodel_result = self._swigobj.scattermodel()
        return _scattermodel_result

    def getimstore(self):
        """
        """
        _getimstore_result = _wrap_synthesisimstore(swig_object=self._swigobj.getimstore())
        return _getimstore_result

    def setimstore(self, imstore=None):
        """
        """
        schema = {'imstore': {'type': 'csynthesisimstoreTool'}}
        doc = {'imstore': imstore}
        assert _pc.validate(doc,schema), str(_pc.errors)
        _setimstore_result = self._swigobj.setimstore(_pc.document['imstore'] if _pc.document['imstore'] is None else _pc.document['imstore']._swigobj)
        return _setimstore_result

    def done(self):
        """
        """
        _done_result = self._swigobj.done()
        return _done_result

